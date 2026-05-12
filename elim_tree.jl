# ------------------------------------------------------------
# 1. Elimination tree from the sparsity pattern of a symmetric matrix
# ------------------------------------------------------------
"""
    elim_tree(A::SparseMatrixCSC{Tv,Ti}) -> Vector{Int}

Return the elimination tree (parent array) of a symmetric positive
definite matrix `A`.  `A` must be square and stored in CSC format.
The algorithm works on the *structure* only – the numeric values are ignored.
"""
function elim_tree(A::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti<:Integer}
    n = size(A, 1)
    @assert n == size(A, 2) "Matrix must be square"

    parent = zeros(Int, n)          # 0 ⇒ root
    ancestor = zeros(Int, n)        # temporary workspace

    # Walk columns left‑to‑right
    for j = 1:n
        # Scan the lower‑triangular part of column j
        for p = A.colptr[j]:(A.colptr[j+1]-1)
            i = A.rowval[p]          # row index
            i == j && continue       # skip diagonal (i == j)
            i > j && continue        # we only need i < j (strictly lower)
            # Follow the ancestor chain and set parent links
            k = i
            while k != 0 && k < j
                if ancestor[k] == 0
                    ancestor[k] = j
                    parent[k] = j
                    k = 0
                else
                    t = ancestor[k]
                    ancestor[k] = j
                    k = t
                end
            end
        end
        # reset ancestor entries for next column
        for p = A.colptr[j]:(A.colptr[j+1]-1)
            i = A.rowval[p]
            i == j && continue
            i > j && continue
            ancestor[i] = 0
        end
    end
    return parent
end

# ------------------------------------------------------------
# 2. Right‑looking sparse Cholesky using the elimination tree
# ------------------------------------------------------------
"""
    right_looking_cholesky(A::SparseMatrixCSC{Tv,Ti}) -> SparseMatrixCSC{Tv,Ti}

Compute the Cholesky factor `L` (lower triangular) of a symmetric
positive‑definite sparse matrix `A` using a right‑looking algorithm
that explicitly exploits the elimination tree.

The result satisfies `A ≈ L * L'` (up to rounding error).
"""
function right_looking_cholesky(A::SparseMatrixCSC{Tv,Ti}) where {Tv<:Real,Ti<:Integer}
    n = size(A, 1)
    @assert n == size(A, 2) "Matrix must be square"
    @assert issymmetric(A) "Matrix must be symmetric (structure only is checked)"

    # --------------------------------------------------------
    # 2.1 Build elimination tree (structure only)
    # --------------------------------------------------------
    parent = elim_tree(A)

    # --------------------------------------------------------
    # 2.2 Allocate workspace for the factor L
    # --------------------------------------------------------
    # We do not know the exact number of nonzeros a priori,
    # but a safe upper bound is `nnz(A) + n` (diagonal) + fill‑in.
    # We'll grow the arrays on the fly.
    Lcolptr = Vector{Ti}(undef, n + 1)
    Lcolptr[1] = 1
    Lrowval = Ti[]
    Lvals   = Tv[]

    # temporary dense vector that holds the current column of L
    x  = zeros(Tv, n)
    # integer stack that records the pattern of x (indices where x ≠ 0)
    xi = Ti[]

    # For each row i we keep the index of the *most recent* column
    # where i appears in L (used for the column‑merge step)
    prev = zeros(Int, n)

    # --------------------------------------------------------
    # 2.3 Main loop: process columns left‑to‑right
    # --------------------------------------------------------
    for k = 1:n
        # -----------------------------------------------------------------
        # 2.3.1 Gather the pattern of column k of A (strictly lower part)
        # -----------------------------------------------------------------
        empty!(xi)                         # clear pattern stack
        # copy the nonzeros of column k of A into x
        for p = A.colptr[k]:(A.colptr[k+1]-1)
            i = A.rowval[p]
            if i < k                         # lower triangular part only
                x[i] = A.nzval[p]
                push!(xi, i)                # remember that i is non‑zero in x
            elseif i == k
                # diagonal entry will be handled later
                x[i] = A.nzval[p]           # keep it for the sqrt step
                # No need to push i into xi – diagonal is not part of the pattern
            end
        end

        # -----------------------------------------------------------------
        # 2.3.2 Walk up the elimination tree to incorporate already computed
        #        columns (the “right‑looking” updates)
        # -----------------------------------------------------------------
        # The pattern of the current column after all updates is exactly the
        # union of the pattern of A(:,k) and the patterns of its ancestors.
        # The classic column‑merge algorithm does this efficiently by
        # following the `prev` links.
        for idx = 1:length(xi)
            i = xi[idx]                     # current row index in pattern
            # Walk up the tree as long as we have a previous column for i
            # that is lower than the current column k.
            while prev[i] != 0 && prev[i] < k
                j = prev[i]                 # column where row i last appeared
                # The value L[i,j] is already stored in L (see later)
                # Use it to eliminate x[i] (i.e., compute the contribution
                # to the trailing submatrix).
                lij = Lvals[ findfirst( (r)->r==i, Lrowval[Lcolptr[j]:(Lcolptr[j+1]-1)] ) ]
                # Actually we can obtain lij directly because we stored it
                # in the dense vector `x` when we processed column j.
                # For clarity we recompute it here:
                lij = x[i] / Lvals[ Lcolptr[j] ]   # L[j,j] is the first entry of column j
                # Update all rows i' that are non‑zero in column j below row i
                for p = Lcolptr[j]+1:(Lcolptr[j+1]-1)
                    i2 = Lrowval[p]
                    x[i2] -= lij * Lvals[p]
                    # If this is the first time we touch i2 in this column,
                    # add it to the pattern stack.
                    if x[i2] != 0 && !(i2 in xi)
                        push!(xi, i2)
                    end
                end
                # After using column j we can discard it from the pattern.
                prev[i] = parent[j]          # jump to the next ancestor
            end
        end

        # -----------------------------------------------------------------
        # 2.3.3 Compute the diagonal entry L[k,k] = sqrt( x[k] )
        # -----------------------------------------------------------------
        d = x[k]
        @assert d > 0 "Matrix is not positive definite (pivot $k = $d ≤ 0)"
        lkk = sqrt(d)
        # store diagonal at the beginning of column k in L
        push!(Lrowval, k)
        push!(Lvals,   lkk)

        # -----------------------------------------------------------------
        # 2.3.4 Store the off‑diagonal entries of column k (rows i < k)
        # -----------------------------------------------------------------
        # The pattern `xi` now contains exactly the rows i < k where L[i,k] ≠ 0.
        # Sort them (required for CSC format) and write the values.
        sort!(xi)          # ascending order
        for i in xi
            # after the updates, the remaining value in x[i] is L[i,k]
            lik = x[i] / lkk
            push!(Lrowval, i)
            push!(Lvals,   lik)

            # Update the `prev` link: the most recent column where row i appears
            prev[i] = k
        end

        # -----------------------------------------------------------------
        # 2.3.5 Finish column k: set column pointer, clean temporary data
        # -----------------------------------------------------------------
        Lcolptr[k+1] = length(Lrowval) + 1
        # Reset the dense workspace for the next column
        for i in xi
            x[i] = 0.0
        end
        x[k] = 0.0            # diagonal already consumed
        empty!(xi)           # pattern stack cleared
    end

    # --------------------------------------------------------
    # 2.4 Build the SparseMatrixCSC for L and return it
    # --------------------------------------------------------
    return SparseMatrixCSC(n, n, Lcolptr, Lrowval, Lvals)
end

# ------------------------------------------------------------
# 3. Simple test harness
# ------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    using LinearAlgebra, Random, Printf

    # Helper to build a random sparse SPD matrix
    function rand_spd_sparse(n::Int; density::Float64 = 0.05, rng = Random.GLOBAL_RNG)
        # generate a random sparse lower‑triangular matrix R
        R = sprand(rng, n, n, density)
        R = tril(R)                     # keep strictly lower part
        for i = 1:n
            R[i,i] = 0.0                # we will add diagonal later
        end
        # make it diagonally dominant → SPD
        d = [sum(abs.(R[i, :])) + 1.0 for i = 1:n]
        R = R + spdiagm(0 => d)
        # Form A = R * Rᵀ (guaranteed SPD)
        return R * transpose(R)
    end

    # --------------------------------------------------------
    # 3.1 Small sanity check (compare with built‑in cholesky)
    # --------------------------------------------------------
    n = 10
    A = rand_spd_sparse(n, density=0.2)
    L_my = right_looking_cholesky(A)

    # built‑in (uses SuiteSparse)
    L_ref = cholesky(A, check = false).L

    @printf "‖L_my - L_ref‖_F = %.3e\n" norm(L_my - L_ref, Frobenius)

    # --------------------------------------------------------
    # 3.2 Larger test (n = 2000, density = 0.001)
    # --------------------------------------------------------
    n = 2000
    A = rand_spd_sparse(n, density=0.001)
    @time L_my = right_looking_cholesky(A)   # our implementation
    @time L_ref = cholesky(A, check = false).L

    @printf "Relative error (Frobenius) = %.3e\n" norm(L_my - L_ref, Frobenius) / norm(L_ref, Frobenius)
end
