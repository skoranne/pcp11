################################################################################
# File    : KaporinDefence.jl
# Author  : Sandeep Koranne and Advay Koranne (C) 2026. All rights reserved.
# Purpose : Experiment for Low rank adaptation using Kaporin condition number
################################################################################

using LinearAlgebra
using Arpack  # High-performance partial SVD
using Random
using Statistics
using Zygote                 # ⟶ automatic differentiation
using Plots                  # (install with `import Pkg; Pkg.add("Plots")`)

"""
    initialize_pissa(W, r)

Executes a rank-r partial SVD on the base weight matrix W to generate 
the initial orthogonalized low-rank factors A and B.
"""
function initialize_pissa(W::AbstractMatrix, r::Int)
    svd_result, _ = svds(W, nsv=r)
    U, S, V = svd_result.U, svd_result.S, svd_result.V
    
    sqrt_S = Diagonal(sqrt.(S))
    B = U * sqrt_S
    A = sqrt_S * V' # Convert right singular vectors to input rows
    
    return A, B, S
end

"""
    lanczos_tridiag(A_matrix, v, m)

Generates an m x m symmetric tridiagonal Krylov projection matrix 
from the symmetric matrix A_matrix and starting vector v.
"""
function lanczos_tridiag(A_matrix, v::Vector{Float64}, m::Int)
    n = length(v)
    α = zeros(Float64, m)
    β = zeros(Float64, m-1)
    
    q_prev = zeros(Float64, n)
    q_curr = v ./ norm(v)
    w = similar(q_curr)
    
    β_curr = 0.0
    
    for j in 1:m
        mul!(w, A_matrix, q_curr)
        if j > 1
            w .-= β_curr .* q_prev
        end
        α[j] = dot(q_curr, w)
        w .-= α[j] .* q_curr
        
        β_curr = norm(w)
        if j < m
            β[j] = β_curr
            if β_curr < 1e-12
                return SymTridiagonal(α[1:j], β[1:j-1])
            end
            q_prev .= q_curr
            q_curr .= w ./ β_curr
        end
    end
    return SymTridiagonal(α, β)
end

"""
    randomized_kaporin(G, k, m)

Computes the matrix-free Randomized Kaporin estimator utilizing 
Hutchinson traces and Stochastic Lanczos Quadrature.
"""
function randomized_kaporin(G, k::Int=10, m::Int=30)
    n = size(G, 1)
    trace_est = 0.0
    logdet_est = 0.0
    
    for _ in 1:k
        v = rand([-1.0, 1.0], n)
        
        # 1. Hutchinson Trace Component
        trace_est += dot(v, G * v)
        
        # 2. Stochastic Lanczos Quadrature Component
        T_m = lanczos_tridiag(G, v, m)
        vals, vecs = eigen(T_m)
        
        safe_vals = log.(max.(vals, 1e-12))
        τ = (vecs[1, :] .^ 2)' * safe_vals
        logdet_est += (norm(v)^2) * τ
    end
    
    trace_est /= k
    logdet_est /= k
    
    return log(trace_est / n) - (logdet_est / n)
end

function exact_kaporin(G::AbstractMatrix)
    r = size(G, 1)
    comp_tr = tr(G)
    comp_logdet = logdet(G)                 # uses a Cholesky internally (stable for SPD)
    return log(comp_tr / r) - comp_logdet / r
end


# ====================================================================
# Execution Context: Simulating a Qwen2.5 Math MLP Up-Project Block
# ====================================================================
d_out, d_in = 18944, 3584
r = 16

println("Instantiating Base Matrix Layer...")
W_base = randn(Float64, d_out, d_in) / sqrt(d_in)

println("Running PiSSA Factorization Engine...")
A, B, S_vals = initialize_pissa(W_base, r)

# Apply Tikhonov regularization shift to insulate the log operator
ε = 1e-4
G_A = Symmetric(A * A' + ε * I)  # Evaluates input row-rank
G_B = Symmetric(B' * B + ε * I)  # Evaluates output column-rank

println("Evaluating Subspace Structural Multipliers...")
kaporin_A = randomized_kaporin(G_A, 10, min(30, r))
kaporin_B = randomized_kaporin(G_B, 10, min(30, r))

println("\n" * "="^30 * "\n  METRIC EXECUTION REPORT\n" * "="^30)
println("Randomized Kaporin (Matrix A Row-Rank Proxy):    ", round(kaporin_A, digits=4))
println("Randomized Kaporin (Matrix B Column-Rank Proxy): ", round(kaporin_B, digits=4))
println("Total Joint Subspace Regularization Value:       ", round(kaporin_A + kaporin_B, digits=4))

eigvals = eigen(G_A).values           # sorted ascending
println("Eigenvalues (first 5 / last 5): ",
        round.(eigvals[1:5]; digits=4), "...",
        round.(eigvals[end-4:end]; digits=4))
println("Mean eigenvalue = ", mean(eigvals))
println("Geometric mean  = ", exp(mean(log.(eigvals))))
println("Condition (tr/geom) = ", mean(eigvals) / exp(mean(log.(eigvals))))

println("‖A‖_F = ", norm(A))
println("‖A‖_2 = ", opnorm(A))

kaporin_exact = exact_kaporin(G_A)
println("Exact Kaporin (no Monte‑Carlo) = ", round(kaporin_exact; digits=6))
"""
Instantiating Base Matrix Layer...
Running PiSSA Factorization Engine...
Evaluating Subspace Structural Multipliers...

==============================
  METRIC EXECUTION REPORT
==============================
Randomized Kaporin (Matrix A Row-Rank Proxy):    0.0
Randomized Kaporin (Matrix B Column-Rank Proxy): 0.0
Total Joint Subspace Regularization Value:       0.0
Eigenvalues (first 5 / last 5): [3.2566, 3.2584, 3.2626, 3.2636, 3.2659]...[3.2817, 3.2839, 3.2871, 3.2906, 3.2948]                              
Mean eigenvalue = 3.2738964633350687
Geometric mean  = 3.273877881717938
Condition (tr/geom) = 1.0000056757209042
‖A‖_F = 7.237454208032069
‖A‖_2 = 1.815125390252534
Exact Kaporin (no Monte‑Carlo) = 6.0e-6
"""


# ------------------------------------------------------------
# 1.  Keep adapters in BF16 (or FP8) – only the *Gram* is FP32
# ------------------------------------------------------------
A_fp = convert.(BFloat16, A)          # or Float8E4M3 if supported
B_fp = convert.(BFloat16, B)

# ------------------------------------------------------------
# 2.  Build Gram matrices on‑the‑fly (GPU‑resident)
# ------------------------------------------------------------
function gram_A(A_fp)
    # A_fp is (r, d_in); we want AAᵀ (r×r) in FP32
    G = Float32.(A_fp) * Float32.(A_fp)'   # CUDA kernels will handle it
    ε = 1e-4f0 * I                         # same dtype as G
    return Symmetric(G + ε)
end

function gram_B(B_fp)
    G = Float32.(B_fp)' * Float32.(B_fp)   # (r×r)
    ε = 1e-4f0 * I
    return Symmetric(G + ε)
end

# ------------------------------------------------------------
# 3.  Exact Kaporin (GPU Cholesky) – cheap for r ≤ 64
# ------------------------------------------------------------
#function kaporin_exact_gpu(G::CuArray{Float32,2})
function kaporin_exact_gpu(G::Array{Float32,2})
    # Cholesky on the device (cuSOLVER) – O(r³) is trivial
    L = cholesky(G).L               # L is lower‑triangular
    trG   = tr(G)
    logdet = 2.0f0 * sum(log.(diag(L)))   # log‑det from Cholesky
    r = size(G, 1)
    return log(trG / r) - logdet / r
end

# ------------------------------------------------------------
# 4.  Integration in the loss
# ------------------------------------------------------------
λ = 0.05f0                     # regularisation weight (tune!)
function my_loss_fn(ŷ, y, A_fp, B_fp)
    ce = crossentropy(ŷ, y)               # task loss
    kA = kaporin_exact_gpu(gram_A(A_fp))
    kB = kaporin_exact_gpu(gram_B(B_fp))
    return ce + λ * (kA + kB)
end

# --------------------------------------------------------------
# 1️⃣ Kaporin proxy + custom adjoint (unchanged – see previous answer)
# --------------------------------------------------------------
function kaporin_proxy(M::AbstractMatrix{T}, ε::Real) where {T<:Real}
    r = size(M, 1)
    G = Symmetric(M * M' + ε * I)
    L = cholesky(G).L
    trG   = tr(G)
    logdetG = 2.0 * sum(log.(diag(L)))
    return log(trG / r) - logdetG / r
end

function kaporin_grad(M::AbstractMatrix{T}, ε::Real) where {T}
    r = size(M, 1)
    G = Symmetric(M * M' + ε * I)
    invG = inv(G)
    trG = tr(G)
    return (2.0 / r) * (M / trG .- invG * M)
end

Zygote.@adjoint function kaporin_proxy(M::AbstractMatrix{T}, ε::Real) where {T}
    y = kaporin_proxy(M, ε)
    function back(Δ)
        return (Δ * kaporin_grad(M, ε), zero(ε))
    end
    return y, back
end

# --------------------------------------------------------------
# 2️⃣ Helper: effective rank (entropy‑based)
# --------------------------------------------------------------
function effective_rank(M::AbstractMatrix)
    s = svdvals(M) .+ eps()
    p = s / sum(s)
    return exp(-sum(p .* log.(p)))
end

# --------------------------------------------------------------
# 3️⃣ Core simulation (gradient call fixed)
# --------------------------------------------------------------
function simulate_rka(; d_out=128, d_in=64, r=8, batch=32,
                     noise_std=0.01, opt=:adam, η=1e-2,
                     β1=0.9, β2=0.999, λ=0.0, ε=1e-4,
                     steps=500, seed=42, record_every=5)

    rng = MersenneTwister(seed)

    # base weight (fixed throughout training)
    W0 = randn(rng, d_out, d_in) / sqrt(d_in)

    # a single mini‑batch we reuse (synthetic regression)
    X = randn(rng, batch, d_in)
    Y = X * W0' .+ noise_std * randn(rng, batch, d_out)

    # low‑rank factors (tiny random init)
    A = randn(rng, r, d_in) * 0.01
    B = randn(rng, d_out, r) * 0.01

    # Adam state (kept even if we later switch to SGD)
    mA = zeros(size(A)); vA = zeros(size(A))
    mB = zeros(size(B)); vB = zeros(size(B))

    # containers for diagnostics
    nrec = div(steps, record_every) + 1
    loss_hist = Vector{Float64}(undef, nrec)
    kapA_hist = Vector{Float64}(undef, nrec)
    kapB_hist = Vector{Float64}(undef, nrec)
    erA_hist  = Vector{Float64}(undef, nrec)
    erB_hist  = Vector{Float64}(undef, nrec)

    # ----------------------------------------------------------
    # loss + auxiliary metrics (closure for Zygote)
    # ----------------------------------------------------------
    function loss_and_metrics(A_, B_)
        ΔW = B_ * A_                     # low‑rank update
        Ŷ  = X * (W0 .+ ΔW)'             # predictions (batch × d_out)
        mse = sum(abs2, Ŷ .- Y) / batch
        kapA = kaporin_proxy(A_, ε)      # row‑wise Gram
        kapB = kaporin_proxy(B_', ε)     # column‑wise Gram (B' is r×d_out)
        total = mse - λ * (kapA + kapB)
        return total, (mse, kapA, kapB)
    end

    # scalar version used for the gradient call
    loss_fn(A_, B_) = loss_and_metrics(A_, B_)[1]

    # ----------------------------------------------------------
    # optimisation loop
    # ----------------------------------------------------------
    rec = 1
    for t in 1:steps
        # 1️⃣ forward + gradient (the only line that changed)
        ∂A, ∂B = Zygote.gradient(loss_fn, A, B)

        # 2️⃣ Adam update (or simple SGD)
        if opt == :adam
            mA .= β1 .* mA .+ (1-β1) .* ∂A
            vA .= β2 .* vA .+ (1-β2) .* (∂A .^ 2)
            m̂A = mA ./ (1 - β1^t)
            v̂A = vA ./ (1 - β2^t)
            A .-= η .* m̂A ./ (sqrt.(v̂A) .+ 1e-8)

            mB .= β1 .* mB .+ (1-β1) .* ∂B
            vB .= β2 .* vB .+ (1-β2) .* (∂B .^ 2)
            m̂B = mB ./ (1 - β1^t)
            v̂B = vB ./ (1 - β2^t)
            B .-= η .* m̂B ./ (sqrt.(v̂B) .+ 1e-8)
        else                # vanilla SGD
            A .-= η .* ∂A
            B .-= η .* ∂B
        end

        # 3️⃣ record diagnostics every `record_every` steps
        if t % record_every == 0 || t == 1
            ℓ, (mse, kapA, kapB) = loss_and_metrics(A, B)
            loss_hist[rec] = ℓ
            kapA_hist[rec] = kapA
            kapB_hist[rec] = kapB
            erA_hist[rec]  = effective_rank(A)
            erB_hist[rec]  = effective_rank(B')
            rec += 1
        end
    end

    # ----------------------------------------------------------
    # pack results ------------------------------------------------
    # ----------------------------------------------------------
    return Dict(
        :loss   => loss_hist,
        :kapA   => kapA_hist,
        :kapB   => kapB_hist,
        :erA    => erA_hist,
        :erB    => erB_hist,
        :A_final=> copy(A),
        :B_final=> copy(B),
        :settings=> (d_out=d_out, d_in=d_in, r=r,
                    λ=λ, ε=ε, η=η, opt=opt, steps=steps)
    )
end

# --------------------------------------------------------------
# 4️⃣ Demo – compare “no regulariser” vs “with Kaporin”
# --------------------------------------------------------------
@inbounds res_no  = simulate_rka(d_out = 18944, d_in=3584,r=16, λ=0.0, opt=:adam, steps=800, record_every=10)
@inbounds res_reg = simulate_rka(d_out = 18944, d_in=3584,r=16, λ=0.1, opt=:adam, steps=800, record_every=10)

# --------------------------------------------------------------
# 5️⃣ Simple visualisation (requires Plots.jl)
# --------------------------------------------------------------
function plot_trajectories(no, yes; title_suffix="")
    p1 = plot(no[:loss], label="no‑reg", lw=2,
              xlabel="record step", ylabel="total loss",
              title="Training loss $title_suffix")
    plot!(p1, yes[:loss], label="Kaporin λ=0.1", lw=2, ls=:dash)

    p2 = plot(no[:kapA], label="kapA (no‑reg)", lw=2,
              xlabel="record step", ylabel="Kaporin proxy",
              title="Kaporin A‑side $title_suffix")
    plot!(p2, yes[:kapA], label="kapA (λ=0.1)", lw=2, ls=:dash)

    p3 = plot(no[:kapB], label="kapB (no‑reg)", lw=2,
              xlabel="record step", ylabel="Kaporin proxy",
              title="Kaporin B‑side $title_suffix")
    plot!(p3, yes[:kapB], label="kapB (λ=0.1)", lw=2, ls=:dash)

    p4 = plot(no[:erA], label="eff‑rank A (no‑reg)", lw=2,
              xlabel="record step", ylabel="effective rank",
              title="Effective rank of A $title_suffix")
    plot!(p4, yes[:erA], label="eff‑rank A (λ=0.1)", lw=2, ls=:dash)

    p5 = plot(no[:erB], label="eff‑rank B (no‑reg)", lw=2,
              xlabel="record step", ylabel="effective rank",
              title="Effective rank of B $title_suffix")
    plot!(p5, yes[:erB], label="eff‑rank B (λ=0.1)", lw=2, ls=:dash)

    plot(p1, p2, p3, p4, p5; layout=(3,2), size=(960,800),
         legend=:bottomright)
end

plot_trajectories(res_no, res_reg, title_suffix="(r=8, Adam)")
