################################################################################
# File     : dp2mat.jl
# Autho    : Sandeep Koranne (C) 2026. All rights reserved
# Purpose  : Dump matrix from GetDP for a EleSta_v problem
# getdp cap.pro -msh SDFFC.msh -sol EleSta_v -pos Map -v 0 -mat_view ascii -vec_view ascii > SDFFC.log
# row 124263: (36163, -943.335)  (36164, -4465.48)
# 
################################################################################

using LinearAlgebra

function convert_to_mtx(input_filename::String, output_filename::String)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    
    # Regex to capture row numbers and (col, val) coordinate pairs
    row_regex = r"row\s+(\d+):"
    pair_regex = r"\(\s*(\d+)\s*,\s*([^\)]+)\)"

    open(input_filename, "r") do file
        for line in eachline(file)
            if isempty(strip(line)) continue end
            
            # 1. Extract the current row index
            row_match = match(row_regex, line)
            if row_match === nothing
                continue
            end
            current_row = 1+parse(Int, row_match.captures[1])
            
            # 2. Extract all (column, value) tuples for this row
            for m in eachmatch(pair_regex, line)
                col = 1+parse(Int, m.captures[1])
                val = parse(Float64, m.captures[2])
                
                push!(rows, current_row)
                push!(cols, col)
                push!(vals, val)
            end
        end
    end

    # Determine structural metadata
    num_entries = length(vals)
    max_row = maximum(rows)
    max_col = maximum(cols)

    # 3. Write Matrix Market File
    open(output_filename, "w") do out
        # Write standard MTX Banner Line
        println(out, "%%MatrixMarket matrix coordinate real general")
        
        # Write Size Line: [Total Rows] [Total Columns] [Total Non-Zero Entries]
        println(out, "$max_row $max_col $num_entries")
        
        # Write Entry Lines
        for i in 1:num_entries
            println(out, "$(rows[i]) $(cols[i]) $(vals[i])")
        end
    end
    
    println("Successfully converted to $output_filename")
    println("Dimensions: $max_row x $max_col with $num_entries non-zeros.")
end

# Example usage:
# convert_to_mtx("getdp_dump.txt", "matrix.mtx")
using LinearAlgebra
using SparseArrays

"""
    CalculateNormUsingProbeVector(A::SparseMatrixCSC, c::SuiteSparse.CHOLMOD.Factor)

Calculates the relative residual error norm of a Cholesky factorization `c` 
against its source sparse matrix `A` using an implicit matrix-vector probe 
to prevent OutOfMemoryError.
"""
function CalculateNormUsingProbeVector(A, c)
    n = size(A, 1)
    
    # 1. Safely extract the lower triangular factor component
    L = sparse(c.L)
    p = c.p
    
    # 2. Generate a random probe vector to test the identity
    x = rand(n)
    
    # 3. Compute y1 = A * x
    y1 = A * x
    
    # 4. Compute y2 = (P' * L * L' * P) * x sequentially 
    # This prevents the creation of a dense matrix.
    Px = x[p]          # Apply forward permutation (P * x)
    LtPx = L' * Px     # Compute (L' * P * x)
    LLtPx = L * LtPx   # Compute (L * L' * P * x)
    
    # Apply inverse permutation to match A's indexing space (P' * LLtPx)
    y2 = zeros(eltype(A), n)
    y2[p] = LLtPx      
    
    # 5. Compute and return the relative norm error
    return norm(y1 - y2) / norm(y1)
end
