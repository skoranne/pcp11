################################################################################
# File   : SolveSparse.jl
# Author : Sandeep Koranne (C) 2026. All rights reserved.
# Purpose: Read a sparse matrix (orig) from GetDP and solve the Ax=b sparse system on GPU

using CUDA
using CUDA.CUSPARSE
using CUDSS
import LinearSolve as LS
import LinearAlgebra as LA

function LoadMatrix(fileName)
    A = mmread(fileName)
    A_gpu = CuSparseMatrixCSR(A)
    x = rand(size(A,1))
    b = A*x
    b_gpu = CuVector(b)
    (A,A_gpu,x,b,b_gpu)
end
const T = Float32

function SolveProblem(A,b)
    prob = LS.LinearProblem(A,b)
    sol  = LS.solve(prob, LS.LUFactorization())
    return sol
end

function SolveProblemFactorizationDoesNOTWORK(A,b)
    prob = LS.LinearProblem(A,b)
    sol=LS.solve(prob, LS.CUSOLVERRFFactorization())
    return sol
end

function SolveProblemKLU(A,b)
    prob = LS.LinearProblem(A,b)
    # Use KLU for symbolic factorization
    sol = LS.solve(prob, LS.CUSOLVERRFFactorization(symbolic = :KLU))
    # Reuse symbolic factorization for better performance
    sol = LS.solve(prob, LS.CUSOLVERRFFactorization(reuse_symbolic = true))
    return sol
end
