module CUDAHelpers

using CUDA
using CurvilinearTSV

export CUDAHelpers

function gpu_allocate!(arr::CuArray, val::Real)
    arr .= val
    return arr
end

function gpu_fill!(arr::CuArray, val::Real)
    arr .= val
    return arr
end

function gpu_mul!(dest::CuArray, src::CuArray)
    dest .= src
    return dest
end

function gpu_axpy!(a::Real, x::CuArray, y::CuArray)
    y .= y .+ a .* x
    return y
end

function gpu_dot(x::CuArray, y::CuArray)
    return dot(x, y)
end

function gpu_norm(x::CuArray)
    return norm(x)
end

function gpu_gemv!(trans::Char, alpha::Real, A::CuMatrix, x::CuVector, beta::Real, y::CuVector)
    y .= alpha * (A * x) + beta * y
    return y
end

function gpu_gemm!(transA::Char, transB::Char, alpha::Real, A::CuMatrix, B::CuMatrix, beta::Real, C::CuMatrix)
    C .= alpha * (A * B) + beta * C
    return C
end

function gpu_lufact!(A::CuMatrix)
    return lu(A)
end

function gpu_linsolve!(A::CuMatrix, b::CuVector)
    return A \ b
end

function gpu_bicgstab!(A, b; kwargs...)
    return bicgstab(A, b; kwargs...)
end

function gpu_cublas_gemm!(transA::Char, transB::Char, alpha::Real, A::CuMatrix, B::CuMatrix, beta::Real, C::CuMatrix)
    C .= alpha * (A * B) + beta * C
    return C
end

function gpu_cublas_gemv!(trans::Char, alpha::Real, A::CuMatrix, x::CuVector, beta::Real, y::CuVector)
    y .= alpha * (A * x) + beta * y
    return y
end

function gpu_cublas_dot(x::CuVector, y::CuVector)
    return dot(x, y)
end

function gpu_cublas_axpy!(alpha::Real, x::CuVector, y::CuVector)
    y .= y .+ alpha .* x
    return y
end

function gpu_cublas_nrm2(x::CuVector)
    return norm(x)
end

function gpu_cublas_copy!(x::CuVector, y::CuVector)
    y .= x
    return y
end

function gpu_cublas_scal!(alpha::Real, x::CuVector)
    x .= alpha .* x
    return x
end

function gpu_cublas_set_matrix!(A::CuMatrix, B)
    A .= B
    return A
end

function gpu_cublas_get_matrix!(A::CuMatrix, B)
    B .= A
    return B
end

function gpu_cublas_set_vector!(x::CuVector, v)
    x .= v
    return x
end

function gpu_cublas_get_vector!(x::CuVector, v)
    v .= x
    return v
end

function gpu_cublas_create!(handle)
    return CuBLAS.cublasCreate()
end

function gpu_cublas_destroy!(handle)
    CuBLAS.cublasDestroy(handle)
    return nothing
end

function gpu_cublas_set_stream!(handle, stream)
    CuBLAS.cublasSetStream(handle, stream)
    return nothing
end

function gpu_cublas_get_stream!(handle)
    return CuBLAS.cublasGetStream(handle)
end

function gpu_cublas_sgemm!(transA::Char, transB::Char, m::Int, n::Int, k::Int, alpha::Real, A::CuMatrix, lda::Int, B::CuMatrix, ldb::Int, beta::Real, C::CuMatrix, ldc::Int)
    C .= alpha * (A * B) + beta * C
    return C
end

function gpu_cublas_sgemv!(trans::Char, m::Int, n::Int, alpha::Real, A::CuMatrix, lda::Int, x::CuVector, incx::Int, beta::Real, y::CuVector, incy::Int)
    y .= alpha * (A * x) + beta * y
    return y
end

function gpu_cublas_sdot!(n::Int, x::CuVector, incx::Int, y::CuVector, incy::Int)
    return dot(x[1:incy:end], y[1:incy:end])
end

function gpu_cublas_saxpy!(n::Int, alpha::Real, x::CuVector, incx::Int, y::CuVector, incy::Int)
    y[1:incy:end] .= y[1:incy:end] .+ alpha .* x[1:incx:end]
    return y
end

function gpu_cublas_snrm2!(n::Int, x::CuVector, incx::Int)
    return norm(x[1:incx:end])
end

function gpu_cublas_scopy!(n::Int, x::CuVector, incx::Int, y::CuVector, incy::Int)
    y[1:incy:end] .= x[1:incx:end]
    return y
end

function gpu_cublas_sscal!(n::Int, alpha::Real, x::CuVector, incx::Int)
    x[1:incx:end] .= alpha .* x[1:incx:end]
    return x
end

function gpu_cublas_ssetmatrix!(m::Int, n::Int, A::CuMatrix, lda::Int, B, ldb::Int)
    A .= B
    return A
end

function gpu_cublas_sgetmatrix!(m::Int, n::Int, A::CuMatrix, lda::Int, B, ldb::Int)
    B .= A
    return B
end

function gpu_cublas_ssetvector!(n::Int, x::CuVector, incx::Int, y, incy::Int)
    y[1:incy:end] .= x[1:incx:end]
    return y
end

function gpu_cublas_sgetvector!(n::Int, x::CuVector, incx::Int, y, incy::Int)
    y[1:incy:end] .= x[1:incx:end]
    return y
end

function gpu_cublas_sgemm_batched!(handle, transA::Char, transB::Char, m::Int, n::Int, k::Int, alpha::Real, A::CuVector{CuMatrix}, lda::Int, B::CuVector{CuMatrix}, ldb::Int, beta::Real, C::CuVector{CuMatrix}, ldc::Int, batchCount::Int)
    for i in 1:batchCount
        C[i] .= alpha * (A[i] * B[i]) + beta * C[i]
    end
    return C
end

function gpu_cublas_sgemm_strided_batched!(handle, transA::Char, transB::Char, m::Int, n::Int, k::Int, alpha::Real, A::CuMatrix, lda::Int, strideA::Int64, B::CuMatrix, ldb::Int, strideB::Int64, beta::Real, C::CuMatrix, ldc::Int, strideC::Int64, batchCount::Int)
    for i in 1:batchCount
        offset = (i-1) * strideA
        Ao = A[offset.+.(1:size(A,1), 1:size(A,2))...]
        Bo = B[offset.+.(1:size(B,1), 1:size(B,2))...]
        Co = C[offset.+.(1:size(C,1), 1:size(C,2))...]
        Co .= alpha * (Ao * Bo) + beta * Co
    end
    return C
end

function gpu_cublas_sgemv_batched!(handle, trans::Char, m::Int, n::Int, alpha::Real, A::CuVector{CuMatrix}, lda::Int, x::CuVector{CuVector}, incx::Int, beta::Real, y::CuVector{CuVector}, incy::Int, batchCount::Int)
    for i in 1:batchCount
        y[i] .= alpha * (A[i] * x[i]) + beta * y[i]
    end
    return y
end

function gpu_cublas_sgemv_strided_batched!(handle, trans::Char, m::Int, n::Int, alpha::Real, A::CuMatrix, lda::Int, strideA::Int64, x::CuVector, incx::Int, strideX::Int64, beta::Real, y::CuVector, incy::Int, strideY::Int64, batchCount::Int)
    for i in 1:batchCount
        offsetA = (i-1) * strideA
        offsetx = (i-1) * strideX
        offsety = (i-1) * strideY
        Ao = A[offsetA.+.(1:size(A,1), 1:size(A,2))...]
        xo = x[offsetx.+.(1:incy:end)]
        yo = y[offsety.+.(1:incy:end)]
        yo .= alpha * (Ao * xo) + beta * yo
    end
    return y
end

function gpu_cublas_sgetrf_batched!(handle, m::Int, n::Int, A::CuVector{CuMatrix}, lda::Int, ipiv::CuVector{CuVector{Int}}, info::CuVector{Int}, batchCount::Int)
    for i in 1:batchCount
        lu!(A[i])
    end
    return A
end

function gpu_cublas_sgetrs_batched!(handle, trans::Char, n::Int, nrhs::Int, A::CuVector{CuMatrix}, lda::Int, ipiv::CuVector{CuVector{Int}}, B::CuVector{CuVector}, ldb::Int, info::CuVector{Int}, batchCount::Int)
    for i in 1:batchCount
        B[i] = A[i] \ B[i]
    end
    return B
end

function gpu_cublas_sgesv!(handle, n::Int, nrhs::Int, A::CuMatrix, lda::Int, ipiv::CuVector{Int}, B::CuVector, ldb::Int, info::CuVector{Int})
    B = A \ B
    return B
end

function gpu_cublas_sgesv_batched!(handle, n::Int, nrhs::Int, A::CuVector{CuMatrix}, lda::Int, ipiv::CuVector{CuVector{Int}}, B::CuVector{CuVector}, ldb::Int, info::CuVector{Int}, batchCount::Int)
    for i in 1:batchCount
        B[i] = A[i] \ B[i]
    end
    return B
end

function gpu_cublas_sgesv_strided_batched!(handle, n::Int, nrhs::Int, A::CuMatrix, lda::Int, strideA::Int64, ipiv::CuVector{Int}, strideP::Int64, B::CuVector, ldb::Int, strideB::Int64, info::CuVector{Int}, batchCount::Int)
    for i in 1:batchCount
        offsetA = (i-1) * strideA
        offsetB = (i-1) * strideB
        Ao = A[offsetA.+.(1:size(A,1), 1:size(A,2))...]
        Bo = B[offsetB.+.(1:ldb:end)]
        Bo = Ao \ Bo
    end
    return B
end

end
