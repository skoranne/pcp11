# Linear solvers (direct LU & Krylov GMRES) dispatched across CPU/CUDA/cuTile backends

function solve_linear_system(A::AbstractMatrix{T}, b::AbstractVector{T},
                             ::CPUBackend=CPUBackend();
                             method::Symbol=:auto, rtol::Float64=1e-12) where {T<:Number}
    n = size(A, 1)
    if method == :direct || (method == :auto && n <= 1500)
        return A \ b
    else
        x, stats = Krylov.gmres(A, b; rtol=rtol, atol=1e-14, itmax=min(n, 400))
        return x
    end
end

function solve_linear_system(A::AbstractMatrix{T}, B::AbstractMatrix{T},
                             ::CPUBackend=CPUBackend()) where {T<:Number}
    return lu(A) \ B
end
