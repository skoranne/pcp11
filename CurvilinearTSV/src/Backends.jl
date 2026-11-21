module Backends

# Backend abstraction for CPU, CUDA.jl, and cuTile.jl (Blackwell sm_121)

abstract type AbstractBackend end

"""
    CPUBackend(; tile_size=32, threaded=true)

Reference tiled CPU backend matching the power-of-two blocking layout of `CuTileBackend`.
"""
struct CPUBackend <: AbstractBackend
    tile_size::Int
    threaded::Bool
end
CPUBackend(; tile_size::Int=32, threaded::Bool=true) = CPUBackend(tile_size, threaded)

"""
    CUDABackend(; block_size=256, use_fp32=false)

CUDA.jl GPU backend activated via package extension `ext/CurvilinearTSVCUDAExt.jl`.
"""
struct CUDABackend <: AbstractBackend
    block_size::Int
    use_fp32::Bool
end
CUDABackend(; block_size::Int=256, use_fp32::Bool=false) = CUDABackend(block_size, use_fp32)

"""
    CuTileBackend(; tile_m=32, tile_n=32)

cuTile.jl (CUDA 13 Tile IR, Blackwell sm_121) backend activated via `ext/CurvilinearTSVcuTileExt.jl`.
"""
struct CuTileBackend <: AbstractBackend
    tile_m::Int
    tile_n::Int
end
CuTileBackend(; tile_m::Int=32, tile_n::Int=32) = CuTileBackend(tile_m, tile_n)

backend_available(::CPUBackend) = true
backend_available(::CUDABackend) = false
backend_available(::CuTileBackend) = false

default_backend() = CPUBackend()

export AbstractBackend, CPUBackend, CUDABackend, CuTileBackend, backend_available, default_backend

end
