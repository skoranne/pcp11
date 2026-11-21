module CurvilinearTSV

using StaticArrays
using SpecialFunctions
using Krylov
using ArnoldiMethod
using FastGaussQuadrature
using JSON3

import LinearAlgebra: mul!, ldiv!
import SparseArrays: spzeros

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
end # module Backends

module Geometry
using ..StaticArrays

# 2D smooth C^∞ curvilinear contours for TSV cross-sections, teardrops, and offset SiO2 liners
# Star-shaped contours with Fourier series have unique solutions for small voltage perturbations
# due to convexity and monotonicity properties (MEMS PDE theory)

abstract type AbstractCurve2D end

"""
    Circle2D(center, R)

Circular contour of radius `R` centered at `center` (CCW parameterization).
"""
struct Circle2D <: AbstractCurve2D
    center::SVector{2,Float64}
    R::Float64
end
Circle2D(center::SVector{2,Float64}, R::Real) = Circle2D(center, Float64(R))
Circle2D(cx::Real, cy::Real, R::Real) = Circle2D(SVector{2,Float64}(cx, cy), Float64(R))
Circle2D(R::Real) = Circle2D(0.0, 0.0, R)

"""
    Ellipse2D(center, a, b, theta0=0.0)

Rotated ellipse with semi-axes `(a, b)` and rotation angle `theta0`.
"""
struct Ellipse2D <: AbstractCurve2D
    center::SVector{2,Float64}
    a::Float64
    b::Float64
    theta0::Float64
end
Ellipse2D(center::SVector{2,Float64}, a::Real, b::Real, theta0::Real=0.0) = 
    Ellipse2D(center, Float64(a), Float64(b), Float64(theta0))
Ellipse2D(cx::Real, cy::Real, a::Real, b::Real, theta0::Real=0.0) = 
    Ellipse2D(SVector{2,Float64}(cx, cy), Float64(a), Float64(b), Float64(theta0))
Ellipse2D(a::Real, b::Real) = Ellipse2D(0.0, 0.0, a, b, 0.0)

"""
    Teardrop2D(center, R, h)

Teardrop-shaped contour with base radius `R` and peak height `h`.
"""
struct Teardrop2D <: AbstractCurve2D
    center::SVector{2,Float64}
    R::Float64
    h::Float64
end
Teardrop2D(center::SVector{2,Float64}, R::Real, h::Real) = 
    Teardrop2D(center, Float64(R), Float64(h))
Teardrop2D(cx::Real, cy::Real, R::Real, h::Real) = 
    Teardrop2D(SVector{2,Float64}(cx, cy), Float64(R), Float64(h))
Teardrop2D(R::Real, h::Real) = Teardrop2D(0.0, 0.0, R, h)

"""
    FourierStar2D(center, R, coefficients)

Star-shaped contour defined by Fourier series coefficients.
Unique solutions exist for small voltage perturbations due to convexity and monotonicity.
"""
struct FourierStar2D <: AbstractCurve2D
    center::SVector{2,Float64}
    R::Float64
    coefficients::SVector{<:Any,Float64}
end
FourierStar2D(center::SVector{2,Float64}, R::Real, coefficients::AbstractVector{<:Real}) = 
    FourierStar2D(center, Float64(R), SVector{length(coefficients),Float64}(coefficients))
FourierStar2D(cx::Real, cy::Real, R::Real, coefficients::AbstractVector{<:Real}) = 
    FourierStar2D(SVector{2,Float64}(cx, cy), Float64(R), SVector{length(coefficients),Float64}(coefficients))
FourierStar2D(R::Real, coefficients::AbstractVector{<:Real}) = 
    FourierStar2D(0.0, 0.0, R, coefficients)

"""
    sample_circle_2d(circle, n)

Sample `n` equally spaced points on a circle.
"""
function sample_circle_2d(circle::Circle2D, n::Int)
    θs = LinRange(0, 2π, n+1)[1:end-1]
    return [SVector{2,Float64}(circle.center[1] + circle.R*cos(θ), 
                              circle.center[2] + circle.R*sin(θ)) for θ in θs]
end

"""
    sample_ellipse_2d(ellipse, n)

Sample `n` equally spaced points on an ellipse.
"""
function sample_ellipse_2d(ellipse::Ellipse2D, n::Int)
    θs = LinRange(0, 2π, n+1)[1:end-1]
    cosθ = cos.(θs)
    sinθ = sin.(θs)
    x = ellipse.center[1] .+ ellipse.a * cosθ * cos(ellipse.theta0) .- ellipse.b * sinθ * sin(ellipse.theta0)
    y = ellipse.center[2] .+ ellipse.a * cosθ * sin(ellipse.theta0) .+ ellipse.b * sinθ * cos(ellipse.theta0)
    return [SVector{2,Float64}(x[i], y[i]) for i in 1:n]
end

"""
    sample_teardrop_2d(teardrop, n)

Sample `n` points on a teardrop contour.
"""
function sample_teardrop_2d(teardrop::Teardrop2D, n::Int)
    θs = LinRange(0, 2π, n+1)[1:end-1]
    pts = SVector{2,Float64}[]
    for θ in θs
        r = teardrop.R + teardrop.h * (1 - cos(θ))/2
        push!(pts, SVector{2,Float64}(teardrop.center[1] + r*cos(θ), 
                                      teardrop.center[2] + r*sin(θ)))
    end
    return pts
end

# 3D surfaces for TSV structures

abstract type AbstractSurface3D end

"""
    Sphere3D(center, R)

Spherical surface of radius `R` centered at `center`.
"""
struct Sphere3D <: AbstractSurface3D
    center::SVector{3,Float64}
    R::Float64
end
Sphere3D(center::SVector{3,Float64}, R::Real) = Sphere3D(center, Float64(R))
Sphere3D(cx::Real, cy::Real, cz::Real, R::Real) = Sphere3D(SVector{3,Float64}(cx, cy, cz), Float64(R))
Sphere3D(R::Real) = Sphere3D(0.0, 0.0, 0.0, R)

"""
    sample_sphere_3d(sphere, n)

Sample `n` points on a sphere surface using Fibonacci lattice.
"""
function sample_sphere_3d(sphere::Sphere3D, n::Int)
    pts = SVector{3,Float64}[]
    golden_angle = π * (3 - sqrt(5))
    for i in 1:n
        z = 1 - 2 * (i - 0.5) / n
        θ = 2π * (i - 0.5) / n
        r = sqrt(1 - z^2)
        push!(pts, SVector{3,Float64}(sphere.center[1] + sphere.R * r * cos(θ),
                                      sphere.center[2] + sphere.R * r * sin(θ),
                                      sphere.center[3] + sphere.R * z))
    end
    return pts
end

export AbstractCurve2D, Circle2D, Ellipse2D, Teardrop2D, FourierStar2D
export sample_circle_2d, sample_ellipse_2d, sample_teardrop_2d
export AbstractSurface3D, Sphere3D, sample_sphere_3d

end # module Geometry

module Quadrature
# Quadrature rules for boundary integral equations
end

module Kernels
# Kernel functions and assembly routines
end

module Physics
# Physics models for electrostatics
end

# Re-export submodules
export Backends, Geometry, Quadrature, Kernels, Physics

end
