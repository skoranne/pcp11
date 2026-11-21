# 2D smooth C^∞ curvilinear contours for TSV cross-sections, teardrops, and offset SiO2 liners

module Geometry

abstract type AbstractCurve2D end

"""
    Circle2D(center, R)

Circular contour of radius `R` centered at `center` (CCW parameterization).
"""
struct Circle2D <: AbstractCurve2D
    center::SVector{2,Float64}
    R::Float64
end
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
Teardrop2D(cx::Real, cy::Real, R::Real, h::Real) = 
    Teardrop2D(SVector{2,Float64}(cx, cy), Float64(R), Float64(h))
Teardrop2D(R::Real, h::Real) = Teardrop2D(0.0, 0.0, R, h)

"""
    FourierStar2D(center, R, coefficients)

Star-shaped contour defined by Fourier series coefficients.
"""
struct FourierStar2D <: AbstractCurve2D
    center::SVector{2,Float64}
    R::Float64
    coefficients::SVector{<:Any,Float64}
end
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
    x = ellipse.center[1] + ellipse.a * cosθ * cos(ellipse.theta0) - ellipse.b * sinθ * sin(ellipse.theta0)
    y = ellipse.center[2] + ellipse.a * cosθ * sin(ellipse.theta0) + ellipse.b * sinθ * cos(ellipse.theta0)
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

export AbstractCurve2D, Circle2D, Ellipse2D, Teardrop2D, FourierStar2D
export sample_circle_2d, sample_ellipse_2d, sample_teardrop_2d

end # module Geometry
