# Fundamental solutions (Green's functions) for 2D/3D Laplace, Biharmonic-Poisson, and Helmholtz PDEs
# Verified symbolically in mathematica/01_greens_functions_2d_3d.wl

const INV_2PI = 1.0 / (2.0 * pi)
const INV_4PI = 1.0 / (4.0 * pi)
const INV_8PI = 1.0 / (8.0 * pi)

# --- 2D Laplace Green's Function: G(x, y) = -1/(2π) ln|x - y| ---

@inline function green_laplace_2d(x::SVector{2,Float64}, y::SVector{2,Float64})
    r = norm(x - y)
    return -INV_2PI * log(r)
end

# ∂G/∂n_y = ((x - y) ⋅ n_y) / (2π |x - y|^2)
@inline function dlp_laplace_2d(x::SVector{2,Float64}, y::SVector{2,Float64}, ny::SVector{2,Float64})
    d = x - y
    r2 = dot(d, d)
    return INV_2PI * dot(d, ny) / r2
end

# ∂G/∂n_x = -((x - y) ⋅ n_x) / (2π |x - y|^2)
@inline function adlp_laplace_2d(x::SVector{2,Float64}, nx::SVector{2,Float64}, y::SVector{2,Float64})
    d = x - y
    r2 = dot(d, d)
    return -INV_2PI * dot(d, nx) / r2
end

# --- 2D Biharmonic Fundamental Solution for Constant-Charge Poisson Volume Reduction ---
# ψ_2D(r) = -r^2 (ln r - 1) / (8π) satisfies Δ_y ψ_2D = G_2D
# ∂ψ_2D/∂n_y = (d ⋅ n_y) * (2 ln r - 1) / (8π)  where d = x - y
@inline function dpsi_dn_2d(x::SVector{2,Float64}, y::SVector{2,Float64}, ny::SVector{2,Float64})
    d = x - y
    r2 = dot(d, d)
    r2 < 1e-30 && return 0.0
    # ψ(r) = -(r^2/8π)(ln r - 1), dψ/dr = -(r/(8π))(2 ln r - 1), grad_y r = -(x-y)/r
    # => ∂ψ/∂n_y = +(d ⋅ ny) * (2 ln r - 1) / (8π) = (d ⋅ ny) * (ln(r^2) - 1) / (8π)
    return INV_8PI * dot(d, ny) * (log(r2) - 1.0)
end

# --- 3D Laplace Green's Function: G(x, y) = 1 / (4π |x - y|) ---

@inline function green_laplace_3d(x::SVector{3,Float64}, y::SVector{3,Float64})
    return INV_4PI / norm(x - y)
end

@inline function kernels_laplace_3d(x::SVector{3,Float64}, nx::SVector{3,Float64},
                                    y::SVector{3,Float64}, ny::SVector{3,Float64})
    d = x - y
    r2 = dot(d, d)
    r = sqrt(r2)
    inv_r = 1.0 / r
    inv_4pi_r3 = INV_4PI * inv_r / r2
    g = INV_4PI * inv_r
    dlp = dot(d, ny) * inv_4pi_r3
    adlp = -dot(d, nx) * inv_4pi_r3
    return g, dlp, adlp
end

# --- 3D Biharmonic Fundamental Solution: ψ_3D(r) = r / (8π), Δ ψ_3D = G_3D ---
# ∂ψ_3D/∂n_y = -(x - y) ⋅ n_y / (8π |x - y|)
@inline function dpsi_dn_3d(x::SVector{3,Float64}, y::SVector{3,Float64}, ny::SVector{3,Float64})
    d = x - y
    r = norm(d)
    r < 1e-15 && return 0.0
    return -INV_8PI * dot(d, ny) / r
end

# --- 2D & 3D Helmholtz Green's Functions (Eddy-Current Skin Effect & 60 GHz Retardation) ---

# 2D Helmholtz: G_k(r) = (i/4) H_0^(1)(k r)  satisfying (Δ + k^2) G_k = -δ
# For eddy-current diffusion (Δ - γ^2) u = -f, set k = i γ => G(r) = 1/(2π) K_0(γ r)!
@inline function green_modified_helmholtz_2d(x::SVector{2,Float64}, y::SVector{2,Float64},
                                             gamma::ComplexF64)
    r = norm(x - y)
    return INV_2PI * besselk(0, gamma * r)
end

# ∂/∂n_y [ (1/2π) K_0(γ r) ] = (γ / 2π) K_1(γ r) * ((x - y) ⋅ n_y) / r
@inline function dlp_modified_helmholtz_2d(x::SVector{2,Float64}, y::SVector{2,Float64},
                                           ny::SVector{2,Float64}, gamma::ComplexF64)
    d = x - y
    r = norm(d)
    return (gamma * INV_2PI) * besselk(1, gamma * r) * (dot(d, ny) / r)
end

@inline function adlp_modified_helmholtz_2d(x::SVector{2,Float64}, nx::SVector{2,Float64},
                                            y::SVector{2,Float64}, gamma::ComplexF64)
    d = x - y
    r = norm(d)
    return -(gamma * INV_2PI) * besselk(1, gamma * r) * (dot(d, nx) / r)
end

# 3D Retarded Helmholtz Green's function: G_k(r) = exp(-i k r) / (4π r)  (exp(+iωt) convention)
@inline function kernels_helmholtz_3d(x::SVector{3,Float64}, nx::SVector{3,Float64},
                                      y::SVector{3,Float64}, ny::SVector{3,Float64},
                                      k::ComplexF64)
    d = x - y
    r = norm(d)
    ikr = im * k * r
    e = exp(-ikr)
    g = e * INV_4PI / r
    fac = (1.0 + ikr) * e * INV_4PI / (r * r * r)
    dlp = dot(d, ny) * fac
    adlp = -dot(d, nx) * fac
    return g, dlp, adlp
end
