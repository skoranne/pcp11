# Spectral 2D periodic quadrature (Kress / Martensen-Kussmaul log-singularity splitting)
# and Fourier-upsampled thin-liner close evaluation.

"""
    kress_log_weights(N::Int) -> Vector{Float64}

Exact Martensen-Kussmaul / Kress trigonometric quadrature weights W_k^(N) (k = 0..N-1) for
    ∫_0^{2π} (-1/(2π) ln(2 |sin((t_i - τ)/2)|)) ψ(τ) dτ ≈ Σ_{j=0}^{N-1} W_{|i-j|}^(N) ψ(τ_j)
for even `N`, giving exponential (spectral) convergence on smooth closed contours.
"""
function kress_log_weights(N::Int)
    iseven(N) || throw(ArgumentError("Kress quadrature requires even N"))
    W = Vector{Float64}(undef, N)
    halfN = N ÷ 2
    for k in 0:N-1
        s = 0.0
        for m in 1:halfN-1
            s += cos(2pi * m * k / N) / m
        end
        W[k+1] = s / N + (iseven(k) ? 1.0 : -1.0) / (N * N)
    end
    return W
end

"""
    fourier_upsample_matrix(N::Int, M_factor::Int) -> Matrix{Float64}

Trigonometric interpolation matrix mapping N equispaced periodic values to N*M_factor nodes,
used to achieve machine-precision close-evaluation across thin SiO2 liners (t_ox << R_0).
"""
function fourier_upsample_matrix(N::Int, M_factor::Int)
    M_factor == 1 && return Matrix{Float64}(I, N, N)
    iseven(N) || throw(ArgumentError("N must be even"))
    Nfine = N * M_factor
    P = Matrix{Float64}(undef, Nfine, N)
    halfN = N ÷ 2
    for i in 0:Nfine-1
        ti = 2pi * i / Nfine
        for j in 0:N-1
            tj = 2pi * j / N
            d = ti - tj
            val = 1.0 + cos(halfN * d)
            for m in 1:halfN-1
                val += 2.0 * cos(m * d)
            end
            P[i+1, j+1] = val / N
        end
    end
    return P
end

"""
    close_eval_upsample_factor(target_pts, src_curve::DiscretizedCurve2D; safety=6.0, max_factor=8) -> Int

Determines required Fourier oversampling factor M so that fine source spacing h_fine < d_min / safety.
"""
function close_eval_upsample_factor(target_pts::AbstractVector{SVector{2,Float64}},
                                    src::DiscretizedCurve2D;
                                    safety::Float64=5.0, max_factor::Int=8)
    h_max = maximum(src.weights)
    d_min = Inf
    for p in target_pts, q in src.pts
        d = norm(p - q)
        if d > 1e-14 && d < d_min
            d_min = d
        end
    end
    if !isfinite(d_min) || d_min >= safety * h_max
        return 1
    end
    fac = ceil(Int, safety * h_max / d_min)
    return clamp(fac, 1, max_factor)
end
