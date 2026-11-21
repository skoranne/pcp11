# Tiled CPU BIE operator assembly & matrix-free matvec (power-of-two tiles mirroring CuTileBackend)

"""
    assemble_slp_2d(tgt, src, backend::CPUBackend; self_curve=false) -> Matrix{Float64}

Assembles the 2D Laplace Single-Layer Potential matrix V_ij using:
- Spectral Martensen-Kussmaul / Kress log-singularity splitting when `self_curve == true`.
- Fourier-upsampled close-evaluation when `tgt` and `src` are closely spaced (thin SiO2 liners).
- Power-of-two tile-blocked evaluation matching `CuTileBackend`.
"""
function assemble_slp_2d(tgt::DiscretizedCurve2D, src::DiscretizedCurve2D,
                         backend::CPUBackend=CPUBackend();
                         self_curve::Bool=(tgt === src))
    M, N = tgt.N, src.N
    V = Matrix{Float64}(undef, M, N)
    if self_curve
        W_kress = kress_log_weights(N)
        dt = 2pi / N
        @inbounds for i in 1:M
            xi = tgt.pts[i]
            ti = tgt.t[i]
            for j in 1:N
                if i == j
                    # Smooth remainder L_smooth(t_i, t_i) = -1/(2π) ln|r'(t_i)|
                    l_smooth = -INV_2PI * log(src.speeds[j])
                    V[i, j] = (W_kress[1] + dt * l_smooth) * src.speeds[j]
                else
                    d = norm(xi - src.pts[j])
                    s_sin = 2.0 * abs(sin(0.5 * (ti - src.t[j])))
                    l_smooth = -INV_2PI * log(d / s_sin)
                    k_idx = mod(i - j, N) + 1
                    V[i, j] = (W_kress[k_idx] + dt * l_smooth) * src.speeds[j]
                end
            end
        end
        return V
    end

    # Distinct curves: check thin-liner proximity for Fourier upsampling
    ups = close_eval_upsample_factor(tgt.pts, src)
    if ups > 1
        src_fine = discretize(src.curve, N * ups)
        P = fourier_upsample_matrix(N, ups)
        V_fine = Matrix{Float64}(undef, M, src_fine.N)
        _assemble_slp_2d_tiled!(V_fine, tgt.pts, src_fine.pts, src_fine.weights, backend.tile_size)
        mul!(V, V_fine, P)
    else
        _assemble_slp_2d_tiled!(V, tgt.pts, src.pts, src.weights, backend.tile_size)
    end
    return V
end

function _assemble_slp_2d_tiled!(V::Matrix{Float64},
                                 tpts::Vector{SVector{2,Float64}},
                                 spts::Vector{SVector{2,Float64}},
                                 swts::Vector{Float64},
                                 B::Int)
    M, N = length(tpts), length(spts)
    @inbounds for i0 in 1:B:M
        i1 = min(i0 + B - 1, M)
        for j0 in 1:B:N
            j1 = min(j0 + B - 1, N)
            for j in j0:j1
                yj = spts[j]
                wj = swts[j]
                for i in i0:i1
                    V[i, j] = green_laplace_2d(tpts[i], yj) * wj
                end
            end
        end
    end
    return V
end

"""
    assemble_adlp_2d(tgt, src, backend::CPUBackend; self_curve=false) -> Matrix{Float64}

Assembles the 2D Laplace Adjoint Double-Layer matrix K^*_ij = (∂G/∂n_x)(x_i, y_j) ds_j,
with exact diagonal curvature limit `-κ(t_i)/(4π) ds_i` and thin-liner Fourier close-evaluation.
"""
function assemble_adlp_2d(tgt::DiscretizedCurve2D, src::DiscretizedCurve2D,
                          backend::CPUBackend=CPUBackend();
                          self_curve::Bool=(tgt === src))
    M, N = tgt.N, src.N
    K = Matrix{Float64}(undef, M, N)
    if self_curve
        @inbounds for i in 1:M
            xi, nxi = tgt.pts[i], tgt.normals[i]
            for j in 1:N
                if i == j
                    K[i, j] = -INV_4PI * tgt.curvatures[i] * src.weights[j]
                else
                    K[i, j] = adlp_laplace_2d(xi, nxi, src.pts[j]) * src.weights[j]
                end
            end
        end
        return K
    end

    ups = close_eval_upsample_factor(tgt.pts, src)
    if ups > 1
        src_fine = discretize(src.curve, N * ups)
        P = fourier_upsample_matrix(N, ups)
        K_fine = Matrix{Float64}(undef, M, src_fine.N)
        _assemble_adlp_2d_tiled!(K_fine, tgt.pts, tgt.normals, src_fine.pts, src_fine.weights, backend.tile_size)
        mul!(K, K_fine, P)
    else
        _assemble_adlp_2d_tiled!(K, tgt.pts, tgt.normals, src.pts, src.weights, backend.tile_size)
    end
    return K
end

function _assemble_adlp_2d_tiled!(K::Matrix{Float64},
                                  tpts::Vector{SVector{2,Float64}},
                                  tnrms::Vector{SVector{2,Float64}},
                                  spts::Vector{SVector{2,Float64}},
                                  swts::Vector{Float64},
                                  B::Int)
    M, N = length(tpts), length(spts)
    @inbounds for i0 in 1:B:M
        i1 = min(i0 + B - 1, M)
        for j0 in 1:B:N
            j1 = min(j0 + B - 1, N)
            for j in j0:j1
                yj = spts[j]
                wj = swts[j]
                for i in i0:i1
                    K[i, j] = adlp_laplace_2d(tpts[i], tnrms[i], yj) * wj
                end
            end
        end
    end
    return K
end

"""
    assemble_dlp_2d(tgt, src, backend::CPUBackend; self_curve=false) -> Matrix{Float64}

Assembles the 2D Laplace Double-Layer matrix K_ij = (∂G/∂n_y)(x_i, y_j) ds_j.
"""
function assemble_dlp_2d(tgt::DiscretizedCurve2D, src::DiscretizedCurve2D,
                         backend::CPUBackend=CPUBackend();
                         self_curve::Bool=(tgt === src))
    M, N = tgt.N, src.N
    K = Matrix{Float64}(undef, M, N)
    if self_curve
        @inbounds for i in 1:M
            xi = tgt.pts[i]
            for j in 1:N
                if i == j
                    K[i, j] = -INV_4PI * src.curvatures[j] * src.weights[j]
                else
                    K[i, j] = dlp_laplace_2d(xi, src.pts[j], src.normals[j]) * src.weights[j]
                end
            end
        end
        return K
    end
    ups = close_eval_upsample_factor(tgt.pts, src)
    if ups > 1
        src_fine = discretize(src.curve, N * ups)
        P = fourier_upsample_matrix(N, ups)
        K_fine = Matrix{Float64}(undef, M, src_fine.N)
        @inbounds for j in 1:src_fine.N, i in 1:M
            K_fine[i, j] = dlp_laplace_2d(tgt.pts[i], src_fine.pts[j], src_fine.normals[j]) * src_fine.weights[j]
        end
        mul!(K, K_fine, P)
    else
        @inbounds for j in 1:N, i in 1:M
            K[i, j] = dlp_laplace_2d(tgt.pts[i], src.pts[j], src.normals[j]) * src.weights[j]
        end
    end
    return K
end

"""
    tiled_slp_matvec_2d(tpts, spts, swts, x, backend::CPUBackend) -> Vector{Float64}

Matrix-free tiled single-layer potential matvec (off-diagonal / distinct curves),
sharing exact tile structure with `CuTileBackend`.
"""
function tiled_slp_matvec_2d(tpts::Vector{SVector{2,Float64}},
                             spts::Vector{SVector{2,Float64}},
                             swts::Vector{Float64},
                             x::Vector{Float64},
                             backend::CPUBackend=CPUBackend())
    M, N = length(tpts), length(spts)
    y = zeros(Float64, M)
    B = backend.tile_size
    @inbounds for i0 in 1:B:M
        i1 = min(i0 + B - 1, M)
        for j0 in 1:B:N
            j1 = min(j0 + B - 1, N)
            for j in j0:j1
                sj = spts[j]
                wx = swts[j] * x[j]
                for i in i0:i1
                    y[i] += green_laplace_2d(tpts[i], sj) * wx
                end
            end
        end
    end
    return y
end

"""
    assemble_bie_3d(tgt::SurfaceMesh3D, src::SurfaceMesh3D, backend::CPUBackend;
                    need_dlp=false, need_adlp=false, n_duffy=8, near_ratio=1.6,
                    k_wave::Union{Nothing,ComplexF64}=nothing)

Assembles 3D curvilinear BIE matrices (SLP, optional DLP, optional Adjoint DLP) using:
- Duffy 4-triangle singularity cancellation on self-panels
- Adaptive sub-panel Gauss-Legendre on near-singular panels (thin 3D liners)
- Direct tensor-product Gauss-Legendre on far-field panels.
"""
function assemble_bie_3d(tgt::SurfaceMesh3D, src::SurfaceMesh3D,
                         backend::CPUBackend=CPUBackend();
                         need_dlp::Bool=false, need_adlp::Bool=false,
                         n_duffy::Int=8, near_ratio::Float64=1.6,
                         k_wave::Union{Nothing,ComplexF64}=nothing)
    M, N = length(tgt.pts), length(src.pts)
    T = k_wave === nothing ? Float64 : ComplexF64
    V = zeros(T, M, N)
    D = need_dlp ? zeros(T, M, N) : nothing
    A = need_adlp ? zeros(T, M, N) : nothing

    same_mesh = (tgt === src)
    p = src.order
    p2 = p * p
    kfn = k_wave === nothing ? kernels_laplace_3d :
          (x, nx, y, ny) -> kernels_helmholtz_3d(x, nx, y, ny, k_wave)

    tmp_v = Vector{T}(undef, p2)
    tmp_d = need_dlp ? Vector{T}(undef, p2) : nothing
    tmp_a = need_adlp ? Vector{T}(undef, p2) : nothing

    for (pidx, patch) in enumerate(src.patches)
        col_start = (pidx - 1) * p2 + 1
        col_end = pidx * p2
        cols = col_start:col_end
        h_patch = patch.char_size
        pcenter = patch.center

        for i in 1:M
            xi = tgt.pts[i]
            nxi = tgt.normals[i]
            if same_mesh && tgt.patch_of_node[i] == pidx
                uv0 = tgt.local_uv[i]
                duffy_self_panel_weights!(tmp_v, tmp_d, tmp_a, patch, src.gl_nodes,
                                          uv0[1], uv0[2], xi, nxi, kfn; n_duffy=n_duffy)
                @inbounds for k in 1:p2
                    V[i, col_start + k - 1] = tmp_v[k]
                    need_dlp && (D[i, col_start + k - 1] = tmp_d[k])
                    need_adlp && (A[i, col_start + k - 1] = tmp_a[k])
                end
            else
                dist = norm(xi - pcenter)
                if dist < near_ratio * h_patch
                    nsub = dist < 0.6 * h_patch ? 3 : 2
                    adaptive_near_panel_weights!(tmp_v, tmp_d, tmp_a, patch, src.gl_nodes,
                                                 xi, nxi, kfn; n_sub=nsub, q_ord=max(p + 2, 6))
                    @inbounds for k in 1:p2
                        V[i, col_start + k - 1] = tmp_v[k]
                        need_dlp && (D[i, col_start + k - 1] = tmp_d[k])
                        need_adlp && (A[i, col_start + k - 1] = tmp_a[k])
                    end
                else
                    @inbounds for j in cols
                        g, dlp, adlp = kfn(xi, nxi, src.pts[j], src.normals[j])
                        wj = src.weights[j]
                        V[i, j] = g * wj
                        need_dlp && (D[i, j] = dlp * wj)
                        need_adlp && (A[i, j] = adlp * wj)
                    end
                end
            end
        end
    end
    return (V=V, D=D, A=A)
end
