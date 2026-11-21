# 3D high-order curvilinear patch quadrature:
# - Far field: tensor-product Gauss-Legendre
# - Self-singular (1/r): 4-triangle Duffy polar singularity cancellation
# - Near-singular (thin 3D liners): recursive 2x2 sub-panel Gauss-Legendre refinement

"""
    lagrange_basis_1d(nodes, x) -> Vector{Float64}

Evaluates all 1D Lagrange interpolation basis polynomials at `x`.
"""
function lagrange_basis_1d(nodes::AbstractVector{Float64}, x::Float64)
    p = length(nodes)
    L = ones(Float64, p)
    @inbounds for j in 1:p
        xj = nodes[j]
        num = 1.0
        den = 1.0
        for m in 1:p
            if m != j
                num *= (x - nodes[m])
                den *= (xj - nodes[m])
            end
        end
        L[j] = num / den
    end
    return L
end

"""
    duffy_self_panel_weights!(row_slp, row_dlp, row_adlp,
                              patch, gl_nodes, u0, v0, x_tgt, n_tgt, kernel_fn; n_duffy=8)

Integrates weakly singular 3D Laplace/Helmholtz kernels over the self-panel [-1, 1]^2
centered at collocation point `(u0, v0)` by splitting into 4 sub-triangles and applying
the exact Duffy map (Jacobian cancels 1/r identically), projecting onto the tensor-product
Lagrange basis of the panel's `order × order` Gauss-Legendre nodes.
"""
function duffy_self_panel_weights!(row_slp::AbstractVector{T},
                                   row_dlp::Union{Nothing,AbstractVector{T}},
                                   row_adlp::Union{Nothing,AbstractVector{T}},
                                   patch::CurvedPatch3D,
                                   gl_nodes::Vector{Float64},
                                   u0::Float64, v0::Float64,
                                   x_tgt::SVector{3,Float64},
                                   n_tgt::SVector{3,Float64},
                                   kernel_eval;
                                   n_duffy::Int=8) where {T}
    p = length(gl_nodes)
    fill!(row_slp, zero(T))
    row_dlp !== nothing && fill!(row_dlp, zero(T))
    row_adlp !== nothing && fill!(row_adlp, zero(T))

    # Map [0, 1] Gauss-Legendre rules
    qx, qw = gausslegendre(n_duffy)
    s_nodes = 0.5 .* (qx .+ 1.0)
    s_weights = 0.5 .* qw

    # 4 sub-triangles from apex P0 = (u0, v0) to edges of [-1, 1]^2
    # Each triangle has vertices P0, P1, P2; Duffy map (a, b) ∈ [0,1]^2:
    #   (u, v) = P0 + a * ((1 - b)*(P1 - P0) + b*(P2 - P0))
    #   area Jacobian d(u,v)/d(a,b) = a * |(P1 - P0) × (P2 - P0)|
    corners = (SVector(-1.0, -1.0), SVector(1.0, -1.0), SVector(1.0, 1.0), SVector(-1.0, 1.0))
    P0 = SVector(u0, v0)

    for tri in 1:4
        P1 = corners[tri]
        P2 = corners[mod1(tri + 1, 4)]
        d1 = P1 - P0
        d2 = P2 - P0
        tri_area2 = abs(d1[1] * d2[2] - d1[2] * d2[1])
        tri_area2 < 1e-15 && continue

        for ia in 1:n_duffy
            a = s_nodes[ia]
            wa = s_weights[ia]
            for ib in 1:n_duffy
                b = s_nodes[ib]
                wb = s_weights[ib]
                uv = P0 + a * ((1.0 - b) * d1 + b * d2)
                u, v = uv[1], uv[2]
                pos, nrm, jac = patch.map_fn(u, v)
                w_tot = wa * wb * a * tri_area2 * jac
                g_val, dlp_val, adlp_val = kernel_eval(x_tgt, n_tgt, pos, nrm)
                Lu = lagrange_basis_1d(gl_nodes, u)
                Lv = lagrange_basis_1d(gl_nodes, v)
                idx = 0
                @inbounds for i in 1:p, j in 1:p
                    idx += 1
                    phi = Lu[i] * Lv[j] * w_tot
                    row_slp[idx] += g_val * phi
                    if row_dlp !== nothing
                        row_dlp[idx] += dlp_val * phi
                    end
                    if row_adlp !== nothing
                        row_adlp[idx] += adlp_val * phi
                    end
                end
            end
        end
    end
    return nothing
end

"""
    adaptive_near_panel_weights!(row_slp, row_dlp, row_adlp,
                                 patch, gl_nodes, x_tgt, n_tgt, kernel_eval; n_sub=2, q_ord=6)

Evaluates near-singular off-panel interactions (crucial for thin 3D SiO2 liners where
target distance d < 1.8 h_panel) via `n_sub × n_sub` sub-panel Gauss-Legendre quadrature
projected onto the panel's Lagrange nodes.
"""
function adaptive_near_panel_weights!(row_slp::AbstractVector{T},
                                      row_dlp::Union{Nothing,AbstractVector{T}},
                                      row_adlp::Union{Nothing,AbstractVector{T}},
                                      patch::CurvedPatch3D,
                                      gl_nodes::Vector{Float64},
                                      x_tgt::SVector{3,Float64},
                                      n_tgt::SVector{3,Float64},
                                      kernel_eval;
                                      n_sub::Int=2, q_ord::Int=6) where {T}
    p = length(gl_nodes)
    fill!(row_slp, zero(T))
    row_dlp !== nothing && fill!(row_dlp, zero(T))
    row_adlp !== nothing && fill!(row_adlp, zero(T))

    qx, qw = gausslegendre(q_ord)
    du = 2.0 / n_sub
    for sx in 0:n_sub-1, sy in 0:n_sub-1
        u_lo = -1.0 + sx * du
        v_lo = -1.0 + sy * du
        for ia in 1:q_ord
            u = u_lo + 0.5 * du * (qx[ia] + 1.0)
            wu = 0.5 * du * qw[ia]
            Lu = lagrange_basis_1d(gl_nodes, u)
            for ib in 1:q_ord
                v = v_lo + 0.5 * du * (qx[ib] + 1.0)
                wv = 0.5 * du * qw[ib]
                Lv = lagrange_basis_1d(gl_nodes, v)
                pos, nrm, jac = patch.map_fn(u, v)
                w_tot = wu * wv * jac
                g_val, dlp_val, adlp_val = kernel_eval(x_tgt, n_tgt, pos, nrm)
                idx = 0
                @inbounds for i in 1:p, j in 1:p
                    idx += 1
                    phi = Lu[i] * Lv[j] * w_tot
                    row_slp[idx] += g_val * phi
                    if row_dlp !== nothing
                        row_dlp[idx] += dlp_val * phi
                    end
                    if row_adlp !== nothing
                        row_adlp[idx] += adlp_val * phi
                    end
                end
            end
        end
    end
    return nothing
end
