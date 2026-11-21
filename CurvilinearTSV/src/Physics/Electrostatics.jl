# 2D & 3D Electrostatics BIE Solvers:
# - Single-region & multi-conductor Laplace capacitance matrix extraction
# - Piecewise-homogeneous dielectric transmission BIE (Cu + SiO2 + Depleted Si + Lossy Si)
# - Constant-charge Poisson BIE via exact biharmonic volume-to-boundary reduction

"""
    solve_capacitance_2d_exterior(curves, eps_bg; N=64, backend=CPUBackend()) -> Matrix{Float64}

Computes the M×M Maxwell capacitance matrix [F/m] for M disjoint 2D curvilinear conductors
in a uniform dielectric `eps_bg` with far-field charge neutrality Σ Q_k = 0 (M >= 2).
Uses spectral Kress SLP quadrature with augmented scalar potential offset φ_∞.
"""
function solve_capacitance_2d_exterior(curves::AbstractVector{<:AbstractCurve2D}, eps_bg::Float64;
                                       N::Int=64, backend::AbstractBackend=CPUBackend())
    M = length(curves)
    M >= 2 || throw(ArgumentError("2D exterior capacitance requires >= 2 conductors for charge neutrality"))
    dcurves = [discretize(c, N) for c in curves]
    ntot = M * N
    A = zeros(Float64, ntot + 1, ntot + 1)

    for a in 1:M
        ra = ((a - 1) * N + 1):(a * N)
        for b in 1:M
            rb = ((b - 1) * N + 1):(b * N)
            A[ra, rb] = assemble_slp_2d(dcurves[a], dcurves[b], backend; self_curve=(a == b))
        end
        A[ra, ntot + 1] .= 1.0
        A[ntot + 1, ra] .= dcurves[a].weights
    end
    A[ntot + 1, ntot + 1] = 0.0

    # Right-hand sides: unit voltage on conductor k, 0 on others
    RHS = zeros(Float64, ntot + 1, M)
    for k in 1:M
        rk = ((k - 1) * N + 1):(k * N)
        RHS[rk, k] .= 1.0
    end
    Sol = solve_linear_system(A, RHS, backend)

    C = zeros(Float64, M, M)
    for k in 1:M, a in 1:M
        ra = ((a - 1) * N + 1):(a * N)
        # Q_a = eps_bg * ∫_{Γ_a} ρ_a ds
        C[a, k] = eps_bg * dot(dcurves[a].weights, @view Sol[ra, k])
    end
    return C
end

"""
    solve_capacitance_2d_coaxial(inner_curve, outer_curve, eps_diel; N=64, backend=CPUBackend())

Computes exact per-unit-length capacitance [F/m] between an inner curvilinear conductor
(at V=1) enclosed inside an outer curvilinear conductor (at V=0) filled with dielectric `eps_diel`.
Works for concentric, eccentric (bipolar benchmark), and curvilinear teardrop coaxial structures.
"""
function solve_capacitance_2d_coaxial(inner_curve::AbstractCurve2D, outer_curve::AbstractCurve2D,
                                      eps_diel::Float64; N::Int=64,
                                      backend::AbstractBackend=CPUBackend())
    dc_in = discretize(inner_curve, N)
    dc_out = discretize(outer_curve, N)
    # Represent φ(x) = φ_0 + V_in ρ_in + V_out ρ_out with ∫ ρ_in + ∫ ρ_out = 0
    A = zeros(Float64, 2N + 1, 2N + 1)
    A[1:N, 1:N] = assemble_slp_2d(dc_in, dc_in, backend; self_curve=true)
    A[1:N, N+1:2N] = assemble_slp_2d(dc_in, dc_out, backend; self_curve=false)
    A[N+1:2N, 1:N] = assemble_slp_2d(dc_out, dc_in, backend; self_curve=false)
    A[N+1:2N, N+1:2N] = assemble_slp_2d(dc_out, dc_out, backend; self_curve=true)
    A[1:2N, 2N+1] .= 1.0
    A[2N+1, 1:N] .= dc_in.weights
    A[2N+1, N+1:2N] .= dc_out.weights

    rhs = zeros(Float64, 2N + 1)
    rhs[1:N] .= 1.0
    sol = solve_linear_system(A, rhs, backend)
    Q_in = eps_diel * dot(dc_in.weights, @view sol[1:N])
    Q_out = eps_diel * dot(dc_out.weights, @view sol[N+1:2N])
    return (C=Q_in, Q_in=Q_in, Q_out=Q_out, rho_in=sol[1:N], rho_out=sol[N+1:2N],
            phi0=sol[2N+1], dc_in=dc_in, dc_out=dc_out)
end

"""
    solve_multidielectric_2d(conductors, voltages, interfaces, eps_in_list, eps_out_list, eps_adj;
                             N=64, backend=CPUBackend(), enclosed=false, Q_fixed=nothing)

Solves the 2D piecewise-homogeneous multi-dielectric Laplace transmission BIE (real or complex ε):
- On conductor k (Dirichlet): V ρ + φ_∞ = V_k
- On dielectric interface m (eps_a inside, eps_b outside):
      ((ε_a + ε_b)/2) ρ_m + (ε_a - ε_b) K^* ρ = Q_f
- Charge neutrality on conductors + bound interfaces: Σ ∫ ρ ds = 0 (or zero on conductors when enclosed).

Returns conductor charges Q_k, densities ρ, offset φ_∞, and potential evaluator.
"""
function solve_multidielectric_2d(conductors::AbstractVector{<:AbstractCurve2D},
                                  voltages::AbstractVector{<:Number},
                                  interfaces::AbstractVector{<:AbstractCurve2D},
                                  eps_in_list::AbstractVector{<:Number},
                                  eps_out_list::AbstractVector{<:Number},
                                  eps_adj::AbstractVector{<:Number};
                                  N::Int=64, backend::AbstractBackend=CPUBackend(),
                                  Q_fixed::Union{Nothing,AbstractVector{<:Number}}=nothing)
    Mc = length(conductors)
    Mi = length(interfaces)
    Mtot = Mc + Mi
    dcs = vcat([discretize(c, N) for c in conductors], [discretize(c, N) for c in interfaces])

    T = promote_type(Float64, eltype(voltages), eltype(eps_in_list), eltype(eps_out_list), eltype(eps_adj))
    ntot = Mtot * N
    A = zeros(T, ntot + 1, ntot + 1)
    rhs = zeros(T, ntot + 1)

    # Conductor rows (Dirichlet potential)
    for a in 1:Mc
        ra = ((a - 1) * N + 1):(a * N)
        for b in 1:Mtot
            rb = ((b - 1) * N + 1):(b * N)
            A[ra, rb] .= assemble_slp_2d(dcs[a], dcs[b], backend; self_curve=(a == b))
        end
        A[ra, ntot + 1] .= one(T)
        rhs[ra] .= T(voltages[a])
    end

    # Dielectric interface rows (flux jump transmission condition)
    for m in 1:Mi
        a = Mc + m
        ra = ((a - 1) * N + 1):(a * N)
        ea = T(eps_in_list[m])
        eb = T(eps_out_list[m])
        qf = Q_fixed === nothing ? zero(T) : T(Q_fixed[m])
        for b in 1:Mtot
            rb = ((b - 1) * N + 1):(b * N)
            Kmat = assemble_adlp_2d(dcs[a], dcs[b], backend; self_curve=(a == b))
            A[ra, rb] .= (ea - eb) .* Kmat
        end
        for j in 1:N
            A[ra[j], ra[j]] += 0.5 * (ea + eb)
        end
        rhs[ra] .= qf
    end

    # Total charge neutrality constraint (eliminates log|x| growth at ∞)
    for b in 1:Mtot
        rb = ((b - 1) * N + 1):(b * N)
        A[ntot + 1, rb] .= dcs[b].weights
    end
    A[ntot + 1, ntot + 1] = zero(T)
    rhs[ntot + 1] = zero(T)

    sol = solve_linear_system(A, rhs, backend)
    Q_cond = Vector{T}(undef, Mc)
    for a in 1:Mc
        ra = ((a - 1) * N + 1):(a * N)
        Q_cond[a] = T(eps_adj[a]) * dot(dcs[a].weights, @view sol[ra])
    end
    return (Q=Q_cond, rho=sol[1:ntot], phi0=sol[ntot + 1], dcurves=dcs, N=N, Mc=Mc, Mi=Mi)
end

"""
    biharmonic_volume_potential_2d(x_pts, outer_dc, inner_dc, rho0_over_eps)

Evaluates the exact 2D constant-charge Poisson volume potential
    φ_V(x) = (ρ_0 / ε) ∫_{Ω_dep} G_2D(x, y) dΩ_y
purely as boundary integrals via Green's second identity with biharmonic ψ_2D (Δψ_2D = G_2D):
    ∫_{Ω_dep} G_2D(x, y) dΩ_y = ∮_{Γ_outer} (∂ψ_2D/∂n_y) ds_y - ∮_{Γ_inner} (∂ψ_2D/∂n_y) ds_y
for annular/curvilinear depletion zone Ω_dep bounded externally by `outer_dc` and internally by `inner_dc`.
"""
function biharmonic_volume_potential_2d(x_pts::AbstractVector{SVector{2,Float64}},
                                        outer_dc::DiscretizedCurve2D,
                                        inner_dc::Union{Nothing,DiscretizedCurve2D},
                                        rho0_over_eps::Float64)
    M = length(x_pts)
    phiV = zeros(Float64, M)
    @inbounds for i in 1:M
        xi = x_pts[i]
        s = 0.0
        for j in 1:outer_dc.N
            s += dpsi_dn_2d(xi, outer_dc.pts[j], outer_dc.normals[j]) * outer_dc.weights[j]
        end
        if inner_dc !== nothing
            for j in 1:inner_dc.N
                s -= dpsi_dn_2d(xi, inner_dc.pts[j], inner_dc.normals[j]) * inner_dc.weights[j]
            end
        end
        phiV[i] = rho0_over_eps * s
    end
    return phiV
end

"""
    biharmonic_volume_normal_grad_2d(tgt_dc, outer_dc, inner_dc, rho0_over_eps)

Evaluates ∂φ_V/∂n_x on `tgt_dc` purely via weakly-singular boundary integrals:
    ∇_x ∫_{Ω_dep} G_2D(x, y) dΩ_y = -∮_{∂Ω_dep} G_2D(x, y) n_y ds_y
using Kress spectral quadrature when `tgt_dc` coincides with `outer_dc` or `inner_dc`.
"""
function biharmonic_volume_normal_grad_2d(tgt_dc::DiscretizedCurve2D,
                                          outer_dc::DiscretizedCurve2D,
                                          inner_dc::Union{Nothing,DiscretizedCurve2D},
                                          rho0_over_eps::Float64;
                                          backend::AbstractBackend=CPUBackend())
    function _signed_contrib(src::DiscretizedCurve2D, sgn::Float64)
        V = assemble_slp_2d(tgt_dc, src, backend; self_curve=(tgt_dc === src))
        # -sgn * (n_x ⋅ ∫ G(x, y) n_y ds_y)
        nx_c = [p[1] for p in src.normals]
        ny_c = [p[2] for p in src.normals]
        gx = V * nx_c
        gy = V * ny_c
        out = Vector{Float64}(undef, tgt_dc.N)
        @inbounds for i in 1:tgt_dc.N
            out[i] = -sgn * (tgt_dc.normals[i][1] * gx[i] + tgt_dc.normals[i][2] * gy[i])
        end
        return out
    end
    g = _signed_contrib(outer_dc, 1.0)
    if inner_dc !== nothing
        g .+= _signed_contrib(inner_dc, -1.0)
    end
    return rho0_over_eps .* g
end

"""
    solve_poisson_2d_annular(inner_curve, outer_curve, v_in, v_out, rho0, eps_diel; N=64)

Solves ∇²φ = -ρ_0 / ε in the 2D curvilinear region bounded by `inner_curve` (φ = v_in)
and `outer_curve` (φ = v_out) using exact biharmonic volume-to-boundary reduction.
"""
function solve_poisson_2d_annular(inner_curve::AbstractCurve2D, outer_curve::AbstractCurve2D,
                                  v_in::Float64, v_out::Float64, rho0::Float64, eps_diel::Float64;
                                  N::Int=64, backend::AbstractBackend=CPUBackend())
    dc_in = discretize(inner_curve, N)
    dc_out = discretize(outer_curve, N)
    c0 = rho0 / eps_diel
    # Particular solution φ_p(r) = -(ρ_0 / (4ε)) |r - r_0|^2 satisfies ∇²φ_p = -ρ_0/ε identically!
    # Also cross-checked against biharmonic_volume_potential_2d in regression tests.
    r0 = hasproperty(inner_curve, :center) ? inner_curve.center : SVector(0.0, 0.0)
    phi_p_in = [-0.25 * c0 * dot(p - r0, p - r0) for p in dc_in.pts]
    phi_p_out = [-0.25 * c0 * dot(p - r0, p - r0) for p in dc_out.pts]

    A = zeros(Float64, 2N + 1, 2N + 1)
    A[1:N, 1:N] = assemble_slp_2d(dc_in, dc_in, backend; self_curve=true)
    A[1:N, N+1:2N] = assemble_slp_2d(dc_in, dc_out, backend; self_curve=false)
    A[N+1:2N, 1:N] = assemble_slp_2d(dc_out, dc_in, backend; self_curve=false)
    A[N+1:2N, N+1:2N] = assemble_slp_2d(dc_out, dc_out, backend; self_curve=true)
    A[1:2N, 2N+1] .= 1.0
    A[2N+1, 1:N] .= dc_in.weights
    A[2N+1, N+1:2N] .= dc_out.weights

    rhs = zeros(Float64, 2N + 1)
    rhs[1:N] .= v_in .- phi_p_in
    rhs[N+1:2N] .= v_out .- phi_p_out
    sol = solve_linear_system(A, rhs, backend)
    rho_in, rho_out, phi0 = sol[1:N], sol[N+1:2N], sol[2N+1]

    # Normal derivatives on outer & inner boundaries (outward from each closed contour)
    K_in_in = assemble_adlp_2d(dc_in, dc_in, backend; self_curve=true)
    K_in_out = assemble_adlp_2d(dc_in, dc_out, backend; self_curve=false)
    K_out_in = assemble_adlp_2d(dc_out, dc_in, backend; self_curve=false)
    K_out_out = assemble_adlp_2d(dc_out, dc_out, backend; self_curve=true)

    # Inside annular domain: exterior to dc_in (∂n φ_h^+ = K^* ρ - 1/2 ρ_in),
    # interior to dc_out (∂n φ_h^- = K^* ρ + 1/2 ρ_out)
    dn_p_in = [-0.5 * c0 * dot(dc_in.pts[i] - r0, dc_in.normals[i]) for i in 1:N]
    dn_p_out = [-0.5 * c0 * dot(dc_out.pts[i] - r0, dc_out.normals[i]) for i in 1:N]
    dn_in = (K_in_in * rho_in + K_in_out * rho_out .- 0.5 .* rho_in) .+ dn_p_in
    dn_out = (K_out_in * rho_in + K_out_out * rho_out .+ 0.5 .* rho_out) .+ dn_p_out

    return (dc_in=dc_in, dc_out=dc_out, rho_in=rho_in, rho_out=rho_out, phi0=phi0,
            dn_in=dn_in, dn_out=dn_out, r0=r0, c0=c0)
end

# --- 3D Electrostatics BIE Solvers (Self-Capacitance, Spherical/TSV Capacitors, Multi-Conductor) ---

"""
    solve_capacitance_3d_isolated(mesh::SurfaceMesh3D, eps_bg; backend=CPUBackend(), n_duffy=8)

Computes the 3D isolated self-capacitance C [F] of a smooth curvilinear 3D conductor
(sphere, prolate spheroid, or 3D tapered teardrop TSV) at V = 1 in background `eps_bg`.
"""
function solve_capacitance_3d_isolated(mesh::SurfaceMesh3D, eps_bg::Float64;
                                       backend::AbstractBackend=CPUBackend(), n_duffy::Int=8)
    ops = assemble_bie_3d(mesh, mesh, backend; n_duffy=n_duffy)
    rhs = ones(Float64, length(mesh.pts))
    rho = solve_linear_system(ops.V, rhs, backend)
    Q = eps_bg * dot(mesh.weights, rho)
    return (C=Q, rho=rho, mesh=mesh)
end

"""
    solve_capacitance_3d_multi(meshes, eps_bg; backend=CPUBackend(), n_duffy=8) -> Matrix{Float64}

Computes the full M×M 3D Maxwell capacitance matrix [F] for M curvilinear 3D conductors
(e.g., concentric spheres, coaxial 3D teardrop TSV + outer liner/ground, or 3D HBM TSV array).
"""
function solve_capacitance_3d_multi(meshes::AbstractVector{SurfaceMesh3D}, eps_bg::Float64;
                                    backend::AbstractBackend=CPUBackend(), n_duffy::Int=8)
    M = length(meshes)
    counts = [length(m.pts) for m in meshes]
    offsets = cumsum(vcat(0, counts))
    ntot = offsets[end]
    V = zeros(Float64, ntot, ntot)
    for a in 1:M
        ra = (offsets[a] + 1):offsets[a + 1]
        for b in 1:M
            rb = (offsets[b] + 1):offsets[b + 1]
            ops = assemble_bie_3d(meshes[a], meshes[b], backend; n_duffy=n_duffy)
            V[ra, rb] = ops.V
        end
    end
    RHS = zeros(Float64, ntot, M)
    for k in 1:M
        rk = (offsets[k] + 1):offsets[k + 1]
        RHS[rk, k] .= 1.0
    end
    Sol = solve_linear_system(V, RHS, backend)
    C = zeros(Float64, M, M)
    for k in 1:M, a in 1:M
        ra = (offsets[a] + 1):offsets[a + 1]
        C[a, k] = eps_bg * dot(meshes[a].weights, @view Sol[ra, k])
    end
    return C
end

"""
    solve_poisson_3d_shell(mesh_in::SurfaceMesh3D, mesh_out::SurfaceMesh3D,
                           v_in, v_out, rho0, eps_diel; backend=CPUBackend(), n_duffy=8)

Solves 3D constant-charge Poisson ∇²φ = -ρ_0/ε inside a 3D curvilinear shell (e.g., 3D spherical
or 3D teardrop depletion shell) via exact particular solution φ_p = -(ρ_0/(6ε))|r - r_0|^2 + 3D BIE.
"""
function solve_poisson_3d_shell(mesh_in::SurfaceMesh3D, mesh_out::SurfaceMesh3D,
                                v_in::Float64, v_out::Float64, rho0::Float64, eps_diel::Float64;
                                backend::AbstractBackend=CPUBackend(), n_duffy::Int=8)
    nin, nout = length(mesh_in.pts), length(mesh_out.pts)
    c0 = rho0 / eps_diel
    r0 = sum(mesh_in.pts) / nin
    phi_p_in = [-c0 / 6.0 * dot(p - r0, p - r0) for p in mesh_in.pts]
    phi_p_out = [-c0 / 6.0 * dot(p - r0, p - r0) for p in mesh_out.pts]

    V = zeros(Float64, nin + nout, nin + nout)
    ops_ii = assemble_bie_3d(mesh_in, mesh_in, backend; need_adlp=true, n_duffy=n_duffy)
    ops_io = assemble_bie_3d(mesh_in, mesh_out, backend; need_adlp=true, n_duffy=n_duffy)
    ops_oi = assemble_bie_3d(mesh_out, mesh_in, backend; need_adlp=true, n_duffy=n_duffy)
    ops_oo = assemble_bie_3d(mesh_out, mesh_out, backend; need_adlp=true, n_duffy=n_duffy)

    V[1:nin, 1:nin] = ops_ii.V
    V[1:nin, nin+1:end] = ops_io.V
    V[nin+1:end, 1:nin] = ops_oi.V
    V[nin+1:end, nin+1:end] = ops_oo.V

    rhs = vcat(v_in .- phi_p_in, v_out .- phi_p_out)
    rho = solve_linear_system(V, rhs, backend)
    rho_in = rho[1:nin]
    rho_out = rho[nin+1:end]

    dn_p_in = [-(c0 / 3.0) * dot(mesh_in.pts[i] - r0, mesh_in.normals[i]) for i in 1:nin]
    dn_p_out = [-(c0 / 3.0) * dot(mesh_out.pts[i] - r0, mesh_out.normals[i]) for i in 1:nout]
    dn_in = (ops_ii.A * rho_in + ops_io.A * rho_out .- 0.5 .* rho_in) .+ dn_p_in
    dn_out = (ops_oi.A * rho_in + ops_oo.A * rho_out .+ 0.5 .* rho_out) .+ dn_p_out
    Q_in = -eps_diel * dot(mesh_in.weights, dn_in)
    return (Q_in=Q_in, dn_in=dn_in, dn_out=dn_out, rho_in=rho_in, rho_out=rho_out)
end
