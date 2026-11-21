# 3D curvilinear surfaces: spheres, prolate spheroids, tapered/scalloped teardrop TSVs,
# offset SiO2 liners, curvilinear quad-patch meshes, and Gmsh .msh import.

abstract type AbstractSurface3D end

struct Sphere3D <: AbstractSurface3D
    center::SVector{3,Float64}
    R::Float64
end
Sphere3D(R::Real; center=(0.0, 0.0, 0.0)) =
    Sphere3D(SVector{3,Float64}(center[1], center[2], center[3]), Float64(R))

struct ProlateSpheroid3D <: AbstractSurface3D
    center::SVector{3,Float64}
    a::Float64   # equatorial semi-axis (x, y)
    c::Float64   # polar semi-axis (z)
end
ProlateSpheroid3D(a::Real, c::Real; center=(0.0, 0.0, 0.0)) =
    ProlateSpheroid3D(SVector{3,Float64}(center[1], center[2], center[3]), Float64(a), Float64(c))

"""
    CurvilinearTSV3D

3D tapered & optionally Bosch-scalloped TSV surface swept from a 2D curvilinear cross-section
`profile` over z ∈ [zmin, zmax], with scale profile:
    s(z) = 1 - taper * ((z - zmin)/h) + scallop_amp * sin(2π * scallop_cycles * (z - zmin)/h).
"""
struct CurvilinearTSV3D{C<:AbstractCurve2D} <: AbstractSurface3D
    profile::C
    zmin::Float64
    zmax::Float64
    taper::Float64
    scallop_amp::Float64
    scallop_cycles::Int
end
CurvilinearTSV3D(profile::AbstractCurve2D, h::Real; zmin::Real=0.0, taper::Real=0.05,
                 scallop_amp::Real=0.0, scallop_cycles::Int=8) =
    CurvilinearTSV3D(profile, Float64(zmin), Float64(zmin + h), Float64(taper),
                     Float64(scallop_amp), scallop_cycles)

@inline function scale_profile(tsv::CurvilinearTSV3D, z::Float64)
    h = tsv.zmax - tsv.zmin
    ζ = (z - tsv.zmin) / h
    s = 1.0 - tsv.taper * ζ + tsv.scallop_amp * sin(2pi * tsv.scallop_cycles * ζ)
    ds = (-tsv.taper + 2pi * tsv.scallop_cycles * tsv.scallop_amp * cos(2pi * tsv.scallop_cycles * ζ)) / h
    return s, ds
end

"""
    CurvedPatch3D

Smooth parametric quadrilateral surface patch mapping (u, v) ∈ [-1, 1]^2 → ℝ^3,
storing analytic evaluator `map_fn(u, v) -> (pos, normal, jac_area)`.
"""
struct CurvedPatch3D{F}
    map_fn::F
    center::SVector{3,Float64}
    char_size::Float64   # approximate patch diameter h for near-singular detection
    region_id::Int
end

"""
    SurfaceMesh3D

Collection of curvilinear quadrilateral patches and their high-order Gauss-Legendre nodes.
"""
struct SurfaceMesh3D
    patches::Vector{CurvedPatch3D}
    order::Int
    gl_nodes::Vector{Float64}
    gl_weights::Vector{Float64}
    pts::Vector{SVector{3,Float64}}
    normals::Vector{SVector{3,Float64}}
    weights::Vector{Float64}
    patch_of_node::Vector{Int}
    local_uv::Vector{SVector{2,Float64}}
    region_of_node::Vector{Int}
end

function build_surface_mesh(patches::Vector{<:CurvedPatch3D}, order::Int)
    glx, glw = gausslegendre(order)
    np = length(patches)
    ntot = np * order * order
    pts = Vector{SVector{3,Float64}}(undef, ntot)
    normals = Vector{SVector{3,Float64}}(undef, ntot)
    weights = Vector{Float64}(undef, ntot)
    patch_of_node = Vector{Int}(undef, ntot)
    local_uv = Vector{SVector{2,Float64}}(undef, ntot)
    region_of_node = Vector{Int}(undef, ntot)
    idx = 0
    for (pidx, p) in enumerate(patches)
        for i in 1:order, j in 1:order
            idx += 1
            u, v = glx[i], glx[j]
            pos, nrm, jac = p.map_fn(u, v)
            pts[idx] = pos
            normals[idx] = nrm
            weights[idx] = glw[i] * glw[j] * jac
            patch_of_node[idx] = pidx
            local_uv[idx] = SVector(u, v)
            region_of_node[idx] = p.region_id
        end
    end
    return SurfaceMesh3D(Vector{CurvedPatch3D}(patches), order, glx, glw,
                         pts, normals, weights, patch_of_node, local_uv, region_of_node)
end

# --- Exact curvilinear cubed-sphere patches (singularity-free poles!) ---

function _cubed_sphere_patch_map(face::Int, u0::Float64, u1::Float64, v0::Float64, v1::Float64,
                                 center::SVector{3,Float64}, rx::Float64, ry::Float64, rz::Float64)
    du = 0.5 * (u1 - u0)
    dv = 0.5 * (v1 - v0)
    um = 0.5 * (u1 + u0)
    vm = 0.5 * (v1 + v0)
    return function (u::Float64, v::Float64)
        a = um + du * u
        b = vm + dv * v
        # Tangent-gnomonic cubed-sphere projection (C^∞ everywhere, no polar crowding)
        ta, tb = tan(a), tan(b)
        sa2, sb2 = sec(a)^2, sec(b)^2
        invr = 1.0 / sqrt(1.0 + ta^2 + tb^2)
        if face == 1      # +x
            q = SVector(1.0, ta, tb); qa = SVector(0.0, sa2, 0.0); qb = SVector(0.0, 0.0, sb2)
        elseif face == 2  # -x
            q = SVector(-1.0, -ta, tb); qa = SVector(0.0, -sa2, 0.0); qb = SVector(0.0, 0.0, sb2)
        elseif face == 3  # +y
            q = SVector(-ta, 1.0, tb); qa = SVector(-sa2, 0.0, 0.0); qb = SVector(0.0, 0.0, sb2)
        elseif face == 4  # -y
            q = SVector(ta, -1.0, tb); qa = SVector(sa2, 0.0, 0.0); qb = SVector(0.0, 0.0, sb2)
        elseif face == 5  # +z
            q = SVector(-tb, ta, 1.0); qa = SVector(0.0, sa2, 0.0); qb = SVector(-sb2, 0.0, 0.0)
        else              # -z
            q = SVector(tb, ta, -1.0); qa = SVector(0.0, sa2, 0.0); qb = SVector(sb2, 0.0, 0.0)
        end
        s_hat = q * invr
        ds_da = (qa - s_hat * dot(s_hat, qa)) * invr
        ds_db = (qb - s_hat * dot(s_hat, qb)) * invr
        scale = SVector(rx, ry, rz)
        pos = center + scale .* s_hat
        ru = (scale .* ds_da) * du
        rv = (scale .* ds_db) * dv
        cr = cross(ru, rv)
        jac = norm(cr)
        nrm = cr / jac
        # Enforce outward normal orientation
        if dot(nrm, pos - center) < 0
            nrm = -nrm
        end
        return pos, nrm, jac
    end
end

function discretize_surface(s::Sphere3D; nsub::Int=2, order::Int=4, region_id::Int=1)
    patches = CurvedPatch3D[]
    edges = range(-pi / 4, pi / 4; length=nsub + 1)
    for face in 1:6, ix in 1:nsub, iy in 1:nsub
        mf = _cubed_sphere_patch_map(face, edges[ix], edges[ix+1], edges[iy], edges[iy+1],
                                     s.center, s.R, s.R, s.R)
        pos0, _, jac0 = mf(0.0, 0.0)
        push!(patches, CurvedPatch3D(mf, pos0, 2.0 * sqrt(jac0), region_id))
    end
    return build_surface_mesh(patches, order)
end

function discretize_surface(s::ProlateSpheroid3D; nsub::Int=2, order::Int=4, region_id::Int=1)
    patches = CurvedPatch3D[]
    edges = range(-pi / 4, pi / 4; length=nsub + 1)
    for face in 1:6, ix in 1:nsub, iy in 1:nsub
        mf = _cubed_sphere_patch_map(face, edges[ix], edges[ix+1], edges[iy], edges[iy+1],
                                     s.center, s.a, s.a, s.c)
        pos0, _, jac0 = mf(0.0, 0.0)
        push!(patches, CurvedPatch3D(mf, pos0, 2.0 * sqrt(jac0), region_id))
    end
    return build_surface_mesh(patches, order)
end

"""
    discretize_surface(tsv::CurvilinearTSV3D; n_theta=8, n_z=6, order=4, capped=true, region_id=1)

Creates a watertight high-order curvilinear patch mesh of a 3D tapered/scalloped teardrop TSV
including sidewall and (optionally) top and bottom teardrop end-caps.
"""
function discretize_surface(tsv::CurvilinearTSV3D; n_theta::Int=8, n_z::Int=6, order::Int=4,
                            capped::Bool=true, region_id::Int=1)
    patches = CurvedPatch3D[]
    t_edges = range(0.0, 2pi; length=n_theta + 1)
    z_edges = range(tsv.zmin, tsv.zmax; length=n_z + 1)
    prof = tsv.profile
    c_xy = hasproperty(prof, :center) ? prof.center : SVector(0.0, 0.0)

    # 1. Sidewall patches
    for it in 1:n_theta, iz in 1:n_z
        t0, t1 = t_edges[it], t_edges[it+1]
        z0, z1 = z_edges[iz], z_edges[iz+1]
        tm, dt = 0.5 * (t0 + t1), 0.5 * (t1 - t0)
        zm, dz = 0.5 * (z0 + z1), 0.5 * (z1 - z0)
        mf = function (u::Float64, v::Float64)
            t = tm + dt * u
            z = zm + dz * v
            r2, dr2, _ = eval_jet(prof, t)
            rel = r2 - c_xy
            s, ds = scale_profile(tsv, z)
            pos = SVector(c_xy[1] + s * rel[1], c_xy[2] + s * rel[2], z)
            ru = SVector(s * dr2[1], s * dr2[2], 0.0) * dt
            rv = SVector(ds * rel[1], ds * rel[2], 1.0) * dz
            cr = cross(ru, rv)
            jac = norm(cr)
            nrm = cr / jac
            # Ensure outward radial orientation
            if nrm[1] * rel[1] + nrm[2] * rel[2] < 0
                nrm = -nrm
            end
            return pos, nrm, jac
        end
        pos0, _, jac0 = mf(0.0, 0.0)
        push!(patches, CurvedPatch3D(mf, pos0, 2.0 * sqrt(jac0), region_id))
    end

    # 2. Top & bottom end-cap patches (polar sectors mapped smoothly from center to boundary)
    if capped
        for (z_cap, nz_sign) in ((tsv.zmin, -1.0), (tsv.zmax, 1.0))
            s_cap, _ = scale_profile(tsv, z_cap)
            for it in 1:n_theta
                t0, t1 = t_edges[it], t_edges[it+1]
                tm, dt = 0.5 * (t0 + t1), 0.5 * (t1 - t0)
                mf = function (u::Float64, v::Float64)
                    t = tm + dt * u
                    ρ = 0.5 * (v + 1.0)   # radial fraction in [0, 1]
                    r2, dr2, _ = eval_jet(prof, t)
                    rel = s_cap * (r2 - c_xy)
                    drel = s_cap * dr2
                    pos = SVector(c_xy[1] + ρ * rel[1], c_xy[2] + ρ * rel[2], z_cap)
                    ru = SVector(ρ * drel[1], ρ * drel[2], 0.0) * dt
                    rv = SVector(0.5 * rel[1], 0.5 * rel[2], 0.0)
                    jac = abs(ru[1] * rv[2] - ru[2] * rv[1])
                    nrm = SVector(0.0, 0.0, nz_sign)
                    return pos, nrm, jac
                end
                pos0, _, jac0 = mf(0.0, 0.0)
                push!(patches, CurvedPatch3D(mf, pos0, 2.0 * sqrt(max(jac0, 1e-30)), region_id))
            end
        end
    end
    return build_surface_mesh(patches, order)
end

surface_area(m::SurfaceMesh3D) = sum(m.weights)
surface_volume(m::SurfaceMesh3D) = sum(dot(m.pts[i], m.normals[i]) * m.weights[i] for i in eachindex(m.pts)) / 3.0

"""
    load_gmsh_msh(filepath::AbstractString) -> (nodes, triangles, physical_tags)

Lightweight reader for Gmsh v2.2 ASCII `.msh` surface/volume meshes for interoperability
with existing `.geo` workflows in this repository.
"""
function load_gmsh_msh(filepath::AbstractString)
    nodes = Vector{SVector{3,Float64}}()
    tris = Vector{SVector{3,Int}}()
    tags = Vector{Int}()
    open(filepath, "r") do io
        while !eof(io)
            line = strip(readline(io))
            if line == "\$Nodes"
                n = parse(Int, strip(readline(io)))
                resize!(nodes, n)
                for _ in 1:n
                    parts = split(strip(readline(io)))
                    id = parse(Int, parts[1])
                    nodes[id] = SVector(parse(Float64, parts[2]),
                                        parse(Float64, parts[3]),
                                        parse(Float64, parts[4]))
                end
            elseif line == "\$Elements"
                ne = parse(Int, strip(readline(io)))
                for _ in 1:ne
                    parts = split(strip(readline(io)))
                    etype = parse(Int, parts[2])
                    ntags = parse(Int, parts[3])
                    if etype == 2 # 3-node triangle
                        phys = ntags >= 1 ? parse(Int, parts[4]) : 0
                        base = 3 + ntags
                        push!(tris, SVector(parse(Int, parts[base+1]),
                                            parse(Int, parts[base+2]),
                                            parse(Int, parts[base+3])))
                        push!(tags, phys)
                    end
                end
            end
        end
    end
    return nodes, tris, tags
end
