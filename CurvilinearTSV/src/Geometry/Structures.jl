# Physical constants, material stacks, cladded TSV structures, and HBM multi-TSV array builders

const EPS0 = 8.8541878128e-12      # Vacuum permittivity [F/m]
const MU0  = 4.0 * pi * 1e-7       # Vacuum permeability [H/m]
const Q_E  = 1.602176634e-19       # Elementary charge [C]
const K_B  = 1.380649e-23          # Boltzmann constant [J/K]
const NI_SI_300K = 1.0e16          # Intrinsic carrier density in Si at 300 K [m^-3]

"""
    TSVMaterialStack

Complete material stack for a modern Cu/Ta-barrier/SiO2/p-Si TSV in HBM 3D-IC technology.
"""
Base.@kwdef struct TSVMaterialStack
    sigma_cu::Float64 = 5.8e7          # Cu core conductivity [S/m]
    sigma_barrier::Float64 = 7.6e6     # Ta/TiN cladding liner conductivity [S/m]
    t_barrier::Float64 = 30.0e-9       # Barrier cladding thickness [m]
    eps_r_ox::Float64 = 3.9            # SiO2 liner relative permittivity
    t_ox::Float64 = 0.2e-6             # SiO2 liner thickness [m]
    eps_r_si::Float64 = 11.7           # Silicon substrate relative permittivity
    sigma_si::Float64 = 10.0           # Bulk p-Si conductivity [S/m] (~10 Ohm-cm)
    N_A::Float64 = 1.0e21              # Acceptor doping [m^-3] (= 1e15 cm^-3)
    V_FB::Float64 = -0.35              # Flatband voltage [V]
    T_kelvin::Float64 = 300.0          # Temperature [K]
    eps_r_ild::Float64 = 3.5           # Inter-layer dielectric / underfill relative permittivity
    Q_f::Float64 = 0.0                 # Fixed oxide interface charge density [C/m^2]
end

eps_ox(m::TSVMaterialStack) = m.eps_r_ox * EPS0
eps_si(m::TSVMaterialStack) = m.eps_r_si * EPS0
eps_ild(m::TSVMaterialStack) = m.eps_r_ild * EPS0
thermal_voltage(m::TSVMaterialStack) = K_B * m.T_kelvin / Q_E
fermi_potential(m::TSVMaterialStack) = thermal_voltage(m) * log(m.N_A / NI_SI_300K)
complex_eps_si(m::TSVMaterialStack, f_hz::Float64) =
    ComplexF64(eps_si(m), -m.sigma_si / (2pi * f_hz))

"""
    CladdedTSV2D

2D cross-section of a single curvilinear Cu TSV with barrier cladding and SiO2 liner.
"""
struct CladdedTSV2D{C<:AbstractCurve2D}
    core::C
    barrier_outer::OffsetCurve2D{C}
    oxide_outer::OffsetCurve2D{C}
    materials::TSVMaterialStack
    role::Symbol   # :signal or :ground
end

function CladdedTSV2D(core::AbstractCurve2D, mat::TSVMaterialStack=TSVMaterialStack(); role::Symbol=:signal)
    ok_b, mf_b = verify_offset_regularity(core, mat.t_barrier)
    ok_o, mf_o = verify_offset_regularity(core, mat.t_barrier + mat.t_ox)
    (ok_b && ok_o) || throw(ArgumentError(
        "Offset liner self-intersects (min regularity factor = $(min(mf_b, mf_o))); reduce alpha or t_ox"))
    return CladdedTSV2D(
        core,
        OffsetCurve2D(core, mat.t_barrier),
        OffsetCurve2D(core, mat.t_barrier + mat.t_ox),
        mat,
        role
    )
end

"""
    HBMTSVArray2D

Multi-TSV 2D cross-sectional array (e.g., Twin SG, GSG, Differential GSSG, 4×4 HBM bundle)
supporting per-TSV teardrop stress/gas-escape orientations.
"""
struct HBMTSVArray2D
    tsvs::Vector{CladdedTSV2D}
    materials::TSVMaterialStack
    die_height::Float64
    num_dies::Int
    bump_capacitance::Float64     # 3D micro-bump + landing pad shunt C per die transition [F]
    bump_inductance::Float64      # 3D micro-bump + RDL series L per die transition [H]
    bump_resistance::Float64      # 3D micro-bump contact resistance [Ohm]
end

"""
    build_hbm_gsg_array(; R0=2.5e-6, pitch=20e-6, alpha=0.35, theta0=0.0, num_dies=1, mat=TSVMaterialStack())

Builds a Ground-Signal-Ground (GSG) curvilinear teardrop HBM TSV array.
"""
function build_hbm_gsg_array(; R0::Float64=2.5e-6, pitch::Float64=20.0e-6,
                             alpha::Float64=0.35, cy::Float64=1.0, theta0::Float64=0.0,
                             die_height::Float64=50.0e-6, num_dies::Int=1,
                             mat::TSVMaterialStack=TSVMaterialStack(),
                             bump_capacitance::Float64=6.0e-15,
                             bump_inductance::Float64=8.0e-12,
                             bump_resistance::Float64=0.04)
    c_g1 = Teardrop2D(R0; cy=cy, alpha=alpha, theta0=theta0 + pi, center=(-pitch, 0.0))
    c_s  = Teardrop2D(R0; cy=cy, alpha=alpha, theta0=theta0,      center=(0.0, 0.0))
    c_g2 = Teardrop2D(R0; cy=cy, alpha=alpha, theta0=theta0,      center=(pitch, 0.0))
    tsvs = CladdedTSV2D[
        CladdedTSV2D(c_g1, mat; role=:ground),
        CladdedTSV2D(c_s,  mat; role=:signal),
        CladdedTSV2D(c_g2, mat; role=:ground)
    ]
    return HBMTSVArray2D(tsvs, mat, die_height, num_dies, bump_capacitance, bump_inductance, bump_resistance)
end

"""
    build_hbm_differential_array(; R0=2.5e-6, pitch=20e-6, alpha=0.35, ...)

Builds a G-S+-S--G differential pair HBM teardrop TSV channel (for NEXT/FEXT & differential S-params).
"""
function build_hbm_differential_array(; R0::Float64=2.5e-6, pitch::Float64=20.0e-6,
                                      alpha::Float64=0.35, cy::Float64=1.0,
                                      die_height::Float64=50.0e-6, num_dies::Int=1,
                                      mat::TSVMaterialStack=TSVMaterialStack(),
                                      bump_capacitance::Float64=6.0e-15,
                                      bump_inductance::Float64=8.0e-12,
                                      bump_resistance::Float64=0.04)
    xs = (-1.5 * pitch, -0.5 * pitch, 0.5 * pitch, 1.5 * pitch)
    roles = (:ground, :signal, :signal, :ground)
    angles = (pi, pi / 4, -pi / 4, 0.0) # Outward radial gas-escape / stress-relief orientation
    tsvs = CladdedTSV2D[
        CladdedTSV2D(Teardrop2D(R0; cy=cy, alpha=alpha, theta0=angles[k], center=(xs[k], 0.0)), mat; role=roles[k])
        for k in 1:4
    ]
    return HBMTSVArray2D(tsvs, mat, die_height, num_dies, bump_capacitance, bump_inductance, bump_resistance)
end

"""
    build_hbm_4x4_staggered_array(; R0=2.5e-6, pitch=20e-6, alpha=0.35, ...)

Builds a 4×4 HBM TSV array where teardrop apexes point radially outward from array center
to relieve 3D thermo-mechanical CTE shear stress and provide outward gas-escape pathways.
"""
function build_hbm_4x4_staggered_array(; R0::Float64=2.5e-6, pitch::Float64=20.0e-6,
                                       alpha::Float64=0.35, cy::Float64=1.0,
                                       die_height::Float64=50.0e-6, num_dies::Int=1,
                                       mat::TSVMaterialStack=TSVMaterialStack(),
                                       bump_capacitance::Float64=6.0e-15,
                                       bump_inductance::Float64=8.0e-12,
                                       bump_resistance::Float64=0.04)
    tsvs = CladdedTSV2D[]
    for ix in 0:3, iy in 0:3
        x = (ix - 1.5) * pitch + (isodd(iy) ? 0.25 * pitch : 0.0)
        y = (iy - 1.5) * pitch
        θ = atan(y, x) # Outward radial stress/outgassing alignment
        role = iseven(ix + iy) ? :ground : :signal
        push!(tsvs, CladdedTSV2D(Teardrop2D(R0; cy=cy, alpha=alpha, theta0=θ, center=(x, y)), mat; role=role))
    end
    return HBMTSVArray2D(tsvs, mat, die_height, num_dies, bump_capacitance, bump_inductance, bump_resistance)
end
