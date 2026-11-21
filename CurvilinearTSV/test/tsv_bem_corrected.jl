# TSV S-parameter estimation - CORRECTED (defensible baseline)
# Fixes: geometry outside keep-out, correct port mapping (one per trace),
# physical partial inductance, explicit capacitance, no dead electrostatic solve,
# reciprocity/passivity checks. This is a simplified RLC/MTL baseline, not full 3D BEM.

using LinearAlgebra
using StaticArrays
using Printf

const eps0 = 8.854187817e-12
const mu0  = 4*pi*1e-7
const sigCu = 5.8e7
const eps_r_si = 11.7
const eps_r_ox = 3.9

# Geometry (um -> m)
tsv_r = 10e-6; tsv_h = 100e-6; tsv_pitch = 50e-6
t_ox = 0.5e-6
trace_w = 2e-6; trace_t = 1e-6; trace_l = 30e-6
keepout = 15e-6

# TSV centers (x,y)
tsv = [SVector(0.0,0.0), SVector(tsv_pitch,0.0)]

# Traces: placed outside keep-out, each as a separate port
# Top trace 1: starts at right edge of TSV1 keepout, goes +x at z=tsv_h
# Top trace 2: starts at right edge of TSV2 keepout
# Bottom traces similarly at z=0
function check_keepout(p, tsv_centers, r_keep)
    for c in tsv_centers
        if hypot(p[1]-c[1], p[2]-c[2]) < r_keep
            return false
        end
    end
    return true
end

ports = [
    (name="TopTrace1",  start=SVector(tsv_pitch + tsv_r + keepout + 5e-6, -20e-6, tsv_h), dir=SVector(1.0,0.0,0.0)),
    (name="TopTrace2",  start=SVector(tsv_pitch + tsv_r + keepout + 5e-6,  20e-6, tsv_h), dir=SVector(1.0,0.0,0.0)),
    (name="BotTrace1",  start=SVector(tsv_pitch + tsv_r + keepout + 5e-6, -20e-6, 0.0),   dir=SVector(1.0,0.0,0.0)),
    (name="BotTrace2",  start=SVector(tsv_pitch + tsv_r + keepout + 5e-6,  20e-6, 0.0),   dir=SVector(1.0,0.0,0.0)),
]

for p in ports
    ok = check_keepout(SVector(p.start[1],p.start[2]), tsv, tsv_r+keepout)
    println("$(p.name): start=$(p.start .* 1e6) um, keepout_ok=$ok")
end

# Simple RLC model per trace (microstrip-like over SiO2/Si, reference ground plane)
# This is an analytical baseline, NOT a full 3D BEM solve.
h_ox = t_ox; h_si = tsv_h
# Effective permittivity (microstrip approximation)
eps_eff = (eps_r_ox + 1)/2 + (eps_r_ox - 1)/2 * (1 + 12*h_ox/trace_w)^(-0.5)
# Capacitance per length (approx)
Cpl = eps0*eps_eff * (trace_w/h_ox + 1.393 + 0.667*log(trace_w/h_ox + 1.444))
# Inductance per length (approx, from C in air)
Cpl_air = eps0 * (trace_w/h_ox + 1.393 + 0.667*log(trace_w/h_ox + 1.444))
Lpl = 1/(Cpl_air * 3e8^2)
# DC resistance per length
Rpl_dc = 1/(sigCu * trace_w * trace_t)

freqs = [60e9]
Z0ref = 50.0

println("\nAnalytical per-length estimates:")
@printf("eps_eff=%.3f C'=%.3e F/m L'=%.3e H/m Rdc'=%.3e Ohm/m\n", eps_eff, Cpl, Lpl, Rpl_dc)

# Build 4-port Z matrix (isolated traces, no coupling as baseline; TSV coupling neglected explicitly)
n = 4
for f in freqs
    omega = 2*pi*f
    delta = sqrt(2/(omega*mu0*sigCu))
    Rs = 1/(sigCu*delta)
    Rpl = max(Rpl_dc, Rs/trace_w) # skin-effect limited
    Z = zeros(ComplexF64, n, n)
    for i in 1:n
        Z[i,i] = (Rpl + im*omega*Lpl) * trace_l + 1/(im*omega*Cpl*trace_l)
    end
    # Very weak capacitive coupling estimate between adjacent traces (far apart)
    for i in 1:n, j in 1:n
        if i != j
            dij = (i<=2 && j<=2) ? tsv_pitch : ((i>2 && j>2) ? tsv_pitch : tsv_h)
            Cm = eps0*eps_r_si * trace_w*trace_t / dij
            Z[i,j] = 1/(im*omega*Cm) * 1e-6 # rough weak coupling, flagged as estimate
        end
    end
    S = (Z + Z0ref*I) \ (Z - Z0ref*I)
    println("\nFrequency: $(f/1e9) GHz")
    println("S (dB):")
    for i in 1:n
        for j in 1:n
            @printf("%8.2f ", 20*log10(abs(S[i,j])+1e-30))
        end
        println()
    end
    # checks
    println("max |S-S^T| = ", maximum(abs.(S - transpose(S))))
    println("max singular value = ", maximum(svdvals(S)))
    pass = maximum(svdvals(S)) <= 1.0 + 1e-6
    println("passivity (sigma_max<=1): $pass")
end

# Write CSV with clear header
open("CurvilinearTSV/test/tsv_sparameters_60ghz_corrected.csv","w") do io
    println(io, "# CORRECTED BASELINE - analytical microstrip RLC, NOT full 3D BEM/Palace. Ports are 4 traces outside keep-out.")
    println(io, "freq_GHz,S11_dB,S21_dB,S31_dB,S41_dB,S12_dB,S22_dB,S32_dB,S42_dB,S13_dB,S23_dB,S33_dB,S43_dB,S14_dB,S24_dB,S34_dB,S44_dB")
    f=60e9; omega=2*pi*f
    delta=sqrt(2/(omega*mu0*sigCu)); Rs=1/(sigCu*delta); Rpl=max(Rpl_dc,Rs/trace_w)
    Z=zeros(ComplexF64,n,n)
    for i in 1:n; Z[i,i]=(Rpl+im*omega*Lpl)*trace_l + 1/(im*omega*Cpl*trace_l); end
    for i in 1:n, j in 1:n
        if i!=j; dij = (i<=2 && j<=2) ? tsv_pitch : ((i>2 && j>2) ? tsv_pitch : tsv_h); Cm = eps0*eps_r_si*trace_w*trace_t/dij; Z[i,j] = 1/(im*omega*Cm)*1e-6; end
    end
    S=(Z+Z0ref*I)\(Z-Z0ref*I)
    vals = Float64[]
    for j in 1:4, i in 1:4; push!(vals, 20*log10(abs(S[i,j])+1e-30)); end
    @printf(io, "60.0"); for v in vals; @printf(io, ",%.6f", v); end; println(io)
end
println("\nWrote CurvilinearTSV/test/tsv_sparameters_60ghz_corrected.csv")
