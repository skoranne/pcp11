# Defensible MTL/RLC extraction for simple2 2-port TSV model (same geometry/materials as palace_tsv_simple2)
# Not full-wave 3D; uses per-unit-length TL + lumped via approximation. Reports passivity/reciprocity.
using LinearAlgebra, Printf

eps0=8.8541878128e-12; mu0=4*pi*1e-7; c0=1/sqrt(eps0*mu0)
sigCu=5.8e7; epsr_si=11.7; sig_si=10.0; epsr_ox=3.9

# Geometry (m)
h_si=50e-6; t_ox=0.5e-6; tsv_r=5e-6; tsv_h=h_si
trace_w=6e-6; trace_t=2e-6; gnd_t=2e-6
# Port-to-via trace lengths (approx, from x=±60 to via at x=0, pad 12um)
L_top = 60e-6 - 6e-6  # 54um
L_bot = 60e-6 - 6e-6

f = 1e9
omega = 2*pi*f

# Microstrip-like per-length on Si (simplified effective medium)
# top trace over Si substrate with air above
h = h_si
eps_eff = (epsr_si+1)/2 + (epsr_si-1)/2 / sqrt(1+12*h/trace_w)
Z0_ms = if trace_w/h <= 1
    60/sqrt(eps_eff)*log(8*h/trace_w + trace_w/(4*h))
else
    120*pi/(sqrt(eps_eff)*(trace_w/h + 1.393 + 0.667*log(trace_w/h+1.444)))
end
v_ms = c0/sqrt(eps_eff)
Cpl = 1/(Z0_ms*v_ms)
Lpl = Z0_ms/v_ms
# conductor loss (skin effect) + substrate loss (simplified)
delta_cu = sqrt(2/(omega*mu0*sigCu))
Rs = 1/(sigCu*delta_cu)
Rpl = Rs*(1/trace_w + 1/(2*h)) # rough
Gpl = omega*Cpl*epsr_si/(eps_eff)* (sig_si/(omega*eps0*epsr_si)) * 0.5 # rough substrate loss
# better: G = omega*C*tand_eff; use q factor
q = (eps_eff-1)/(epsr_si-1)
Gpl2 = omega*eps0*epsr_si*q * (sig_si/(omega*eps0*epsr_si))

function tl_abcd(Rp,Lp,Gp,Cp,om,len)
    gamma = sqrt((Rp+im*om*Lp)*(Gp+im*om*Cp))
    Zc = sqrt((Rp+im*om*Lp)/(Gp+im*om*Cp))
    ch = cosh(gamma*len); sh = sinh(gamma*len)
    [ch Zc*sh; sh/Zc ch]
end

A_top = tl_abcd(Rpl,Lpl,Gpl2,Cpl,omega,L_top)
A_bot = tl_abcd(Rpl,Lpl,Gpl2,Cpl,omega,L_bot)

# TSV via lumped model: Rvia + jωLvia in series, Cvia to ground (oxide)
Rvia = tsv_h/(sigCu*pi*tsv_r^2)
Lvia = mu0*tsv_h/(2*pi)*(log(2*tsv_h/tsv_r)-0.75)
Cox = 2*pi*eps0*epsr_ox*tsv_h/log((tsv_r+t_ox)/tsv_r)
# simple pi: C/2 - series RL - C/2
Zs_via = Rvia + im*omega*Lvia
Ysh = im*omega*Cox/2
A_via = [1+Zs_via*Ysh Zs_via; 2*Ysh+Zs_via*Ysh^2 1+Zs_via*Ysh]

A_total = A_top * A_via * A_bot

# ABCD to S
Z0ref=50.0
A,B,Cc,D = A_total[1,1],A_total[1,2],A_total[2,1],A_total[2,2]
den = A + B/Z0ref + Cc*Z0ref + D
S11 = (A + B/Z0ref - Cc*Z0ref - D)/den
S21 = 2/den
S12 = 2*(A*D-B*Cc)/den
S22 = (-A + B/Z0ref - Cc*Z0ref + D)/den
S = [S11 S12; S21 S22]

println("=== MTL/RLC extraction (simple2, 1 GHz) ===")
@printf("eps_eff=%.4f Z0_ms=%.2f Ohm\n", eps_eff, Z0_ms)
@printf("Rpl=%.3e Ohm/m Lpl=%.3e H/m Gpl=%.3e S/m Cpl=%.3e F/m\n", Rpl,Lpl,Gpl2,Cpl)
@printf("Rvia=%.3e Ohm Lvia=%.3e H Cox=%.3e F\n", Rvia,Lvia,Cox)
@printf("S11=%.6f dB ang=%.3f deg\n", 20*log10(abs(S11)), rad2deg(angle(S11)))
@printf("S21=%.6f dB ang=%.3f deg\n", 20*log10(abs(S21)), rad2deg(angle(S21)))
@printf("S12=%.6f dB S22=%.6f dB\n", 20*log10(abs(S12)), rad2deg(angle(S22)))
println("max|S-S^T|=", maximum(abs.(S-transpose(S))))
println("sigma_max(S)=", maximum(svdvals(S)))
println("passive=", maximum(svdvals(S)) <= 1+1e-9)
