# 2D quasi-static BEM (constant panels, Laplace/ln Green) for per-unit-length C/L extraction.
# Defensible as 2D BEM; not 3D full-wave. Same 4-trace geometry as tsv_bem_corrected.jl.
using LinearAlgebra, Printf, Statistics

const eps0 = 8.8541878128e-12
const mu0  = 4*pi*1e-7

# Geometry (m) - same as corrected baseline
trace_w = 2e-6; trace_t = 1e-6; h_ox = 0.5e-6; h = h_ox + trace_t/2  # trace center above ground
keepout = 15e-6; tsv_r = 10e-6; tsv_pitch = 50e-6
# 4 traces: top1,top2,bot1,bot2 mapped to 4 conductors at same z=h for cross-section baseline
# Place traces outside keep-out
x0 = tsv_pitch + tsv_r + keepout + 5e-6
ys = [-20e-6, 20e-6, -20e-6+1e-3, 20e-6+1e-3] # separate for clarity; use 4 traces
# Use the 4 ports from corrected script: 2 top (y=-20,20), 2 bottom (y=-20,20) but bottom at different z.
# For 2D cross-section BEM, model all 4 as conductors above ground at their respective z.
traces = [
    (x=x0, y=-20e-6, z=h, w=trace_w, t=trace_t, name="T1"),
    (x=x0, y= 20e-6, z=h, w=trace_w, t=trace_t, name="T2"),
    (x=x0, y=-20e-6, z=h+100e-6, w=trace_w, t=trace_t, name="B1_offplane"), # different plane - not same 2D section
]

# Better: 2D BEM for top 2 traces over ground at z=0, and separately note coupling limits.
# Defensible 2D BEM extraction for the 2-trace cross-section.
traces2 = [
    (x=x0, y=-20e-6, z=h_ox + trace_t/2, w=trace_w, t=trace_t, name="P1"),
    (x=x0, y= 20e-6, z=h_ox + trace_t/2, w=trace_w, t=trace_t, name="P2"),
]

function rect_panels(tr; npw=6, npt=3)
    # panels around rectangle centered at (y,z), width w, thickness t. Return (p1,p2,center,normal,len)
    y0,z0,w,t = tr.y,tr.z,tr.w,tr.t
    corners = [(y0-w/2,z0-t/2),(y0+w/2,z0-t/2),(y0+w/2,z0+t/2),(y0-w/2,z0+t/2)]
    segs = [(corners[1],corners[2],npw),(corners[2],corners[3],npt),(corners[3],corners[4],npw),(corners[4],corners[1],npt)]
    panels = NamedTuple[]
    for (a,b,n) in segs
        for k in 1:n
            t0=(k-1)/n; t1=k/n
            p1=(a[1]+t0*(b[1]-a[1]), a[2]+t0*(b[2]-a[2]))
            p2=(a[1]+t1*(b[1]-a[1]), a[2]+t1*(b[2]-a[2]))
            c=((p1[1]+p2[1])/2,(p1[2]+p2[2])/2)
            dy=p2[1]-p1[1]; dz=p2[2]-p1[2]; L=hypot(dy,dz)
            nrm=(dz/L, -dy/L) # outward for CCW? check
            push!(panels,(p1=p1,p2=p2,c=c,n=nrm,L=L))
        end
    end
    panels
end

# 2D Green: G(r,r') = -1/(2*pi*eps) * ln|r-r'|, potential from line charge density sigma: V = -1/(2*pi*eps)*∫ sigma ln|r-r'| dl'
# With ground plane at z=0, use image: G = -1/(2*pi*eps)*(ln|r-r'| - ln|r-r'_im|)
function pot_coeff(p, q, eps)
    # constant element: integrate ln|r-r'| over source panel q, evaluated at center of p
    # simple midpoint approximation for off-diagonal; analytical for self
    cy,cz = p.c
    qy,qz = q.c
    ry = cy-qy; rz = cz-qz
    r = hypot(ry,rz)
    # image of q
    qyi,qzi = qy, -qz
    riy = cy-qyi; riz = cz-qzi
    ri = hypot(riy,riz)
    if r < 1e-15
        # self term: -1/(2*pi*eps)*∫_{-L/2}^{L/2} ln|s| ds = -1/(2*pi*eps)*L*(ln(L/2)-1)
        val = -q.L/(2*pi*eps)*(log(q.L/2)-1)
    else
        val = -q.L/(2*pi*eps)*log(r)
    end
    # image contribution (grounded at z=0)
    if ri < 1e-15
        val_im = -q.L/(2*pi*eps)*(log(q.L/2)-1)
    else
        val_im = -q.L/(2*pi*eps)*log(ri)
    end
    return val - val_im
end

function bem_capacitance(traces; epsr=1.0, npw=6, npt=3)
    eps = eps0*epsr
    all_panels = NamedTuple[]
    owner = Int[]
    for (i,tr) in enumerate(traces)
        ps = rect_panels(tr; npw=npw, npt=npt)
        for p in ps; push!(all_panels,p); push!(owner,i); end
    end
    N = length(all_panels)
    m = length(traces)
    Pmat = zeros(N,N)
    for i in 1:N, j in 1:N
        Pmat[i,j] = pot_coeff(all_panels[i], all_panels[j], eps)
    end
    # Solve for each conductor excited to 1V, others 0V
    C = zeros(m,m)
    for k in 1:m
        rhs = [owner[i]==k ? 1.0 : 0.0 for i in 1:N]
        # Least squares / collocation: P*sigma = V
        sigma = Pmat \ rhs
        for j in 1:m
            Q = 0.0
            for i in 1:N
                if owner[i]==j; Q += sigma[i]*all_panels[i].L; end
            end
            C[j,k] = Q
        end
    end
    # Symmetrize
    C = 0.5*(C+C')
    return C
end

C0 = bem_capacitance(traces2; epsr=1.0, npw=8, npt=4)
Cox = bem_capacitance(traces2; epsr=3.9, npw=8, npt=4)
# L from C0: L*C0 = mu0*eps0*I
L = mu0*eps0*inv(C0)

println("=== 2D BEM (Laplace, constant panels) ===")
println("C0 (F/m):"); display(C0)
println("Cox (F/m):"); display(Cox)
println("L (H/m):"); display(L)
println("C0 symmetric err: ", maximum(abs.(C0-C0')))
println("L symmetric err: ", maximum(abs.(L-L')))

# 60 GHz MTL ABCD
f=60e9; omega=2*pi*f
sigCu=5.8e7; delta=sqrt(2/(omega*mu0*sigCu)); Rs=1/(sigCu*delta)
# R per length approx
Rpl = Rs*(2/(trace_w)+2/(trace_t)) # rough all sides
R = diagm(fill(Rpl,2))
G = zeros(2,2) # lossless dielectric baseline
Z0ref=50.0; l=30e-6

function mtl_s(Z0ref, R,L,G,C,om,len)
    n=size(R,1)
    Zs = R + im*om*L
    Ys = G + im*om*C
    # eigen decomposition for MTL
    Gamma2 = Zs*Ys
    ev = eigen(Gamma2)
    gam = sqrt.(ev.values)
    # simple for 2-conductor symmetric; use direct ABCD via matrix exp
    # ABCD = expm([0 Zs; Ys 0]*len)
    M = [zeros(ComplexF64,n,n) Zs; Ys zeros(ComplexF64,n,n)] * len
    E = exp(M)
    A=E[1:n,1:n]; Bm=E[1:n,n+1:2n]; Cm2=E[n+1:2n,1:n]; Dm=E[n+1:2n,n+1:2n]
    I2=Matrix{ComplexF64}(I,n,n)
    den = (A + Bm/Z0ref + Cm2*Z0ref + Dm)
    S11 = den \ (A + Bm/Z0ref - Cm2*Z0ref - Dm)
    S21 = den \ (2I2)
    return S11,S21
end

S11,S21 = mtl_s(Z0ref,R,L,zeros(2,2),Cox,omega,l)
println("\n60 GHz 2-port MTL (from 2D BEM C/L):")
@printf("S11 dB: %.3f %.3f\n", 20*log10(abs(S11[1,1])), 20*log10(abs(S11[2,2])))
@printf("S21 dB: %.3f %.3f\n", 20*log10(abs(S21[1,1])), 20*log10(abs(S21[2,1])))
println("S21 full:"); display(20 .*log10.(abs.(S21) .+ 1e-30))
