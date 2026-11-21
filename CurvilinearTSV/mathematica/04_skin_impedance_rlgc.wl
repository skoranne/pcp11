(* ::Package:: *)
(* File: 04_skin_impedance_rlgc.wl *)
(* Exact multi-layer cylindrical Cu + Ta/TiN barrier cladding internal surface impedance Z_int(omega), *)
(* slow-wave complex Si permittivity, and 0.1 - 60 GHz analytical coaxial/twin-TSV RLGC benchmarks. *)

VerifySkinAndRLGC[] := Module[
  {zSolid, zSolidDC, zSolidHF, zClad, zCladDCcheck,
   mu0, eps0, sigCu, sigTa, epsOx, epsSi, sigSi, rCu, tBar, rOut, tox, rox, pPitch,
   freqs, skinTable, rlgcTable},

  (* 1. Solid cylindrical conductor exact Bessel internal impedance:
     Z_int(omega) = (gamma / (2 Pi r sigma)) * I0(gamma r) / I1(gamma r) *)
  zSolid[w_, r_, sig_, mu_] := Module[{g = Sqrt[I w mu sig]},
    (g/(2 Pi r sig)) * BesselI[0, g r]/BesselI[1, g r]
  ];
  zSolidDC = FullSimplify[Series[zSolid[w, r, sig, mu], {w, 0, 1}], r > 0 && sig > 0 && mu > 0];
  (* Note: R_DC = 1/(Pi r^2 sig), L_int_DC = mu / (8 Pi) *)

  (* 2. Exact Two-Layer Cladded Cylinder: Cu core (0 <= r <= r1, sig1) + Barrier liner (r1 <= r <= r2, sig2)
     In core: E1(r) = A I0(g1 r), H1(r) = (g1/(I w mu)) A I1(g1 r)
     In cladding: E2(r) = B I0(g2 r) + C K0(g2 r), H2(r) = (g2/(I w mu)) (B I1(g2 r) - C K1(g2 r))
     Matching E1(r1) = E2(r1) and H1(r1) = H2(r1) gives B/A and C/A in closed form! *)
  zClad[w_, r1_, r2_, sig1_, sig2_, mu_] := Module[
    {g1, g2, bA, cA, num, den},
    g1 = Sqrt[I w mu sig1];
    g2 = Sqrt[I w mu sig2];
    (* Wronskian I0(z) K1(z) + I1(z) K0(z) = 1/z *)
    bA = g2 r1 (BesselI[0, g1 r1] BesselK[1, g2 r1] + (g1/g2) BesselI[1, g1 r1] BesselK[0, g2 r1]);
    cA = g2 r1 (BesselI[0, g1 r1] BesselI[1, g2 r1] - (g1/g2) BesselI[1, g1 r1] BesselI[0, g2 r1]);
    num = bA BesselI[0, g2 r2] + cA BesselK[0, g2 r2];
    den = bA BesselI[1, g2 r2] - cA BesselK[1, g2 r2];
    (g2/(2 Pi r2 sig2)) * (num/den)
  ];

  (* Verify exact DC resistance of cladded cylinder: R_DC = 1 / (Pi (sig1 r1^2 + sig2 (r2^2 - r1^2))) *)
  mu0 = 4 Pi * 10^-7;
  eps0 = 8.8541878128*^-12;
  sigCu = 5.8*^7;        (* Cu conductivity [S/m] *)
  sigTa = 7.6*^6;        (* Ta barrier liner conductivity [S/m] *)
  epsOx = 3.9 eps0;
  epsSi = 11.7 eps0;
  sigSi = 10.0;          (* p-Si substrate conductivity [S/m] (~10 Ohm-cm) *)
  rCu = 2.5*^-6;
  tBar = 30.0*^-9;       (* 30 nm Ta barrier cladding *)
  rOut = rCu + tBar;
  tox = 0.2*^-6;
  rox = rOut + tox;
  pPitch = 20.0*^-6;

  With[{rDCexact = 1/(Pi (sigCu rCu^2 + sigTa (rOut^2 - rCu^2))),
        zLow = N[zClad[2 Pi * 1.0, rCu, rOut, sigCu, sigTa, mu0], 25]},
    zCladDCcheck = Abs[Re[zLow] - rDCexact]/rDCexact;
  ];

  freqs = {1.0*^8, 1.0*^9, 5.0*^9, 10.0*^9, 20.0*^9, 30.0*^9, 45.0*^9, 60.0*^9};
  skinTable = Table[
    Module[{w = 2 Pi f, zc, zs, dCu, lossTanSi},
      zc = N[zClad[w, rCu, rOut, sigCu, sigTa, mu0], 18];
      zs = N[zSolid[w, rOut, sigCu, mu0], 18];
      dCu = N[1/Sqrt[Pi f mu0 sigCu], 15];
      lossTanSi = N[sigSi/(w epsSi), 15];
      <|
        "f_Hz" -> f,
        "delta_Cu_m" -> dCu,
        "Si_loss_tangent" -> lossTanSi,
        "R_clad_Ohm_per_m" -> Re[zc],
        "L_int_clad_H_per_m" -> Im[zc]/w,
        "R_solidCu_Ohm_per_m" -> Re[zs],
        "L_int_solidCu_H_per_m" -> Im[zs]/w
      |>
    ],
    {f, freqs}
  ];

  (* Analytical Twin-TSV (Signal + Ground return at pitch p) RLGC at V_bias = 0 (rdep = 3.4 um sample) *)
  With[{rdepSample = 3.4*^-6},
    rlgcTable = Table[
      Module[{w = 2 Pi f, zInt, lExt, zTot, cOx, cDep, cIns, ySi, yTot},
        zInt = 2 * N[zClad[w, rCu, rOut, sigCu, sigTa, mu0], 18]; (* signal + return *)
        lExt = (mu0/Pi) ArcCosh[pPitch/(2 rOut)];
        zTot = zInt + I w lExt;
        (* Series insulator + depletion per conductor, half for two conductors in series *)
        cOx = (2 Pi epsOx)/Log[rox/rOut];
        cDep = (2 Pi epsSi)/Log[rdepSample/rox];
        cIns = 0.5 / (1/cOx + 1/cDep); (* signal + ground MOS capacitance in series *)
        (* Bulk lossy Si admittance between two cylinders of radius rdepSample separated by pPitch *)
        ySi = (Pi (sigSi + I w epsSi))/ArcCosh[pPitch/(2 rdepSample)];
        yTot = 1/(1/(I w cIns) + 1/ySi);
        <|
          "f_Hz" -> f,
          "R_Ohm_per_m" -> Re[zTot],
          "L_H_per_m" -> Im[zTot]/w,
          "G_S_per_m" -> Re[yTot],
          "C_F_per_m" -> Im[yTot]/w
        |>
      ],
      {f, freqs}
    ];
  ];

  <|
    "SolidDC_R_exact" -> FullSimplify[SeriesCoefficient[zSolidDC, 0] - 1/(Pi r^2 sig), r > 0 && sig > 0],
    "SolidDC_Lint_exact" -> FullSimplify[SeriesCoefficient[zSolidDC, 1]/(I) - mu/(8 Pi), r > 0 && sig > 0 && mu > 0],
    "CladDC_rel_err" -> N[zCladDCcheck, 15],
    "Benchmarks" -> <|
      "skin_impedance" -> skinTable,
      "twin_tsv_rlgc" -> rlgcTable
    |>
  |>
];
