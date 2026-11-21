(* ::Package:: *)
(* File: 05_sparameters_hbm.wl *)
(* Symbolic verification of MTL ABCD <-> Y <-> Z <-> S transformations, reciprocity, passivity, *)
(* and 1-Hi / 8-Hi HBM 60 GHz TSV S-parameter reference benchmarks (inductive vs purely dissipative). *)

VerifySParametersHBM[] := Module[
  {abcdFromRLGC, sFromABCD, recipCheck, unitCheckLossless,
   mu0, eps0, sigCu, sigTa, epsOx, epsSi, sigSi, rCu, tBar, rOut, tox, rox, pPitch, z0,
   zCladNum, computeS, freqs, s1HiInd, s8HiInd, s8HiDiss},

  (* Symbolic 2-port transmission line ABCD matrix for length h *)
  abcdFromRLGC[Z_, Y_, h_] := Module[{gam = Sqrt[Z Y], Zc = Sqrt[Z/Y]},
    {{Cosh[gam h], Zc Sinh[gam h]}, {Sinh[gam h]/Zc, Cosh[gam h]}}
  ];

  (* Standard ABCD -> S transformation for port reference impedance Z0 *)
  sFromABCD[{{A_, B_}, {C_, D_}}, Z0_] := Module[{den = A + B/Z0 + C Z0 + D},
    {
      {(A + B/Z0 - C Z0 - D)/den, 2 (A D - B C)/den},
      {2/den, (-A + B/Z0 - C Z0 + D)/den}
    }
  ];

  (* Verify exact reciprocity (S12 == S21) and unitarity (S^H S == I) for lossless line *)
  With[{M = abcdFromRLGC[I w L, I w C0, h]},
    Module[{S = sFromABCD[M, Z0]},
      recipCheck = FullSimplify[S[[1, 2]] - S[[2, 1]], w > 0 && L > 0 && C0 > 0 && h > 0 && Z0 > 0];
      unitCheckLossless = FullSimplify[
        ConjugateTranspose[S].S - IdentityMatrix[2],
        w > 0 && L > 0 && C0 > 0 && h > 0 && Z0 > 0
      ];
    ]
  ];

  (* Numerical HBM benchmarks: 1-Hi (50 um) and 8-Hi (400 um) twin-TSV channel *)
  mu0 = 4 Pi * 10^-7;
  eps0 = 8.8541878128*^-12;
  sigCu = 5.8*^7;
  sigTa = 7.6*^6;
  epsOx = 3.9 eps0;
  epsSi = 11.7 eps0;
  sigSi = 10.0;
  rCu = 2.5*^-6;
  tBar = 30.0*^-9;
  rOut = rCu + tBar;
  tox = 0.2*^-6;
  rox = rOut + tox;
  pPitch = 20.0*^-6;
  z0 = 50.0;

  zCladNum[w_] := Module[{g1 = Sqrt[I w mu0 sigCu], g2 = Sqrt[I w mu0 sigTa], bA, cA, num, den},
    bA = g2 rCu (BesselI[0, g1 rCu] BesselK[1, g2 rCu] + (g1/g2) BesselI[1, g1 rCu] BesselK[0, g2 rCu]);
    cA = g2 rCu (BesselI[0, g1 rCu] BesselI[1, g2 rCu] - (g1/g2) BesselI[1, g1 rCu] BesselI[0, g2 rCu]);
    num = bA BesselI[0, g2 rOut] + cA BesselK[0, g2 rOut];
    den = bA BesselI[1, g2 rOut] - cA BesselK[1, g2 rOut];
    (g2/(2 Pi rOut sigTa)) * (num/den)
  ];

  computeS[f_, h_, rdep_, inductiveQ_] := Module[
    {w = 2 Pi f, zInt, lExt, Z, cOx, cDep, cIns, ySi, Y, M, S, svdMax},
    zInt = 2 N[zCladNum[w], 18];
    lExt = If[inductiveQ, (mu0/Pi) ArcCosh[pPitch/(2 rOut)], 0.0];
    Z = If[inductiveQ, zInt + I w lExt, Re[zInt]];
    cOx = (2 Pi epsOx)/Log[rox/rOut];
    cDep = (2 Pi epsSi)/Log[rdep/rox];
    cIns = 0.5 / (1/cOx + 1/cDep);
    ySi = (Pi (sigSi + I w epsSi))/ArcCosh[pPitch/(2 rdep)];
    Y = 1/(1/(I w cIns) + 1/ySi);
    M = abcdFromRLGC[Z, Y, h];
    S = N[sFromABCD[M, z0], 16];
    svdMax = Max[SingularValueList[S]];
    <|
      "f_Hz" -> f,
      "S11_re" -> Re[S[[1, 1]]], "S11_im" -> Im[S[[1, 1]]],
      "S21_re" -> Re[S[[2, 1]]], "S21_im" -> Im[S[[2, 1]]],
      "S21_dB" -> 20 Log10[Abs[S[[2, 1]]]],
      "S11_dB" -> 20 Log10[Abs[S[[1, 1]]]],
      "passive_sigma_max" -> svdMax
    |>
  ];

  freqs = {1.0*^9, 10.0*^9, 30.0*^9, 60.0*^9};
  s1HiInd = Table[computeS[f, 50.0*^-6, 3.4*^-6, True], {f, freqs}];
  s8HiInd = Table[computeS[f, 400.0*^-6, 3.4*^-6, True], {f, freqs}];
  s8HiDiss = Table[computeS[f, 400.0*^-6, 3.4*^-6, False], {f, freqs}];

  <|
    "Reciprocity_zero" -> (recipCheck === 0),
    "LosslessUnitarity_zero" -> (unitCheckLossless === {{0, 0}, {0, 0}}),
    "Benchmarks" -> <|
      "S_1Hi_inductive" -> s1HiInd,
      "S_8Hi_inductive" -> s8HiInd,
      "S_8Hi_dissipative" -> s8HiDiss
    |>
  |>
];
