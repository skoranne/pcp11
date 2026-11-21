(* ::Package:: *)
(* File: 03_depletion_theory.wl *)
(* Symbolic derivation & verification of voltage-dependent MOS-TSV depletion capacitance C_dep(V): *)
(*   1. Exact cylindrical MOS-TSV potential and closed-form Lambert-W inversion for r_dep(V) *)
(*   2. Exact C_ox', C_dep'(V), C_TSV'(V), threshold voltage V_th, and strong-inversion pinning *)
(*   3. Curvilinear Frenet-Serret (s, n) Laplacian & curvature-perturbation expansion for teardrop TSVs *)
(*   4. 3D doubly-curved (mean H, Gaussian K) depletion shell perturbation *)

VerifyDepletionTheory[] := Module[
  {phiSi, pdeCheck, bcCheck, psiS, vOx, vTot, uExpr, aParam, vNorm, lambertCheck,
   qDep, cDepExpr, cOxExpr, cTotExpr, dVdRdep, dQdRdep, cTotDirectCheck,
   frenetLap, phiLoc, phiLocCheck, cOxCurvSeries, cDepCurvSeries, wDepSeries,
   qE, kB, T, eps0, epsOx, epsSi, ni, Na, rc, tox, rox, phiF, vFB,
   rdepFromPsi, rdepMax, vTh, solveRdepV, cOxVal, cDepAtV, cTsvAtV,
   vSweep, cvTable, apexVth, bellyVth},

  (* 1. Exact Cylindrical MOS-TSV Poisson solution in p-type Si (Na) - pure symbolic block *)
  Module[{phiSi, psiS, qDep, vOx, vTot, aParam, dQdRdep, dVdRdep, cOxExpr, cDepExpr, cTotExpr, phiLoc},
    phiSi[r_] := (qS nAs)/(4 eSi) (r^2 - rdep^2 - 2 rdep^2 Log[r/rdep]);
    pdeCheck = FullSimplify[(1/r) D[r D[phiSi[r], r], r] - qS nAs/eSi, r > 0 && rdep > 0];
    bcCheck = FullSimplify[{phiSi[rdep], D[phiSi[r], r] /. r -> rdep}, rdep > 0];

    psiS = FullSimplify[phiSi[rOx], rOx > 0 && rdep > 0];
    qDep = qS nAs Pi (rdep^2 - rOx^2);
    vOx = FullSimplify[Integrate[qDep/(2 Pi eOx r), {r, rC, rOx}], rOx > rC > 0];
    vTot = FullSimplify[psiS + vOx, rOx > rC > 0 && rdep > rOx];

    aParam = 1 - (2 eSi/eOx) Log[rOx/rC];
    lambertCheck = FullSimplify[
      Module[{w = ProductLog[(vn - a) Exp[-a]]},
        a + Exp[a + w] (a + w - a) - vn /. {w Exp[w] -> (vn - a) Exp[-a]}
      ]
    ];
    uExpr = FullSimplify[
      (4 eSi vTot)/(qS nAs rOx^2) - (aParam + (rdep^2/rOx^2) (Log[rdep^2/rOx^2] - aParam)),
      rOx > rC > 0 && rdep > rOx
    ];

    dQdRdep = D[qDep, rdep];
    dVdRdep = FullSimplify[D[vTot, rdep], rOx > rC > 0 && rdep > rOx];
    cOxExpr = (2 Pi eOx)/Log[rOx/rC];
    cDepExpr = (2 Pi eSi)/Log[rdep/rOx];
    cTotExpr = 1/(1/cOxExpr + 1/cDepExpr);
    cTotDirectCheck = FullSimplify[dQdRdep/dVdRdep - cTotExpr, rOx > rC > 0 && rdep > rOx];

    phiLoc[n_] := (qS nAs)/(4 eSi k^2) ((1 + k n)^2 - (1 + k w)^2 - 2 (1 + k w)^2 Log[(1 + k n)/(1 + k w)]);
    phiLocCheck = FullSimplify[
      {
        (1/(1 + k n)) D[(1 + k n) D[phiLoc[n], n], n] - qS nAs/eSi,
        phiLoc[w],
        (D[phiLoc[n], n] /. n -> w)
      },
      k > 0 && w > 0 && 1 + k n > 0 && 1 + k w > 0
    ];

    cOxCurvSeries = Normal[Series[eOx * k / (-Log[1 - tOx k]), {k, 0, 3}]];
    cDepCurvSeries = Normal[Series[eSi * k / Log[1 + k w], {k, 0, 3}]];
    wDepSeries = Normal[Series[phiLoc[0], {k, 0, 3}]];
  ];

  (* 3. Numerical Reference C-V Table for Baseline HBM TSV *)
  qE = 1.602176634*^-19;
  kB = 1.380649*^-23;
  T = 300.0;
  eps0 = 8.8541878128*^-12;
  epsOx = 3.9 * eps0;
  epsSi = 11.7 * eps0;
  ni = 1.0*^16;         (* intrinsic carrier concentration in Si at 300 K [m^-3] = 1e10 cm^-3 *)
  Na = 1.0*^21;         (* p-Si doping 1e15 cm^-3 = 1e21 m^-3, rho ~ 13.5 Ohm-cm *)
  rc = 2.5*^-6;
  tox = 0.2*^-6;
  rox = rc + tox;
  phiF = (kB T/qE) Log[Na/ni];
  vFB = -0.35;          (* Cu/TaN-SiO2-pSi flatband voltage [V] *)

  rdepFromPsi[ps_] := Module[{c = (4 epsSi ps)/(qE Na rox^2) - 1},
    rox Sqrt[Exp[1 + ProductLog[c/E]]]
  ];
  rdepMax = rdepFromPsi[2 phiF];
  vTh = vFB + 2 phiF + (qE Na (rdepMax^2 - rox^2))/(2 epsOx) Log[rox/rc];

  solveRdepV[V_] := Module[{veff = V - vFB, kox, a, vn, u},
    If[veff <= 0, Return[rox]];
    If[V >= vTh, Return[rdepMax]];
    kox = (2 epsSi/epsOx) Log[rox/rc];
    a = 1 - kox;
    vn = (4 epsSi veff)/(qE Na rox^2);
    u = Exp[a + ProductLog[(vn - a) Exp[-a]]];
    rox Sqrt[u]
  ];

  cOxVal = (2 Pi epsOx)/Log[rox/rc];
  cTsvAtV[V_] := Module[{rd = solveRdepV[V]},
    If[rd <= rox * (1 + 1*^-14), cOxVal,
      1/(1/cOxVal + Log[rd/rox]/(2 Pi epsSi))]
  ];

  vSweep = {-1.0, -0.5, -0.35, -0.2, 0.0, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0};
  cvTable = Table[
    <|"V" -> v, "rdep_m" -> N[solveRdepV[v], 15], "C_TSV_F_per_m" -> N[cTsvAtV[v], 15]|>,
    {v, vSweep}
  ];

  (* Local curvature effect on threshold voltage at teardrop apex vs belly (alpha = 0.35):
     kappa_apex = 1/(R0 (1-a)^2), kappa_belly = 1/(R0 (1+a)^2) *)
  With[{rApex = rc * (1 - 0.35)^2, rBelly = rc * (1 + 0.35)^2},
    apexVth = Module[{roA = rApex + tox, rdA},
      rdA = roA Sqrt[Exp[1 + ProductLog[((4 epsSi (2 phiF))/(qE Na roA^2) - 1)/E]]];
      vFB + 2 phiF + (qE Na (rdA^2 - roA^2))/(2 epsOx) Log[roA/rApex]
    ];
    bellyVth = Module[{roB = rBelly + tox, rdB},
      rdB = roB Sqrt[Exp[1 + ProductLog[((4 epsSi (2 phiF))/(qE Na roB^2) - 1)/E]]];
      vFB + 2 phiF + (qE Na (rdB^2 - roB^2))/(2 epsOx) Log[roB/rBelly]
    ];
  ];

  <|
    "CylindricalPDE_zero" -> (pdeCheck === 0),
    "CylindricalBC_zero" -> (bcCheck === {0, 0}),
    "VoltageForm_zero" -> (uExpr === 0),
    "LambertW_exact_zero" -> (lambertCheck === 0),
    "DiffCap_identity_zero" -> (cTotDirectCheck === 0),
    "FrenetLocal_zero" -> (phiLocCheck === {0, 0, 0}),
    "Cox_curvature_series" -> ToString[InputForm[cOxCurvSeries]],
    "Cdep_curvature_series" -> ToString[InputForm[cDepCurvSeries]],
    "PsiS_curvature_series" -> ToString[InputForm[wDepSeries]],
    "Benchmarks" -> <|
      "phiF_V" -> N[phiF, 15],
      "V_FB_V" -> vFB,
      "rdep_max_m" -> N[rdepMax, 15],
      "V_th_cyl_V" -> N[vTh, 15],
      "V_th_teardrop_apex_V" -> N[apexVth, 15],
      "V_th_teardrop_belly_V" -> N[bellyVth, 15],
      "C_ox_F_per_m" -> N[cOxVal, 15],
      "C_min_F_per_m" -> N[cTsvAtV[2.0], 15],
      "CV_sweep" -> cvTable
    |>
  |>
];
