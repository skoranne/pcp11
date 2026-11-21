(* ::Package:: *)
(* File: 06_mems_esposito_ghoussoub_guo.wl *)
(* Formal symbolic & high-precision numerical verification of Esposito, Ghoussoub & Guo (2010), *)
(* "Mathematical Analysis of Partial Differential Equations Modeling Electrostatic MEMS" *)
(* (AMS Courant Lecture Notes, Vol. 20): Chapters 1-12 & Appendix A. *)

VerifyMEMSBook[] := Module[
  {muB1N2, muB1N3, lamUpperN2, lamUpperN3,
   uSing2nd, singRes2nd, lamSing2nd, hardyConst, critPoly2nd, nCrit2ndAlpha0, nCrit2ndAlpha1, nCrit2ndAlpha2,
   uSing4th, singRes4th, rellichConst, critPoly4th, nCrit4thAlpha0,
   boggio2D, boggio3D, boggio2DCheck, boggio3DCheck, boggio2DPos, boggio3DPos,
   quenchScalingCheck, selfSimC0, outerODECheck,
   shoot2ndN2, shoot2ndN3, shoot4thNavierN2, shoot4thClampedN2,
   branch2ndN2, branch2ndN3, fold2ndN2, fold2ndN3,
   fold4thNavierN2, fold4thClampedN2, branch4thNavierN2, branch4thClampedN2,
   hardyWeightedConst, rellichWeightedConst, quenchBench},

  (* === Chapter 2: Principal Dirichlet eigenvalue & pull-in voltage upper bound === *)
  muB1N2 = N[BesselJZero[0, 1]^2, 16];
  muB1N3 = N[Pi^2, 16];
  lamUpperN2 = (4/27) muB1N2;
  lamUpperN3 = (4/27) muB1N3;

  (* === Chapters 3-7 & Appendix A: Singular solution, Hardy inequality, & critical dimension N_crit(alpha) === *)
  (* 2nd-order singular solution on B_1: U_s(r) = 1 - r^((2+alpha)/3)
     Verify -Delta U_s = lam_s * r^alpha / (1 - U_s)^2 where lam_s = (2+alpha)(3N - 4 + alpha)/9 *)
  uSing2nd[r_] := 1 - r^((2 + alpha)/3);
  singRes2nd = FullSimplify[
    -(D[uSing2nd[r], {r, 2}] + ((NN - 1)/r) D[uSing2nd[r], r]) * (1 - uSing2nd[r])^2 / r^alpha,
    r > 0
  ];
  lamSing2nd = FullSimplify[((2 + alpha) (3 NN - 4 + alpha))/9];

  (* Linearized stability of U_s via Hardy inequality: 2 lam_s <= (N-2)^2 / 4 *)
  hardyConst = (NN - 2)^2 / 4;
  critPoly2nd = FullSimplify[hardyConst - 2 lamSing2nd];
  nCrit2ndAlpha0 = Max[NN /. Solve[(critPoly2nd /. alpha -> 0) == 0, NN]];
  nCrit2ndAlpha1 = Max[NN /. Solve[(critPoly2nd /. alpha -> 1) == 0, NN]];
  nCrit2ndAlpha2 = Max[NN /. Solve[(critPoly2nd /. alpha -> 2) == 0, NN]];

  (* === Part 3 (Ch. 11-12) & Appendix A: 4th-order singular solution, Rellich inequality, & critical dimension === *)
  uSing4th[r_] := 1 - r^((4 + alpha)/3);
  Module[{lapU, bilapU},
    lapU = D[uSing4th[r], {r, 2}] + ((NN - 1)/r) D[uSing4th[r], r];
    bilapU = D[lapU, {r, 2}] + ((NN - 1)/r) D[lapU, r];
    singRes4th = FullSimplify[bilapU * (1 - uSing4th[r])^2 / r^alpha, r > 0];
  ];
  rellichConst = NN^2 (NN - 4)^2 / 16;
  critPoly4th = FullSimplify[rellichConst - 2 singRes4th];
  nCrit4thAlpha0 = Max[Select[NN /. NSolve[(critPoly4th /. alpha -> 0) == 0, NN, Reals], # > 4 &]];

  (* Weighted Hardy & Rellich constants (Appendix A) *)
  hardyWeightedConst = ((NN - 2 + beta)/2)^2;
  rellichWeightedConst = ((NN - 4 + beta) (NN + beta)/4)^2;

  (* === Chapter 11: Boggio's Biharmonic Clamped Green's Function on Unit Ball B_1 === *)
  boggio2D[r_] := (1/(16 Pi)) (1 - r^2 + 2 r^2 Log[r]);
  Module[{lapB, bilapB, fluxB, val0},
    lapB = FullSimplify[(1/r) D[r D[boggio2D[r], r], r], r > 0];
    bilapB = FullSimplify[(1/r) D[r D[lapB, r], r], r > 0];
    fluxB = Simplify[2 Pi r * D[-lapB, r] /. r -> eps, eps > 0];
    val0 = Limit[boggio2D[r], r -> 0, Direction -> "FromAbove"];
    boggio2DCheck = {bilapB, boggio2D[1], (D[boggio2D[r], r] /. r -> 1), fluxB};
    boggio2DPos = (val0 === 1/(16 Pi)) && (MinValue[{boggio2D[r], 10^-6 <= r <= 1}, r] >= 0);
  ];

  boggio3D[r_] := (1/(16 Pi)) (1 - r)^2;
  Module[{lapB3, bilapB3, fluxB3},
    lapB3 = FullSimplify[(1/r^2) D[r^2 D[boggio3D[r], r], r], r > 0];
    bilapB3 = FullSimplify[(1/r^2) D[r^2 D[lapB3, r], r], r > 0];
    fluxB3 = Simplify[4 Pi r^2 * D[-lapB3, r] /. r -> eps, eps > 0];
    boggio3DCheck = {bilapB3, boggio3D[1], (D[boggio3D[r], r] /. r -> 1), fluxB3};
    boggio3DPos = (boggio3D[0] === 1/(16 Pi)) && (MinValue[{boggio3D[r], 0 <= r <= 1}, r] == 0);
  ];

  (* === Chapters 8-10: Parabolic Quenching Asymptotics === *)
  Module[{solQ},
    solQ = DSolve[{u'[t] == lam f0 / (1 - u[t])^2, u[0] == 0}, u[t], t];
    quenchScalingCheck = FullSimplify[(1 - (u[t] /. solQ[[1]]))^3 - (1 - 3 lam f0 t), 1 - 3 lam f0 t > 0];
  ];
  (* Center-manifold Hermite projection & outer matched self-similar ODE W(z) = (3 lam (1 + c0 z^2))^(1/3) *)
  selfSimC0 = Module[{wt = r Exp[-r^2/4], H2 = r^2 - 4, num, den},
    num = Integrate[H2^3 * wt, {r, 0, Infinity}];
    den = Integrate[H2^2 * wt, {r, 0, Infinity}];
    <|"num_over_den" -> Simplify[num/den], "a_linear_coeff" -> Simplify[den/num], "c0_cube_coeff" -> Simplify[3 den/num]|>
  ];
  outerODECheck = Module[{W},
    W[z_] := (3 lam f0 (1 + c0 z^2))^(1/3);
    FullSimplify[-(1/2) z W'[z] + (1/3) W[z] - lam f0 / W[z]^2, lam > 0 && f0 > 0 && c0 > 0 && z > 0]
  ];

  (* === Exact Radial Shooting Benchmarks on Unit Ball B_1 (2nd Order, N=2, 3) === *)
  shoot2ndN2[v0_?NumericQ] := Module[{v0r = SetPrecision[v0, 25], eps0 = 10^-9, rhoStar = 0},
    NDSolveValue[
      {v''[rho] + (1/rho) v'[rho] == 1/v[rho]^2,
       v[eps0] == v0r + eps0^2/(4 v0r^2), v'[eps0] == eps0/(2 v0r^2),
       WhenEvent[v[rho] == 1, rhoStar = rho; "StopIntegration"]},
      v, {rho, eps0, 50}, WorkingPrecision -> 20, AccuracyGoal -> 12, PrecisionGoal -> 12
    ];
    {N[rhoStar^2, 15], N[1 - v0, 15]}
  ];

  shoot2ndN3[v0_?NumericQ] := Module[{v0r = SetPrecision[v0, 25], eps0 = 10^-9, rhoStar = 0},
    NDSolveValue[
      {v''[rho] + (2/rho) v'[rho] == 1/v[rho]^2,
       v[eps0] == v0r + eps0^2/(6 v0r^2), v'[eps0] == eps0/(3 v0r^2),
       WhenEvent[v[rho] == 1, rhoStar = rho; "StopIntegration"]},
      v, {rho, eps0, 50}, WorkingPrecision -> 20, AccuracyGoal -> 12, PrecisionGoal -> 12
    ];
    {N[rhoStar^2, 15], N[1 - v0, 15]}
  ];

  branch2ndN2 = Table[shoot2ndN2[v0], {v0, Join[Range[95/100, 10/100, -5/100], {5/100, 2/100, 1/100}]}];
  branch2ndN3 = Table[shoot2ndN3[v0], {v0, Join[Range[95/100, 10/100, -5/100], {5/100, 2/100, 1/100}]}];

  fold2ndN2 = Module[{f2, res},
    f2[v0_?NumericQ] := shoot2ndN2[v0][[1]];
    res = FindMaximum[{f2[v0], 0.45 < v0 < 0.75}, {v0, 0.56}, AccuracyGoal -> 9, PrecisionGoal -> 9];
    <|"lambda_star" -> N[res[[1]], 12], "u_star_inf" -> N[1 - (v0 /. res[[2]]), 12]|>
  ];
  fold2ndN3 = Module[{f3, res},
    f3[v0_?NumericQ] := shoot2ndN3[v0][[1]];
    res = FindMaximum[{f3[v0], 0.35 < v0 < 0.70}, {v0, 0.50}, AccuracyGoal -> 9, PrecisionGoal -> 9];
    <|"lambda_star" -> N[res[[1]], 12], "u_star_inf" -> N[1 - (v0 /. res[[2]]), 12]|>
  ];

  (* === Fourth-Order Navier & Clamped Radial Benchmarks on Unit Disk B_1 (N = 2) === *)
  (* Fast bracketed secant/Brent in p0 (~10 NDSolve calls per v0 to 1e-10 accuracy) *)
  shoot4thNavierN2[v0_?NumericQ] := Module[{evalNav, lo = 0.01, hi = 1.5, fl, fh, mid, fm, rFinal},
    evalNav[p0_?NumericQ] := Module[{eps0 = 10^-6, rS = 0.0, vEnd = 5.0},
      Quiet@NDSolveValue[
        {v''[rho] + (1/rho) v'[rho] == p[rho],
         p''[rho] + (1/rho) p'[rho] == -1/v[rho]^2,
         v[eps0] == v0 + p0 eps0^2/4, v'[eps0] == p0 eps0/2,
         p[eps0] == p0 - eps0^2/(4 v0^2), p'[eps0] == -eps0/(2 v0^2),
         WhenEvent[p[rho] == 0, rS = rho; vEnd = v[rho]; "StopIntegration"],
         WhenEvent[v[rho] >= 5.0, rS = rho; vEnd = 5.0; "StopIntegration"]},
        {v, p}, {rho, eps0, 25}, AccuracyGoal -> 10, PrecisionGoal -> 10
      ];
      {vEnd - 1.0, rS}
    ];
    fl = evalNav[lo][[1]];
    fh = evalNav[hi][[1]];
    While[fh < 0 && hi < 50.0, hi = 1.8 hi; fh = evalNav[hi][[1]]];
    Do[
      mid = If[Abs[fh - fl] > 10^-14 && Mod[k, 3] != 0 && fh < 3.0,
        Clip[(lo fh - hi fl)/(fh - fl), {lo + 0.15 (hi - lo), hi - 0.15 (hi - lo)}],
        (lo + hi)/2
      ];
      {fm, rFinal} = evalNav[mid];
      If[Abs[fm] < 10^-10, Break[]];
      If[fm < 0, lo = mid; fl = fm, hi = mid; fh = fm],
      {k, 1, 35}
    ];
    {N[rFinal^4, 12], N[1 - v0, 12]}
  ];

  shoot4thClampedN2[v0_?NumericQ] := Module[{evalClamp, lo = 0.01, hi = 1.5, fl, fh, mid, fm, rFinal},
    evalClamp[p0_?NumericQ] := Module[{eps0 = 10^-6, rS = 0.0, vEnd = 5.0},
      Quiet@NDSolveValue[
        {v''[rho] + (1/rho) v'[rho] == p[rho],
         p''[rho] + (1/rho) p'[rho] == -1/v[rho]^2,
         v[eps0] == v0 + p0 eps0^2/4, v'[eps0] == p0 eps0/2,
         p[eps0] == p0 - eps0^2/(4 v0^2), p'[eps0] == -eps0/(2 v0^2),
         WhenEvent[v'[rho] == 0 && rho > 10^-3, rS = rho; vEnd = v[rho]; "StopIntegration"],
         WhenEvent[v[rho] >= 5.0, rS = rho; vEnd = 5.0; "StopIntegration"]},
        {v, p}, {rho, eps0, 25}, AccuracyGoal -> 10, PrecisionGoal -> 10
      ];
      {vEnd - 1.0, rS}
    ];
    fl = evalClamp[lo][[1]];
    fh = evalClamp[hi][[1]];
    While[fh < 0 && hi < 50.0, hi = 1.8 hi; fh = evalClamp[hi][[1]]];
    Do[
      mid = If[Abs[fh - fl] > 10^-14 && Mod[k, 3] != 0 && fh < 3.0,
        Clip[(lo fh - hi fl)/(fh - fl), {lo + 0.15 (hi - lo), hi - 0.15 (hi - lo)}],
        (lo + hi)/2
      ];
      {fm, rFinal} = evalClamp[mid];
      If[Abs[fm] < 10^-10, Break[]];
      If[fm < 0, lo = mid; fl = fm, hi = mid; fh = fm],
      {k, 1, 35}
    ];
    {N[rFinal^4, 12], N[1 - v0, 12]}
  ];

  branch4thNavierN2 = Table[shoot4thNavierN2[v0], {v0, Range[0.95, 0.25, -0.05]}];
  branch4thClampedN2 = Table[shoot4thClampedN2[v0], {v0, Range[0.95, 0.25, -0.05]}];

  (* Fast golden-section search on v0 for exact 4th-order folds *)
  Module[{goldenMax},
    goldenMax[fn_, a0_, b0_, iters_] := Module[{a = a0, b = b0, gr = (Sqrt[5] - 1)/2, c, d, fc, fd},
      c = b - gr (b - a); d = a + gr (b - a);
      fc = fn[c]; fd = fn[d];
      Do[
        If[fc < fd,
          a = c; c = d; fc = fd; d = a + gr (b - a); fd = fn[d],
          b = d; d = c; fd = fc; c = b - gr (b - a); fc = fn[c]
        ],
        {iters}
      ];
      With[{vBest = (a + b)/2}, <|"lambda_star" -> N[fn[vBest], 10], "u_star_inf" -> N[1 - vBest, 10]|>]
    ];
    fold4thNavierN2 = goldenMax[shoot4thNavierN2[#][[1]] &, 0.50, 0.75, 22];
    fold4thClampedN2 = goldenMax[shoot4thClampedN2[#][[1]] &, 0.50, 0.75, 22];
  ];

  (* === Chapter 9: Quenching Time Bounds on B_1 (N=2) for lambda > lambda* === *)
  (* Lower bound (ODE comparison): T*(lambda) >= 1/(3 lambda)
     Upper bound (first-eigenfunction Jensen projection, EGG Theorem 9.1):
     Let y(t) = Int u phi1 / Int phi1 -> y' >= -mu1 y + lam / (1 - y)^2 -> T* <= Int_0^1 dy / (lam/(1-y)^2 - mu1 y) *)
  quenchBench = Table[
    Module[{tLower, tUpper},
      tLower = 1/(3 lamVal);
      tUpper = NIntegrate[1/(lamVal/(1 - y)^2 - muB1N2 y), {y, 0, 1}, WorkingPrecision -> 15];
      <|"lambda" -> lamVal, "T_lower" -> N[tLower, 12], "T_upper" -> N[tUpper, 12]|>
    ],
    {lamVal, {1.0, 1.5, 2.5, 5.0, 10.0}}
  ];

  <|
    "Sing2nd_exact_zero" -> (FullSimplify[singRes2nd - lamSing2nd] === 0),
    "N_crit_2nd_alpha0" -> N[nCrit2ndAlpha0, 12],
    "N_crit_2nd_alpha1" -> N[nCrit2ndAlpha1, 12],
    "N_crit_2nd_alpha2" -> N[nCrit2ndAlpha2, 12],
    "N_crit_4th_alpha0" -> N[nCrit4thAlpha0, 12],
    "Sing4th_lam_formula" -> ToString[InputForm[singRes4th]],
    "Boggio2D_check" -> boggio2DCheck,
    "Boggio2D_positive" -> boggio2DPos,
    "Boggio3D_check" -> boggio3DCheck,
    "Boggio3D_positive" -> boggio3DPos,
    "QuenchODE_zero" -> (quenchScalingCheck === 0),
    "OuterSelfSimilarODE_zero" -> (outerODECheck === 0),
    "SelfSimilar_c0" -> selfSimC0,
    "Hardy_const" -> ToString[InputForm[hardyConst]],
    "Rellich_const" -> ToString[InputForm[rellichConst]],
    "Hardy_weighted" -> ToString[InputForm[hardyWeightedConst]],
    "Rellich_weighted" -> ToString[InputForm[rellichWeightedConst]],
    "Benchmarks" -> <|
      "mu_B1_N2" -> muB1N2,
      "mu_B1_N3" -> muB1N3,
      "lambda_upper_N2" -> lamUpperN2,
      "lambda_upper_N3" -> lamUpperN3,
      "fold_2nd_B1_N2" -> fold2ndN2,
      "fold_2nd_B1_N3" -> fold2ndN3,
      "fold_4th_Navier_B1_N2" -> fold4thNavierN2,
      "fold_4th_Clamped_B1_N2" -> fold4thClampedN2,
      "branch_2nd_B1_N2" -> N[branch2ndN2, 12],
      "branch_2nd_B1_N3" -> N[branch2ndN3, 12],
      "branch_4th_Navier_B1_N2" -> N[branch4thNavierN2, 12],
      "branch_4th_Clamped_B1_N2" -> N[branch4thClampedN2, 12],
      "quench_bounds_B1_N2" -> quenchBench
    |>
  |>
];
