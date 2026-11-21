(* ::Package:: *)
(* File: 02_curvilinear_geometry.wl *)
(* Symbolic verification of smooth 2D teardrop contours, curvature, offset SiO2 liner regularity, *)
(* Gauss-Bonnet turning number, and 3D tapered/scalloped teardrop TSV surface metric & curvatures. *)

VerifyCurvilinearGeometry[] := Module[
  {rTD, drExpr, d2rExpr, speedExpr, kappaExpr, kappa0, kappaPi, areaTD, turningCheck,
   convexBound, alphaConvexMax, offsetSpeedSq,
   R3D, Ru, Rz, g11, g12, g22, detG,
   areaNum, perimNum, kMaxNum, kMinNum, offsetMinNum, vol3DNum},

  (* Smooth C^infinity Teardrop contour (unrotated, tip along +x gas-escape/stress-relief axis):
     x(t) = R0 * Cos[t]
     y(t) = R0 * cy * Sin[t] * (1 - alpha * Cos[t])
     For 0 <= alpha < 1/2, the contour is strictly convex (kappa(t) > 0 for all t).
     For 1/2 <= alpha < 1, it forms a pinched/necked teardrop with sharp tip at t = 0. *)
  rTD[t_] := {R0 Cos[t], R0 cy Sin[t] (1 - alpha Cos[t])};
  drExpr = D[rTD[t], t];
  d2rExpr = D[rTD[t], {t, 2}];
  speedExpr = Sqrt[Simplify[drExpr.drExpr]];
  kappaExpr = Simplify[(drExpr[[1]] d2rExpr[[2]] - drExpr[[2]] d2rExpr[[1]]) / (drExpr.drExpr)^(3/2)];

  kappa0 = FullSimplify[kappaExpr /. t -> 0, R0 > 0 && cy > 0 && 0 < alpha < 1];
  kappaPi = FullSimplify[kappaExpr /. t -> Pi, R0 > 0 && cy > 0 && 0 < alpha < 1];

  (* Exact Green's theorem area: (1/2) Int (x y' - y x') dt = Pi R0^2 cy (independent of alpha!) *)
  areaTD = FullSimplify[(1/2) Integrate[rTD[t][[1]] drExpr[[2]] - rTD[t][[2]] drExpr[[1]], {t, 0, 2 Pi}]];

  (* Exact convexity threshold: numerator factor = 1 - alpha (3 Cos[t] - 2 Cos[t]^3) > 0 for all t iff alpha < 1/Sqrt[2] *)
  convexBound = FullSimplify[(drExpr[[1]] d2rExpr[[2]] - drExpr[[2]] d2rExpr[[1]])/(R0^2 cy)];
  alphaConvexMax = First[alpha /. Solve[MinValue[{1 - alpha (3 u - 2 u^3), -1 <= u <= 1}, u] == 0 && alpha > 0, alpha]];

  (* Offset curve r_ox(t) = r_c(t) + t_ox * n(t):
     Verify |r_ox'(t)| = |r_c'(t)| * (1 + t_ox * kappa(t)) *)
  offsetSpeedSq = FullSimplify[
    Module[{nvec, rox, drox},
      nvec = {y'[t], -x'[t]}/Sqrt[x'[t]^2 + y'[t]^2];
      rox = {x[t], y[t]} + tox nvec;
      drox = D[rox, t];
      Simplify[(drox.drox) - (x'[t]^2 + y'[t]^2) * (1 + tox * (x'[t] y''[t] - y'[t] x''[t])/(x'[t]^2 + y'[t]^2)^(3/2))^2]
    ],
    x'[t]^2 + y'[t]^2 > 0
  ];

  (* 3D Tapered + Scalloped Teardrop TSV surface: R(t, z) = (s(z) x(t), s(z) y(t), z) *)
  R3D[t_, z_] := {s[z] Cos[t], s[z] cy Sin[t] (1 - alpha Cos[t]), z};
  Ru = D[R3D[t, z], t];
  Rz = D[R3D[t, z], z];
  g11 = Simplify[Ru.Ru];
  g12 = Simplify[Ru.Rz];
  g22 = Simplify[Rz.Rz];
  detG = Simplify[g11 g22 - g12^2];

  (* Numerical benchmarks (exact rationals for 20-digit precision) *)
  With[{R0v = 5/2 * 10^-6, cyv = 1, av = 7/20, toxv = 2/10 * 10^-6, hv = 50 * 10^-6, tap = 1/20},
    areaNum = N[Pi R0v^2 cyv, 16];
    perimNum = NIntegrate[Evaluate[speedExpr /. {R0 -> R0v, cy -> cyv, alpha -> av}], {t, 0, 2 Pi}, WorkingPrecision -> 20];
    turningCheck = NIntegrate[Evaluate[(kappaExpr speedExpr) /. {R0 -> R0v, cy -> cyv, alpha -> av}], {t, 0, 2 Pi}, WorkingPrecision -> 20]/(2 Pi);
    kMaxNum = N[kappa0 /. {R0 -> R0v, cy -> cyv, alpha -> av}, 16];
    kMinNum = N[kappaPi /. {R0 -> R0v, cy -> cyv, alpha -> av}, 16];
    offsetMinNum = N[1 + toxv * kMinNum, 16];
    vol3DNum = N[Integrate[Pi (R0v (1 - tap z/hv))^2 cyv, {z, 0, hv}], 16];
  ];

  <|
    "Area_exact" -> areaTD,
    "Kappa_apex_t0" -> kappa0,
    "Kappa_belly_tPi" -> kappaPi,
    "Kappa_numerator_factor" -> TrigReduce[convexBound],
    "Alpha_convex_max" -> alphaConvexMax,
    "OffsetSpeedIdentity_zero" -> (offsetSpeedSq === 0),
    "DetG_3D_positive_at_alpha0" -> FullSimplify[detG /. {alpha -> 0, cy -> 1}],
    "GaussBonnet_TurningNumber" -> N[turningCheck, 15],
    "Benchmarks" -> <|
      "R0_m" -> 2.5*^-6,
      "alpha" -> 0.35,
      "alpha_convex_max" -> N[alphaConvexMax, 15],
      "Area_2D_m2" -> N[areaNum, 15],
      "Perimeter_2D_m" -> N[perimNum, 15],
      "Kappa_apex_1_per_m" -> N[kMaxNum, 15],
      "Kappa_belly_1_per_m" -> N[kMinNum, 15],
      "OffsetFactor_min" -> N[offsetMinNum, 15],
      "Volume_3D_tapered_m3" -> N[vol3DNum, 15]
    |>
  |>
];
