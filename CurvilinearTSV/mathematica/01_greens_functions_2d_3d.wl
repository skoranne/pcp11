(* ::Package:: *)
(* File: 01_greens_functions_2d_3d.wl *)
(* Symbolic verification of 2D/3D Laplace, Biharmonic-Poisson, and Helmholtz Green's functions *)
(* plus exact analytical capacitance formulas for concentric, eccentric, twin, sphere, and spheroid benchmarks. *)

VerifyGreensFunctions[] := Module[
  {g2d, g3d, lap2d, lap3d, flux2d, flux3d, dlp2dInt, dlp2dBnd, dlp2dExt,
   dlp3dInt, dlp3dBnd, dlp3dExt, psi2d, psi3d, bih2dCheck, bih3dCheck,
   phiP2d, phiP3d, partCheck2d, partCheck3d,
   g2dH, g3dH, helm2dCheck, helm3dCheck,
   cEcc, cEccLimit, cProlate, cProlateLimit,
   eps0, cEccRef, cTwinRef, cSphereRef, cProlateRef},

  (* 1. 2D and 3D Laplace Green's functions *)
  g2d = -Log[Sqrt[x^2 + y^2]]/(2 Pi);
  g3d = 1/(4 Pi Sqrt[x^2 + y^2 + z^2]);
  lap2d = FullSimplify[Laplacian[g2d, {x, y}], x^2 + y^2 > 0];
  lap3d = FullSimplify[Laplacian[g3d, {x, y, z}], x^2 + y^2 + z^2 > 0];

  (* Distributional flux across epsilon-ball confirms Laplacian[G] = -Delta *)
  flux2d = Simplify[Integrate[D[-Log[r]/(2 Pi), r] * r, {theta, 0, 2 Pi}], r > 0];
  flux3d = Simplify[Integrate[D[1/(4 Pi r), r] * r^2 Sin[theta], {theta, 0, Pi}, {phi, 0, 2 Pi}], r > 0];

  (* Double-layer jump relations on unit circle: Int_{S^1} (dG/dn') ds' *)
  dlp2dInt = Simplify[Integrate[(rho Cos[t] - 1)/(2 Pi (1 - 2 rho Cos[t] + rho^2)), {t, 0, 2 Pi}, Assumptions -> 0 < rho < 1]];
  dlp2dBnd = Simplify[Integrate[-1/(4 Pi), {t, 0, 2 Pi}]]; (* since (Cos[t]-1)/(2-2 Cos[t]) == -1/2 identically *)
  dlp2dExt = Simplify[Integrate[(rho Cos[t] - 1)/(2 Pi (1 - 2 rho Cos[t] + rho^2)), {t, 0, 2 Pi}, Assumptions -> rho > 1]];

  (* Double-layer jump relations on unit sphere via u = Cos[th], du = -Sin[th] dth *)
  dlp3dInt = Simplify[Integrate[(rho u - 1)/(2 (1 - 2 rho u + rho^2)^(3/2)), {u, -1, 1}, Assumptions -> 0 < rho < 1]];
  dlp3dBnd = Simplify[Integrate[(u - 1)/(2 (2 - 2 u)^(3/2)), {u, -1, 1}]];
  dlp3dExt = Simplify[Integrate[(rho u - 1)/(2 (1 - 2 rho u + rho^2)^(3/2)), {u, -1, 1}, Assumptions -> rho > 1]];

  (* 2. Biharmonic fundamental solutions psi: Laplacian[psi] = G (exact volume-to-boundary Poisson reduction) *)
  psi2d = -r^2 (Log[r] - 1)/(8 Pi);
  psi3d = r/(8 Pi);
  bih2dCheck = FullSimplify[(1/r) D[r D[psi2d, r], r] - (-Log[r]/(2 Pi)), r > 0];
  bih3dCheck = FullSimplify[(1/r^2) D[r^2 D[psi3d, r], r] - (1/(4 Pi r)), r > 0];

  (* Particular solutions for constant charge Laplacian[phi_p] = q Na / epsSi *)
  phiP2d = (q Na)/(4 epsSi) (x^2 + y^2);
  phiP3d = (q Na)/(6 epsSi) (x^2 + y^2 + z^2);
  partCheck2d = FullSimplify[Laplacian[phiP2d, {x, y}] - q Na/epsSi];
  partCheck3d = FullSimplify[Laplacian[phiP3d, {x, y, z}] - q Na/epsSi];

  (* 3. 2D & 3D Helmholtz Green's functions *)
  g2dH = (I/4) HankelH1[0, k r];
  g3dH = Exp[I k r]/(4 Pi r);
  helm2dCheck = FullSimplify[(1/r) D[r D[g2dH, r], r] + k^2 g2dH, r > 0];
  helm3dCheck = FullSimplify[(1/r^2) D[r^2 D[g3dH, r], r] + k^2 g3dH, r > 0];

  (* 4. Analytical Capacitance Identities *)
  cEcc[eps_, R1_, R2_, d_] := (2 Pi eps)/ArcCosh[(R1^2 + R2^2 - d^2)/(2 R1 R2)];
  cEccLimit = FullSimplify[cEcc[eps, R1, R2, 0] - (2 Pi eps)/Log[R2/R1], R2 > R1 > 0];

  cProlate[eps_, a_, c_] := Module[{e = Sqrt[1 - a^2/c^2]}, (4 Pi eps c e)/ArcTanh[e]];
  cProlateLimit = FullSimplify[Limit[cProlate[eps, a, c], a -> c, Direction -> "FromBelow"], c > 0];

  (* Numerical benchmark constants (SI units, F/m or F) *)
  eps0 = 8.8541878128*^-12;
  cEccRef = N[cEcc[3.9 eps0, 2.5*^-6, 5.0*^-6, 1.0*^-6], 16];
  cTwinRef = N[(Pi * 11.7 eps0)/ArcCosh[20.0*^-6 / (2 * 2.5*^-6)], 16];
  cSphereRef = N[(4 Pi * 3.9 eps0 * 2.5*^-6 * 3.5*^-6)/(3.5*^-6 - 2.5*^-6), 16];
  cProlateRef = N[cProlate[eps0, 1.0, 2.5], 16];

  <|
    "Lap2D_zero" -> (lap2d === 0),
    "Lap3D_zero" -> (lap3d === 0),
    "Flux2D" -> flux2d,
    "Flux3D" -> flux3d,
    "DLP2D_Jumps" -> {dlp2dInt, dlp2dBnd, dlp2dExt},
    "DLP3D_Jumps" -> {dlp3dInt, dlp3dBnd, dlp3dExt},
    "Biharmonic2D_zero" -> (bih2dCheck === 0),
    "Biharmonic3D_zero" -> (bih3dCheck === 0),
    "Particular2D_zero" -> (partCheck2d === 0),
    "Particular3D_zero" -> (partCheck3d === 0),
    "Helmholtz2D_zero" -> (helm2dCheck === 0),
    "Helmholtz3D_zero" -> (helm3dCheck === 0),
    "EccentricLimit_zero" -> (cEccLimit === 0),
    "ProlateLimit_sphere" -> FullSimplify[cProlateLimit/(4 Pi eps c), c > 0],
    "Benchmarks" -> <|
      "C_ecc_2D_F_per_m" -> cEccRef,
      "C_twin_2D_F_per_m" -> cTwinRef,
      "C_sphere_3D_F" -> cSphereRef,
      "C_prolate_3D_eps0_a1_c2p5_F" -> cProlateRef
    |>
  |>
];
