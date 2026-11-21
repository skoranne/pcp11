(* ::Package:: *)
(* Automated Wolfram Language Test Suite (.wlt) for CurvilinearTSV & Esposito-Ghoussoub-Guo (2010) *)

dir = "/Users/skoranne/GITHUB/pcp11/CurvilinearTSV/mathematica";
Get[FileNameJoin[{dir, "01_greens_functions_2d_3d.wl"}]];
Get[FileNameJoin[{dir, "02_curvilinear_geometry.wl"}]];
Get[FileNameJoin[{dir, "03_depletion_theory.wl"}]];
Get[FileNameJoin[{dir, "04_skin_impedance_rlgc.wl"}]];
Get[FileNameJoin[{dir, "05_sparameters_hbm.wl"}]];
Get[FileNameJoin[{dir, "06_mems_esposito_ghoussoub_guo.wl"}]];

r1 = VerifyGreensFunctions[];
r2 = VerifyCurvilinearGeometry[];
r3 = VerifyDepletionTheory[];
r4 = VerifySkinAndRLGC[];
r5 = VerifySParametersHBM[];
r6 = VerifyMEMSBook[];

VerificationTest[r1["Lap2D_zero"] && r1["Lap3D_zero"], True, TestID -> "Greens_Laplace_2D_3D"]
VerificationTest[r1["DLP2D_Jumps"], {-1, -1/2, 0}, TestID -> "Greens_DLP2D_Jump"]
VerificationTest[r1["DLP3D_Jumps"], {-1, -1/2, 0}, TestID -> "Greens_DLP3D_Jump"]
VerificationTest[r1["Biharmonic2D_zero"] && r1["Biharmonic3D_zero"], True, TestID -> "Greens_Biharmonic_Poisson"]
VerificationTest[r1["Helmholtz2D_zero"] && r1["Helmholtz3D_zero"], True, TestID -> "Greens_Helmholtz_2D_3D"]
VerificationTest[r1["EccentricLimit_zero"], True, TestID -> "Cap_Eccentric_Limit"]
VerificationTest[r2["OffsetSpeedIdentity_zero"], True, TestID -> "Geom_OffsetSpeed_Identity"]
VerificationTest[r2["Alpha_convex_max"], 1/Sqrt[2], TestID -> "Geom_Teardrop_ConvexityBound"]
VerificationTest[Abs[r2["GaussBonnet_TurningNumber"] - 1.0] < 10^-13, True, TestID -> "Geom_GaussBonnet"]
VerificationTest[r3["CylindricalPDE_zero"] && r3["CylindricalBC_zero"], True, TestID -> "Dep_Cylindrical_PDE_BC"]
VerificationTest[r3["LambertW_exact_zero"], True, TestID -> "Dep_Exact_LambertW"]
VerificationTest[r3["DiffCap_identity_zero"], True, TestID -> "Dep_Differential_Capacitance"]
VerificationTest[r3["FrenetLocal_zero"], True, TestID -> "Dep_Frenet_Local_PDE"]
VerificationTest[r4["SolidDC_R_exact"] === 0 && r4["SolidDC_Lint_exact"] === 0, True, TestID -> "Skin_Solid_DC"]
VerificationTest[r4["CladDC_rel_err"] < 10^-14, True, TestID -> "Skin_Cladded_DC"]
VerificationTest[r5["Reciprocity_zero"] && r5["LosslessUnitarity_zero"], True, TestID -> "SParam_Reciprocity_Unitarity"]
VerificationTest[r6["Sing2nd_exact_zero"], True, TestID -> "EGG_2ndOrder_SingularSolution"]
VerificationTest[Abs[r6["N_crit_2nd_alpha0"] - (14 + 4 Sqrt[6])/3] < 10^-10, True, TestID -> "EGG_2ndOrder_Ncrit"]
VerificationTest[r6["N_crit_4th_alpha0"] > 8 && r6["N_crit_4th_alpha0"] < 9, True, TestID -> "EGG_4thOrder_Ncrit_8"]
VerificationTest[r6["Boggio2D_positive"] && r6["Boggio3D_positive"], True, TestID -> "EGG_Boggio_Positivity"]
VerificationTest[r6["QuenchODE_zero"] && r6["OuterSelfSimilarODE_zero"], True, TestID -> "EGG_Quenching_Asymptotics"]
VerificationTest[r6["Benchmarks"]["fold_2nd_B1_N2"]["lambda_star"] < r6["Benchmarks"]["lambda_upper_N2"], True, TestID -> "EGG_PullIn_UpperBound_N2"]
VerificationTest[r6["Benchmarks"]["fold_4th_Clamped_B1_N2"]["lambda_star"] > r6["Benchmarks"]["fold_4th_Navier_B1_N2"]["lambda_star"], True, TestID -> "EGG_Clamped_gt_Navier"]
