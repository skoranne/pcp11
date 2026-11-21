(* ::Package:: *)
(* File: run_all_checks.wl *)
(* Master symbolic verification & reference benchmark generator for CurvilinearTSV *)

RunAllCurvilinearTSVChecks[dir_String] := Module[
  {r1, r2, r3, r4, r5, r6, allRes, outFile},
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

  allRes = <|
    "01_greens" -> r1,
    "02_geometry" -> r2,
    "03_depletion" -> r3,
    "04_skin_rlgc" -> r4,
    "05_sparameters" -> r5,
    "06_mems_egg2010" -> r6
  |>;

  (* Export purely numerical/string-safe reference JSON for Julia tests *)
  outFile = FileNameJoin[{dir, "reference_data.json"}];
  Export[outFile,
    <|
      "greens" -> r1["Benchmarks"],
      "geometry" -> r2["Benchmarks"],
      "depletion" -> r3["Benchmarks"],
      "skin_rlgc" -> r4["Benchmarks"],
      "sparameters" -> r5["Benchmarks"],
      "mems_egg2010" -> Join[
        r6["Benchmarks"],
        <|
          "N_crit_2nd_alpha0" -> r6["N_crit_2nd_alpha0"],
          "N_crit_2nd_alpha1" -> r6["N_crit_2nd_alpha1"],
          "N_crit_2nd_alpha2" -> r6["N_crit_2nd_alpha2"],
          "N_crit_4th_alpha0" -> r6["N_crit_4th_alpha0"]
        |>
      ]
    |>,
    "JSON"
  ];
  allRes
];
