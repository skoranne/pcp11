(* TSV Structure with Semiconductor Traces - 3D Boundary Element Method *)
(* 60 GHz S-parameter extraction with keep-out regions *)

ClearAll["Global`*"]

(* Material properties *)
ε0 = 8.854187817*^-12;
σCu = 5.8*10^7;
μ0 = 4*π*10^-7;

(* TSV dimensions *)
tsvRadius = 10*10^-6;
tsvHeight = 100*10^-6;
tsvPitch = 50*10^-6;
numTSVs = 2;

(* Semiconductor layer *)
traceHeight = 5*10^-6;
traceWidth = 2*10^-6;
keepOutDistance = 15*10^-6;

(* Positions *)
topTracePositions = Table[{tsvPitch/2 + i*tsvPitch + 20*10^-6, tsvPitch/2, traceHeight}, {i, 0, numTSVs - 1}];
bottomTracePositions = Table[{tsvPitch/2 + i*tsvPitch + 20*10^-6, tsvPitch/2, -traceHeight}, {i, 0, numTSVs - 1}];
tsvCenters = Table[{i*tsvPitch, tsvPitch/2, tsvHeight/2}, {i, 0, numTSVs - 1}];

Print["TSV Centers (μm): ", N[tsvCenters/10^-6]];
Print["Top Trace Positions (μm): ", N[topTracePositions/10^-6]];
Print["Bottom Trace Positions (μm): ", N[bottomTracePositions/10^-6]];

(* Green's function *)
G3D[r_] := 1/(4*π*r);
dist[p1_, p2_] := Sqrt[(p1[[1]] - p2[[1]])^2 + (p1[[2]] - p2[[2]])^2 + (p1[[3]] - p2[[3]])^2];

(* Panelization *)
nTheta = 8;
nZ = 6;
nTrace = 6;

nTSV = numTSVs * nTheta * nZ;
nTop = numTSVs * nTrace * nTrace;
nBottom = nTop;
totalPanels = nTSV + nTop + nBottom;

Print["Panels - TSV: ", nTSV, ", Top: ", nTop, ", Bottom: ", nBottom];

(* Generate TSV panel coordinates *)
tsvPoints = {};
Do[
  thetaI = (i - 0.5)*2*π/nTheta;
  zI = (j - 0.5)*tsvHeight/nZ;
  AppendTo[tsvPoints, {
    tsvCenters[[t]][[1]] + tsvRadius*Cos[thetaI],
    tsvCenters[[t]][[2]] + tsvRadius*Sin[thetaI], zI
  }],
  {t, 1, numTSVs}, {i, 1, nTheta}, {j, 1, nZ}
];

(* Generate top trace panel coordinates *)
topTracePoints = {};
Do[
  xI = topTracePositions[[t]][[1]] + (i - 0.5)*traceWidth/nTrace;
  yI = topTracePositions[[t]][[2]] + (j - 0.5)*traceWidth/nTrace;
  AppendTo[topTracePoints, {xI, yI, traceHeight}],
  {t, 1, numTSVs}, {i, 1, nTrace}, {j, 1, nTrace}
];

(* Generate bottom trace panel coordinates *)
bottomTracePoints = {};
Do[
  xI = bottomTracePositions[[t]][[1]] + (i - 0.5)*traceWidth/nTrace;
  yI = bottomTracePositions[[t]][[2]] + (j - 0.5)*traceWidth/nTrace;
  AppendTo[bottomTracePoints, {xI, yI, -traceHeight}],
  {t, 1, numTSVs}, {i, 1, nTrace}, {j, 1, nTrace}
];

Print["Total panels: ", Length[tsvPoints] + Length[topTracePoints] + Length[bottomTracePoints]];

(* Assemble potential matrix using nested loops *)
nT = Length[tsvPoints];
nTt = Length[topTracePoints];
nBt = Length[bottomTracePoints];

potentialMatrix = Table[0., {i, 1, nT + nTt + nBt}, {j, 1, nT + nTt + nBt}];

(* TSV-TSV interactions *)
Do[
  ri = tsvPoints[[i]];
  Do[
    rj = tsvPoints[[j]];
    d = dist[ri, rj];
    potentialMatrix[[i, j]] = If[i == j, 1/(4*π*(tsvRadius/10)), G3D[d]],
    {j, 1, nT}
  ],
  {i, 1, nT}
];

(* TSV-top trace interactions *)
Do[
  ri = tsvPoints[[i]];
  Do[
    rj = topTracePoints[[j]];
    potentialMatrix[[i, nT + j]] = G3D[dist[ri, rj]],
    {j, 1, nTt}
  ],
  {i, 1, nT}
];

(* Top trace-TSV interactions *)
Do[
  ri = topTracePoints[[i]];
  Do[
    rj = tsvPoints[[j]];
    potentialMatrix[[nT + i, j]] = G3D[dist[ri, rj]],
    {j, 1, nT}
  ],
  {i, 1, nTt}
];

(* Top trace-top trace interactions *)
Do[
  ri = topTracePoints[[i]];
  Do[
    rj = topTracePoints[[j]];
    d = dist[ri, rj];
    potentialMatrix[[nT + i, nT + j]] = If[i == j, 1/(4*π*(traceWidth/10)), G3D[d]],
    {j, 1, nTt}
  ],
  {i, 1, nTt}
];

(* Top trace-bottom trace coupling *)
Do[
  ri = topTracePoints[[i]];
  Do[
    rj = bottomTracePoints[[j]];
    potentialMatrix[[nT + nTt + i, nT + nTt + j]] = G3D[dist[ri, rj]],
    {j, 1, nBt}
  ],
  {i, 1, nTt}
];

(* Bottom trace-bottom trace interactions *)
Do[
  ri = bottomTracePoints[[i]];
  Do[
    rj = bottomTracePoints[[j]];
    d = dist[ri, rj];
    potentialMatrix[[nT + nTt + i, nT + nTt + j]] = If[i == j, 1/(4*π*(traceWidth/10)), G3D[d]],
    {j, 1, nBt}
  ],
  {i, 1, nBt}
];

Print["Matrix assembled: ", Dimensions[potentialMatrix]];

(* Boundary conditions *)
bvector = ConstantArray[0., nT + nTt + nBt];
Do[bvector[[nT + i]] = 1.0, {i, 1, nTt}];
Do[bvector[[nT + nTt + i]] = -1.0, {i, 1, nBt}];

Print["Solving BEM system with ", nT + nTt + nBt, " equations..."];
sigma = LinearSolve[potentialMatrix, bvector];
Print["Solution computed with ", Length[sigma], " unknowns"];

(* Compute capacitance matrix (first 8x8) *)
Print["Computing capacitance matrix (first 8x8)"];
capacitanceMatrix = Table[0., {i, 1, 8}, {j, 1, 8}];

Do[
  vj = ConstantArray[0., nT + nTt + nBt];
  vj[[j]] = 1.;
  sigmaJ = LinearSolve[potentialMatrix, vj];
  
  (* Charge on panel p due to panel j potential *)
  Do[
    capacitanceMatrix[[p, j]] += sigmaJ[[p]] * (2*π*tsvRadius*tsvHeight/(nTheta*nZ)),
    {p, 1, Min[4, nT]}
  ],
  {j, 1, 8}
];

Print["Capacitance matrix (8x8) in fF:"];
MatrixForm[N[capacitanceMatrix*10^15]];

(* Frequency domain analysis *)
ω = 2*π*60*10^9;
δ = Sqrt[2/(ω*μ0*σCu)];
RSkin = 1/(σCu*traceWidth*δ);

Print["Skin depth: ", N[δ*10^6], " μm"];
Print["Sheet resistance: ", N[RSkin], " Ω"];

(* Impedance matrix for top traces *)
zMatrix = Table[
  If[i == j, RSkin + I*ω*10^-12, I*ω*10^-12],
  {i, 1, nTt}, {j, 1, nTt}
];

Print["Impedance matrix (first 4x4) in milliohms:"];
MatrixForm[N[zMatrix[[1 ;; 4, 1 ;; 4]]*1000]];

(* S-parameter calculation *)
Z0 = 50.;
SMatrix = Inverse[zMatrix + Z0*IdentityMatrix[nTt]].(zMatrix - Z0*IdentityMatrix[nTt]);

Print["S-matrix magnitude:"];
MatrixForm[Abs[SMatrix]];

S_dB = 20*Log[10, Abs[SMatrix]];
Print["S-parameters (dB):"];
MatrixForm[N[S_dB, 3]];

(* Export results *)
Export["tsv_sparameters_60ghz.csv",
  Flatten[Table[
    {i, j, Re[SMatrix[[i, j]]], Im[SMatrix[[i, j]]], 20*Log[10, Abs[SMatrix[[i, j]]]]},
    {i, 1, Min[4, nTt]}, {j, 1, Min[4, nTt]}
  ], 1],
  "CSV"];

Print["S-parameters saved to tsv_sparameters_60ghz.csv"];

(* 3D visualization *)
Print["Generating 3D visualization"];

tsvGraphics = Table[
  Graphics3D[Cylinder[{{tsvCenters[[i + 1]][[1]], tsvCenters[[i + 1]][[2]], 0},
    {tsvCenters[[i + 1]][[1]], tsvCenters[[i + 1]][[2]], tsvHeight}}, tsvRadius],
    PlotStyle -> Directive[Gray, Opacity[0.3]], Boxed -> False],
  {i, 0, numTSVs - 1}];

topTraceGraphics = Table[
  Graphics3D[Cuboid[
    {topTracePositions[[i + 1]][[1]] - traceWidth/2, topTracePositions[[i + 1]][[2]] - traceWidth/2,
     traceHeight - 0.1*10^-6},
    {topTracePositions[[i + 1]][[1]] + traceWidth/2, topTracePositions[[i + 1]][[2]] + traceWidth/2,
     traceHeight + 0.1*10^-6}],
    PlotStyle -> Directive[Orange, Opacity[0.8]], Boxed -> False],
  {i, 0, numTSVs - 1}];

bottomTraceGraphics = Table[
  Graphics3D[Cuboid[
    {bottomTracePositions[[i + 1]][[1]] - traceWidth/2, bottomTracePositions[[i + 1]][[2]] - traceWidth/2,
     -traceHeight - 0.1*10^-6},
    {bottomTracePositions[[i + 1]][[1]] + traceWidth/2, bottomTracePositions[[i + 1]][[2]] + traceWidth/2,
     -traceHeight + 0.1*10^-6}],
    PlotStyle -> Directive[Orange, Opacity[0.8]], Boxed -> False],
  {i, 0, numTSVs - 1}];

keepOutGraphics = Table[
  Graphics3D[Cylinder[{{tsvCenters[[i + 1]][[1]], tsvCenters[[i + 1]][[2]], 0},
    {tsvCenters[[i + 1]][[1]], tsvCenters[[i + 1]][[2]], tsvHeight}}, keepOutDistance],
    PlotStyle -> Directive[Red, Opacity[0.2]], Boxed -> False],
  {i, 0, numTSVs - 1}];

combinedPlot = Show[tsvGraphics, topTraceGraphics, bottomTraceGraphics, keepOutGraphics,
  PlotRange -> All, AxesLabel -> {"x (μm)", "y (μm)", "z (μm)"},
  Boxed -> True, Lighting -> "Neutral", ImageSize -> Large];

Export["tsv_structure.png", combinedPlot, ImageSize -> 1200];
Print["3D structure saved to tsv_structure.png"];

Print[""];
Print["=== TSV BEM Analysis Summary ==="];
Print["TSVs: ", numTSVs, " (radius = ", N[tsvRadius*10^6], " μm)"];
Print["Trace width: ", N[traceWidth*10^6], " μm"];
Print["Total panels: ", nT + nTt + nBt];
Print["Analysis frequency: 60 GHz"];
Print["S-parameters computed for ", nTt, " traces"];
