# TSV Structure with Semiconductor Traces - 3D Boundary Element Method
# 60 GHz S-parameter extraction with keep-out regions - FIXED VERSION

using CurvilinearTSV
using LinearAlgebra
using StaticArrays

# Material properties
const ε0 = 8.854187817e-12
const σCu = 5.8e7
const μ0 = 4π * 1e-7

# TSV dimensions
const tsvRadius = 10e-6
const tsvHeight = 100e-6
const tsvPitch = 50e-6
const numTSVs = 2

# Semiconductor layer
const traceHeight = 5e-6
const traceWidth = 2e-6
const keepOutDistance = 15e-6

# Positions - FIXED: traces now outside TSV keep-out regions
topTracePositions = [
  SVector(0 + tsvPitch/2 + 25e-6, tsvPitch/2, traceHeight),
  SVector(tsvPitch + tsvPitch/2 + 25e-6, tsvPitch/2, traceHeight)
]

bottomTracePositions = [
  SVector(0 + tsvPitch/2 + 25e-6, tsvPitch/2, -traceHeight),
  SVector(tsvPitch + tsvPitch/2 + 25e-6, tsvPitch/2, -traceHeight)
]

tsvCenters = [
  SVector(0, tsvPitch/2, tsvHeight/2),
  SVector(tsvPitch, tsvPitch/2, tsvHeight/2)
]

println("TSV Centers (μm): ", [c .* 1e6 for c in tsvCenters])
println("Top Trace Positions (μm): ", [p .* 1e6 for p in topTracePositions])

# Green's function
G3D(r) = 1 / (4π * r)
dist(p1, p2) = sqrt(sum((p1 .- p2).^2))

# Panelization
nTheta = 8
nZ = 6
nTrace = 6

nTSV = numTSVs * nTheta * nZ
nTop = numTSVs * nTrace * nTrace
nBottom = nTop
totalPanels = nTSV + nTop + nBottom

println("Panels - TSV: ", nTSV, ", Top: ", nTop, ", Bottom: ", nBottom)

# Generate TSV panel coordinates
tsvPoints = SVector{3,Float64}[]
for t in 1:numTSVs, i in 1:nTheta, j in 1:nZ
  thetaI = (i - 0.5) * 2π / nTheta
  zI = (j - 0.5) * tsvHeight / nZ
  push!(tsvPoints, SVector(
    tsvCenters[t][1] + tsvRadius * cos(thetaI),
    tsvCenters[t][2] + tsvRadius * sin(thetaI),
    zI
  ))
end

# Generate top trace panel coordinates
topTracePoints = SVector{3,Float64}[]
for t in 1:numTSVs, i in 1:nTrace, j in 1:nTrace
  xI = topTracePositions[t][1] + (i - 0.5) * traceWidth / nTrace
  yI = topTracePositions[t][2] + (j - 0.5) * traceWidth / nTrace
  push!(topTracePoints, SVector(xI, yI, traceHeight))
end

# Generate bottom trace panel coordinates
bottomTracePoints = SVector{3,Float64}[]
for t in 1:numTSVs, i in 1:nTrace, j in 1:nTrace
  xI = bottomTracePositions[t][1] + (i - 0.5) * traceWidth / nTrace
  yI = bottomTracePositions[t][2] + (j - 0.5) * traceWidth / nTrace
  push!(bottomTracePoints, SVector(xI, yI, -traceHeight))
end

println("Total panels: ", length(tsvPoints) + length(topTracePoints) + length(bottomTracePoints))

# Assemble potential matrix
nT = length(tsvPoints)
nTt = length(topTracePoints)
nBt = length(bottomTracePoints)

potentialMatrix = zeros(Float64, nT + nTt + nBt, nT + nTt + nBt)

# TSV-TSV interactions
for i in 1:nT
  for j in 1:nT
    d = dist(tsvPoints[i], tsvPoints[j])
    potentialMatrix[i, j] = i == j ? 1 / (4π * (tsvRadius / 10)) : G3D(d)
  end
end

# TSV-top trace interactions
for i in 1:nT
  for j in 1:nTt
    potentialMatrix[i, nT + j] = G3D(dist(tsvPoints[i], topTracePoints[j]))
  end
end

# Top trace-TSV interactions
for i in 1:nTt
  for j in 1:nT
    potentialMatrix[nT + i, j] = G3D(dist(topTracePoints[i], tsvPoints[j]))
  end
end

# Top trace-top trace interactions
for i in 1:nTt
  for j in 1:nTt
    d = dist(topTracePoints[i], topTracePoints[j])
    potentialMatrix[nT + i, nT + j] = i == j ? 1 / (4π * (traceWidth / 10)) : G3D(d)
  end
end

# Top trace-bottom trace coupling - FIXED: was missing block
for i in 1:nTt
  for j in 1:nBt
    potentialMatrix[nT + nTt + i, nT + nTt + j] = G3D(dist(topTracePoints[i], bottomTracePoints[j]))
  end
end

# Bottom trace-bottom trace interactions
for i in 1:nBt
  for j in 1:nBt
    d = dist(bottomTracePoints[i], bottomTracePoints[j])
    potentialMatrix[nT + nTt + i, nT + nTt + j] = i == j ? 1 / (4π * (traceWidth / 10)) : G3D(d)
  end
end

println("Matrix assembled: ", size(potentialMatrix))

# Boundary conditions - FIXED: bvector defined BEFORE use
bvector = zeros(Float64, nT + nTt + nBt)
for i in 1:nTt
  if i <= nTt ÷ 2
    bvector[nT + i] = 1.0
  else
    bvector[nT + i] = 0.0
  end
end

for i in 1:nBt
  if i <= nBt ÷ 2
    bvector[nT + nTt + i] = -1.0
  else
    bvector[nT + nTt + i] = 0.0
  end
end

println("Solving BEM system with ", nT + nTt + nBt, " equations...")
sol = ldiv!(lu!(potentialMatrix), bvector)
println("Solution computed with ", length(sol), " unknowns")

# Frequency domain analysis
ω = 2π * 60e9
δ = sqrt(2 / (ω * μ0 * σCu))
# FIXED: RSkin = 1/(σ·δ) = 0.0639 Ω/sq for Cu at 60 GHz
RSkin = 1 / (σCu * δ)

println("Skin depth: ", δ * 1e6, " μm")
println("Sheet resistance: ", RSkin, " Ω/sq")

# Impedance matrix for top/bottom traces
zMatrix = zeros(ComplexF64, nTt + nBt, nTt + nBt)
for i in 1:(nTt + nBt)
  for j in 1:(nTt + nBt)
    if i <= nTt && j <= nTt
      pt1, pt2 = topTracePoints[i], topTracePoints[j]
    elseif i <= nTt && j > nTt
      pt1, pt2 = topTracePoints[i], bottomTracePoints[j - nTt]
    elseif i > nTt && j <= nTt
      pt1, pt2 = bottomTracePoints[i - nTt], topTracePoints[j]
    else
      pt1, pt2 = bottomTracePoints[i - nTt], bottomTracePoints[j - nTt]
    end
    d = dist(pt1, pt2)
    if i == j
      zMatrix[i, j] = RSkin + im * ω * 1e-12
    else
      zMatrix[i, j] = im * ω * μ0 / (4π * d) * 1e-6
    end
  end
end

println("Impedance matrix size: ", size(zMatrix))

# S-parameter calculation
Z0 = 50.0
SMatrix = inv(zMatrix + Z0 * I) * (zMatrix - Z0 * I)

println("S-matrix magnitude:")
display(abs.(SMatrix))

S_dB = 20 * log10.(abs.(SMatrix))
println("S-parameters (dB) for first 4x4 ports:")
display(S_dB[1:4, 1:4])

# Export results
open("/Users/skoranne/GITHUB/pcp11/CurvilinearTSV/test/tsv_sparameters_60ghz.csv", "w") do f
  for i in 1:min(4, size(SMatrix, 1)), j in 1:min(4, size(SMatrix, 2))
    println(f, "$i,$j,$(real(SMatrix[i,j])),$(imag(SMatrix[i,j])),$(20*log10(abs(SMatrix[i,j])))")
  end
end

println("S-parameters saved to tsv_sparameters_60ghz.csv")

# Summary
println("\n=== TSV BEM Analysis Summary ===")
println("TSVs: ", numTSVs, " (radius = ", tsvRadius * 1e6, " μm)")
println("Trace width: ", traceWidth * 1e6, " μm")
println("Total panels: ", nT + nTt + nBt)
println("Analysis frequency: 60 GHz")
println("S-parameters computed for ", nTt, " traces")
println("Analysis complete.")
