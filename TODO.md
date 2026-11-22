# TODO - 60GHz TSV S-Parameter Validation / Paper

**Status: NOT PUBLICATION READY**  
**Date:** 2026-11-21  
**Repos:** `pcp11` (this repo), `/Users/skoranne/GITHUB/VLSI_LAYOUT`

## Objective
Address review-identified problems for 60 GHz TSV/trace S-parameter validation and regenerate a comprehensive, honest LaTeX paper.

## Completed

### gds2mesh / layoutdb (VLSI_LAYOUT)
- [x] Fix micron/OpenCASCADE scaling: `unit_scale=1e6` default.
- [x] Strict Gmsh exit code; physical-tag element assertions in tests.
- [x] Hierarchy flattening: `ldb_iterate_boxes_flat` C API, Julia `get_boxes_flat`, `flatten=true`.
- [x] Palace physical-surface tagging in 3D extrusion: `absorbing_surfaces`, `pec_volume_names`.
- [x] Polygon-vertex C API scaffold: `ldb_iterate_polygons/flat`, `LdbPoint/LdbPolyElem`, `LdbCell.polys`.
- [x] GDSII BOUNDARY/PATH/BOX XY vertex storage in `gdsii.cpp`.
- [x] Julia polygon wrappers: `PolyElem`, `get_polygons/flat`, `write_geo_2d_polygons`.
- [x] Tests: **34/34 pass**.

**Evidence:** `cd /Users/skoranne/GITHUB/VLSI_LAYOUT/src/gds2mesh && julia --project -e 'push!(LOAD_PATH,"src"); include("src/GDS2Mesh.jl"); using .GDS2Mesh; include("test/runtests.jl")'`  
**Commits:** `03ec55a`, `4982ecf`, `18ca3bc`, `8b9e99e`

### Julia baselines / BEM mitigation (pcp11)
- [x] Corrected analytical 4-port RLC baseline: `CurvilinearTSV/test/tsv_bem_corrected.jl`, `tsv_sparameters_60ghz_corrected.csv`. Passive/reciprocal.
- [x] 2D quasi-static Laplace BEM+MTL extraction: `CurvilinearTSV/test/bem2d_capacitance.jl`. 60GHz `S11=-22.329 dB, S21=-0.501 dB, cross=-85.857 dB` (30 µm, 2-trace cross-section).
- [x] Simple2 MTL baseline: `CurvilinearTSV/test/tsv_mtl_simple2.jl`. 60GHz `S11=-4.30 dB, S21=-2.39 dB`, passive/reciprocal.

### Palace debugging / evidence (pcp11)
- [x] v2 debug: fixed PEC-on-interface, bottom-gap misclassification. Result invalid: `S11=-7.532e-5 dB, S21=-257.492 dB`, err norm `0.909`.
- [x] v3 edge-port redesign: blocked by Gmsh OCC/MFEM mesh topology (`Invalid mesh topology. Interior triangular face found...`).
- [x] simple2 model built with Si/SiO2/Cu/PEC/absorbing/lumped ports. 1GHz passive/reciprocal: `max|S-S^T|=1.3e-9, sigma_max=0.632`.
- [x] simple2 60GHz convergence sweep: O1 err 0.188; O2 25k err 0.0246; O2 42k err 0.0205; O2+AMR(1) err 0.0188; O2+AMR(2) err 0.0140. S21 stable ~`-30.71 dB`.
- [x] Port redesign attempts: z-restricted trace-to-ground port and stub extension; no valid outer-boundary port faces in corridor after OCC BooleanFragments. Reverted; documented.

**Best Palace evidence:** `CurvilinearTSV/test/postpro_tsv_simple2_amr2/port-S.csv`, `error-indicators.csv`  
`S11=-4.241096582075 dB, S21=-30.71097600287 dB, err norm=0.01397105348127`

### Papers (pcp11)
- [x] Regenerated comprehensive honest papers: `CurvilinearTSV/paper/main.tex/.pdf`, `comparison_report.tex/.pdf`.
- [x] Removed internal MCP/Wolfram references from publication artifacts.
- [x] Included audit, convergence trend, S21 discrepancy analysis, gds2mesh status, reproducibility artifacts.
- [x] Recompiled PDFs (main 9pp, comparison 4pp).

## Unresolved / Blocked

### P0 - Palace validated 60GHz
- [ ] Build boundary-meshing-safe 2-port model with correct trace-to-ground lumped port on same geometry/materials.
- [ ] Achieve err norm `<1e-2` with mesh/order/AMR convergence.
- [ ] Demonstrate passivity, reciprocity, energy/loss sanity at 60GHz.
- [ ] Resolve S21 discrepancy (`-30.71 dB` Palace vs `-2.39 dB` MTL) — currently attributed to port/model mismatch; S21 stable under refinement.
- [ ] Replace/retire invalid v2/v3 evidence in final validated comparison.

**Blocker:** OCC BooleanFragments drop/mis-mesh outer-boundary port-corridor faces; z-restricted port and stub extension produced no valid port surfaces. Needs boundary-meshing-safe CAD/mesh strategy (non-OCC-fragment workflow or manual mesh repair).

### P0 - Defensible 3D BEM/MTL on same validated full-wave geometry
- [ ] Implement defensible 3D BEM/MTL extraction on same geometry/materials/ports as validated Palace model (Maxwell C/R/L, PEEC, or RLGC->MTL->ABCD).
- [ ] Produce 60GHz cross-check table/plots with passivity/reciprocity checks.
- [ ] Include only validated comparison evidence in paper.

**Blocker:** Depends on validated Palace geometry/port definition. 2D quasi-static BEM exists but is not 3D full-wave same-geometry cross-check.

### P0 - Validated 60GHz cross-check
- [ ] Achieve geometry-matched, converged, err<1e-2 Palace-vs-Julia comparison.
- [ ] Add validated comparison evidence to `main.tex` and `comparison_report.tex`.

**Blocker:** Depends on P0 Palace validation + P0 3D BEM/MTL.

### P1 - gds2mesh polygon-vertex preservation complete
- [ ] Add OASIS XY/vertex parsing/storage in `src/db/oasis.cpp`.
- [ ] Complete polygon-preserving 3D extrusion flow (not just 2D polygon writer).
- [ ] Add robust Palace port-surface tagging workflow.
- [ ] Add tests for OASIS polygon round-trip / 3D polygon flow.

**Current state:** C API/Julia/GDSII vertex storage done; 34/34 tests pass; OASIS/full 3D flow incomplete.

## Key Files

### pcp11
- `CurvilinearTSV/paper/main.tex`, `main.pdf` — comprehensive paper (NOT PUBLICATION READY)
- `CurvilinearTSV/paper/comparison_report.tex`, `comparison_report.pdf` — adversarial audit
- `CurvilinearTSV/test/tsv_bem_corrected.jl` — corrected analytical RLC baseline
- `CurvilinearTSV/test/bem2d_capacitance.jl` — 2D Laplace BEM+MTL
- `CurvilinearTSV/test/tsv_mtl_simple2.jl` — simple2 MTL baseline
- `CurvilinearTSV/test/palace_tsv_simple2*.geo/msh/json` — Palace simple2 models
- `CurvilinearTSV/test/postpro_tsv_simple2_amr2/` — best AMR2 Palace result
- `CurvilinearTSV/test/postpro_tsv_2port_v2/` — invalid v2 result
- `CurvilinearTSV/test/tsv_sparameters_60ghz_corrected.csv` — corrected 4-port CSV

### VLSI_LAYOUT
- `src/db/layoutdb.h`, `layoutdb.cpp`, `layoutdb_internal.hpp` — C API, flat/polygon iteration
- `src/db/gdsii.cpp` — GDSII vertex storage
- `src/db/oasis.cpp` — **needs OASIS vertex parsing**
- `src/gds2mesh/src/GDS2Mesh.jl`, `layoutdb.jl`, `mesh2d.jl`, `mesh3d.jl`
- `src/gds2mesh/test/runtests.jl` — 34/34 pass

## Git State
- `pcp11`: main, recent commits through `0b9322f`/`8e2bbf6` (papers, AMR2, BEM, port attempts). Untracked local test artifacts present.
- `VLSI_LAYOUT`: main, `8b9e99e` (polygon C API/tests). Unrelated unstaged Fortran/db files present.

## Reproduction Commands
```bash
# gds2mesh tests
cd /Users/skoranne/GITHUB/VLSI_LAYOUT/src/gds2mesh
julia --project -e 'push!(LOAD_PATH,"src"); include("src/GDS2Mesh.jl"); using .GDS2Mesh; include("test/runtests.jl")'

# corrected analytical baseline
cd /Users/skoranne/GITHUB/pcp11
julia --project=CurvilinearTSV CurvilinearTSV/test/tsv_bem_corrected.jl

# 2D BEM+MTL
julia --project=CurvilinearTSV CurvilinearTSV/test/bem2d_capacitance.jl

# simple2 MTL baseline
julia --project=CurvilinearTSV CurvilinearTSV/test/tsv_mtl_simple2.jl

# compile papers
cd /Users/skoranne/GITHUB/pcp11/CurvilinearTSV/paper
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error comparison_report.tex
```
