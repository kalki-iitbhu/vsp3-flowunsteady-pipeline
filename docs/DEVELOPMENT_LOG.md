# Development & Validation Log — OpenVSP → FLOWUnsteady Pipeline

This document records, in order, everything done to build, debug, and validate the
`bem_vsp3_to_flowunsteady.py` pipeline and its accompanying comparison tools. For
each phase: what was attempted, where it went wrong (kept brief — these are dead
ends, not the point of the story), and what the actual correct fix was (explained
in full, since this is the part worth understanding and reusing).

---

## Phase 0 — Project Goal

Starting question: *does FLOWUnsteady use a 3D mesh for its propeller, and can we
get it?* This led to understanding that FLOWUnsteady (FLOW Lab, BYU) is a
**meshless** solver — a reformulated vortex particle method (rVPM) — so there is
no CFD volume mesh in the traditional sense. It does write one real geometric
VTK file per blade per timestep:

- `..._Rotor_Blade#_vlm_*.vtk` — for the default actuator-*line* model (used
  throughout this project), this is a **flat** representation: at each span
  station it stores just two points (leading edge, trailing edge), connected
  into quad panels along the span. It has no camber or thickness — it's a
  twisted, tapered flat ribbon, not a curved or solid surface. (An earlier
  draft of this log incorrectly stated a separate `_loft_*.vtk` "full-thickness"
  file also exists; this was never actually verified against a real file and
  was wrong — confirmed directly by listing a real run's output, which showed
  only `_vlm_` files. FLOWUnsteady's optional actuator-*surface* model,
  `vlm_vortexsheet=true`, would instead place many points across the chord
  following the true camber line — but this project never ran with that
  setting enabled, and switching it on changes the aerodynamic solve itself,
  not just the VTK output.)

This reframed the actual goal: **validate that FLOWUnsteady's geometry and
aerodynamics, generated from an OpenVSP model via this pipeline, actually match
OpenVSP's own geometry and VSPAERO's own aerodynamic prediction.**

---

## Phase 1 — Understanding FLOWUnsteady's Internals

Before touching any code, the following was established by reading FLOWUnsteady's
own example scripts and documentation:

- The six-stage structure every FLOWUnsteady script follows: Vehicle definition →
  Maneuver definition → Simulation definition → Monitors → `run_simulation()` →
  Output files.
- `generate_rotor()` builds the blade from `chorddist`/`twistdist`/`sweepdist`/
  `heightdist` tables (each a function of `r/R`), lofts the airfoil sections, and
  (optionally) calls XFOIL per station to get real lift/drag polars.
- The wake is built from vortex particles shed at each blade element boundary,
  driven by two mechanisms: a **trailing vortex** (from the spanwise circulation
  gradient between adjacent elements) and an **unsteady/starting vortex** (from
  circulation changing in time, governed by `unsteady_shedcrit`).
- The default actuator-line model represents each blade element as a flat
  leading-edge-to-trailing-edge panel (no camber curvature in the VTK); the
  actuator-*surface* model (`vlm_vortexsheet=true`) instead builds a full
  chordwise lattice with camber, using the mean-camber-line — the same "mid of
  upper and lower surface" construction OpenVSP's own DegenGeom `PLATE` export
  uses, which became the basis for the geometry comparison.

---

## Phase 2 — Building the Geometry Comparison Pipeline

**Approach chosen**: rather than deep-diving into FLOWUnsteady's internal spline
math, compare its *output geometry* directly against OpenVSP's own geometry
export (DegenGeom), since both ultimately produce a mean-camber surface that can
be extracted as real 3D points.

- **`openvsp_prop_to_vtk.py`** — parses an OpenVSP DegenGeom CSV's `PLATE`
  section, pulling the `xxCamber/xyCamber/xzCamber` columns (the actual 3D
  camber-surface location, not the flattened `x,y,z` reference point also
  present in the file) and writes a VTK `STRUCTURED_GRID` mesh per blade.
- **`compare_vsp_flowunsteady_vtk.py`** — parses both the OpenVSP mesh and
  FLOWUnsteady's `_vlm_*.vtk` file generically (no hardcoded station counts),
  auto-detects which coordinate axis is "radius" in each file (largest
  coordinate range), reduces each cross-section to a chord vector, and compares
  chord length/twist angle at matching span stations.

**❌ Failure point (concise)**: a parser bug assumed the FLOWUnsteady VTK's two
edge-point arrays were each `n_span - 1` elements, when they were actually
`n_span` each — an off-by-one that made every computed chord length come out
constant and wrong. Found immediately by sanity-checking against hand-traced
values; fixed by correcting the index arithmetic in `parse_fu_vlm_vtk()`.

**✅ First real result**: with the parser fixed, radius range matched exactly
(0.03–0.15 m in both files) and chord length agreed to ~1–2% almost everywhere —
strong early evidence the two tools describe the same propeller. Twist appeared
to mismatch wildly station-by-station until it was noticed that
`angle_FU = -180° - angle_VSP` held *exactly* at every single station — a
constant relationship, which is the signature of a coordinate-frame reflection
difference, not a real twist disagreement. Once corrected for, twist matched
too. The one genuine difference found at this stage: OpenVSP's pitch axis sits
at **mid-chord**, while FLOWUnsteady's sits at the **leading edge** — a labeling
convention difference, not a shape difference (chord midpoint was exactly `(0,0)`
in one file, LE was exactly `(0,0)` in the other, every station, with the
midpoint offset in the FU file equal to exactly half the local chord length).

---

## Phase 3 — Fixing the Alignment Methodology

To overlay the two meshes visually in ParaView, the OpenVSP mesh needed
rotating/reflecting into FLOWUnsteady's coordinate frame.

**❌ Failure point (concise)**: the first alignment implementation rebased every
span station's points *relative to that station's own leading-edge point*
before rotating. Since this subtracts a point from itself for the LE, the LE
always mapped to exactly `(0,0)` by construction — producing a suspiciously
**perfectly straight leading edge** in the aligned output, regardless of the
blade's real shape. This was caught by direct visual inspection in ParaView
("one of the edges is a straight line").

**✅ Correct fix (detail)**: the per-station rebase was removed entirely. OpenVSP's
raw camber-surface coordinates are *already* expressed relative to the true
rotation axis (that's precisely why the mid-chord line came out at `(0,0)` for
every station in the first place — a real geometric fact about this propeller's
feather-axis placement, not an artifact). The correct alignment is therefore a
**pure rotation + reflection + scale**, applied directly to the raw coordinates,
with *no translation removed*. After this fix, the aligned mesh's leading edge
showed genuine twist-induced curvature (varying magnitude and direction station
to station, e.g. moving from `(-0.0043, 0.0041)` at the hub to `(-0.0072,
0.0121)` mid-span), while the mid-chord line correctly stayed near `(0,0)`
throughout — the geometrically honest result.

---

## Phase 4 — Diagnosing and Fixing the Conversion Script

With a trustworthy comparison tool in hand, attention turned to
`bem_vsp3_to_flowunsteady.py` itself (the actual geometry conversion pipeline).

### Bug 1: Feather-axis reference mismatch

OpenVSP's `.bem` export gives `Rake`/`Skew` columns described in OpenVSP's own
documentation as relative to the local chord orientation (i.e. referenced to
the propeller's **feather axis**, a user-set chord fraction — confirmed directly
in the `.vsp3` XML as `<FeatherAxisXoC Value="0.5"/>`, meaning mid-chord for this
propeller). The pipeline was copying these values verbatim into FLOWUnsteady's
`sweepdist`/`heightdist`, which `generate_rotor` interprets as offsets of the
**leading edge**. Since `Rake`/`Skew` were all exactly `0.0` in this propeller's
table, the bug's effect was to force FLOWUnsteady's leading edge onto the
rotation axis — exactly reproducing the straight-LE symptom independently
diagnosed in Phase 3, now traced to its root cause in the actual pipeline code.

**✅ Fix**: `apply_feather_axis_correction()` was added, converting the
feather-axis-referenced offset into a leading-edge-referenced one using the
known feather-axis fraction and the local twist angle (rotating the
chord-relative offset into the global sweep/height frame).

**❌ Failure point (concise)**: the first version of this correction had the
*sign backwards*. Verified empirically once a real FLOWUnsteady VLM output was
available: the resulting leading-edge point was found to sit at almost exactly
`+0.5 × (chord vector)` at every station (confirmed to 4 significant figures
across multiple stations) — the mirror image of the intended `-0.5 × chord`.
Flipping the sign (`+le_to_feather` instead of `-le_to_feather`) resolved it;
re-verified afterward to show `FU pitch-axis convention: mid-chord (50%)`,
correctly matching OpenVSP's own convention.

### Bug 2: Julia float/int formatting

The `_fmt()` helper used Python's `%.10g` format, which strips the decimal point
from whole numbers (`1.0 → "1"`). Written into the generated `.jl` script, Julia
parses `1` as an `Int64`, not `Float64` — causing
`TypeError: ... expected Float64, got Int64` the first time a CLI option
(`--blade-r 1`) happened to produce an exact whole number.

**✅ Fix**: `_fmt()` now checks for the absence of a decimal point/exponent and
appends `.0` — guaranteeing every numeric value written into the Julia script is
unambiguously a `Float64` literal.

### Bug 3: Chord-curve parsing failure (wrong Reynolds numbers)

A separate code path, used only to estimate each XFOIL polar's Reynolds number,
searched for a `<ParmContainer Name="Chord">` spline-curve block in the `.vsp3`
that simply does not exist in this file's layout — silently falling back to a
flat `chord/R = 0.15` assumption for *every* station, when the real chord
ranges from `0.012–0.03 m` (5–12× smaller). This meant every XFOIL polar in
every prior run was generated at a drastically overestimated Reynolds number.

**✅ Fix**: the real per-station `Chord` value is read directly from each XSec
block in the `.vsp3` (the same block already being parsed for `Camber`,
`ThickChord`, etc.) — no curve reconstruction needed, and confirmed to exactly
match the true tip chord independently established from the `.bem` table
(`0.0195 m` at `r/R = 1.0`).

### Bug 4: Unwanted chord-table smoothing

`generate_rotor` fits a smooth curve through the `.bem` table's chord/twist/
sweep/height points (as a function of span position), then reads values off
that curve for each of the `n` blade-element stations. `spline_s` controls how
tightly that fit is forced through the real data versus how much it's allowed
to round off sharp features for smoothness — and the default is nonzero.

**This fit, and therefore this bug, applies along the entire span, to all four
distributions — it is not inherently a "tip-only" issue.** It only showed up
as a tip-region error in this propeller's specific comparison because the mid-
span chord values are gently varying (nothing sharp for smoothing to distort),
while the tip is where this blade's real taper gets steep — the one place a
smoothing spline actually had something sharp to round off. A blade with a
sharp feature elsewhere along its span (e.g. a sudden sweep change near the
root) would show this same bug there instead.

**❌ Failure point (concise, and a dead end worth recording)**: increasing `n`
(blade-element count) was expected to *reduce* this tip error through better
resolution. Instead, the tip-region error **grew** with more elements (`n=20`:
8.5% off; `n=40`: 16.5%; `n=80`: 21.2%) — the opposite of normal convergence
behavior. A side hypothesis (element clustering via `blade_r` concentrating
resolution badly near the tip) was tested and **ruled out**: `n=40` with
clustered spacing and `n=40` with perfectly uniform spacing produced the
*identical* tip chord value to 4 decimal places, proving `blade_r` was not the
cause. (A second false lead — suspecting a literal zero-chord tip point — was
also chased and found to be based on a misread of the airfoil-contour file
rather than the real blade table, which actually has a finite tip chord of
`0.13 × R`.)

**✅ Correct fix**: setting `spline_s=0.0` (forcing the fit to pass exactly
through every real table point, no smoothing) was tested directly and produced
a **perfect match** against OpenVSP — `0.00%` difference at every single
station, fit residual at machine precision (`~1e-10`). This confirmed smoothing,
not clustering or a geometric singularity, was the actual cause, and it is now
the pipeline's default.

**Outcome of Phase 4**: with all four bugs fixed, FLOWUnsteady's rotor geometry
matches OpenVSP's DegenGeom export to within numerical precision.

---

## Phase 5 — OpenVSP Mesh Convergence Study

Before fully trusting the OpenVSP side of the comparison, a convergence check
was built to confirm the DegenGeom export itself wasn't under-resolved.

- **`openvsp_convergence_sweep.py`** drives OpenVSP's own Python API to sweep
  spanwise tessellation (`Tess_U`) and export a DegenGeom CSV at each level.
- **`convergence_analysis.py`** tracks a chosen metric (chord length at a given
  `r/R`) across the sweep and reports when it stops changing.

**❌ Failure point (concise)**: the first version of the sweep script used
OpenVSP's Analysis-Manager API (`SetStringAnalysisInput("DegenGeom",
"CSVFileName", ...)`) to set the output filename — this input simply doesn't
exist on the user's OpenVSP version (`GetAnalysisInputNames` confirmed only
`WriteCSVFlag` was present, no filename control).

**✅ Fix**: switched to the lower-level, version-stable API pair
`vsp.SetComputationFileName(vsp.DEGEN_GEOM_CSV_TYPE, path)` +
`vsp.ComputeDegenGeom(vsp.SET_ALL, vsp.DEGEN_GEOM_CSV_TYPE)`, confirmed against
OpenVSP's own documented example usage.

**Result**: chord length at `r/R = 0.95` converged to under 0.5% change by
`Tess_U ≈ 19`, with a brief overshoot around `Tess_U = 17` before settling —
`Tess_U ≈ 30–40` was adopted as a safe working value with margin.

---

## Phase 6 — Aerodynamic Cross-Validation Against VSPAERO

With geometry validated, the comparison moved to aerodynamics: does
FLOWUnsteady's predicted thrust/torque/efficiency agree with VSPAERO's, for the
same propeller and operating point (`J = 0.4`, 9,200 RPM, `V∞ = 18.4 m/s`)?

**❌ Failure point (concise)**: VSPAERO's raw `History` CSV output stores results
**per blade panel, per solver iteration** — a label like `CT_h` appears hundreds
of times in the file, almost all of them per-panel breakdowns. The first attempt
blindly matched every occurrence of each label, producing a huge, meaningless
dump. The real whole-rotor aggregate time-history for each quantity turned out
to be specifically the **final occurrence** of that label in the file (with a
row length matching the true iteration count, `200`, rather than the per-panel
block length of `183`) — found by cross-checking row lengths and confirming the
resulting advance ratio (`J = 0.4`) matched the known operating condition
exactly.

**First comparison (VSPAERO defaults)**: CT off by −16.7%, CQ off by −22.1%,
efficiency off by +7.0%.

---

## Phase 7 — The Reynolds-Number Fix, Round Two

Even after Bug 3 (Phase 4) was fixed for the *geometry* table, a **second,
independent instance of the same root problem** was found: a separate parsing
path specifically for XFOIL's Reynolds-number estimate was still falling back
to a flat, wrong chord assumption, because it had its own copy of the broken
`ParmContainer` search. This was fixed the same way — reading the real
per-station chord directly.

**Result after this fix**: CQ moved substantially closer to VSPAERO (−22.1% →
−9.2%) and efficiency improved dramatically (+7.0% → −9.0%), confirming the
Reynolds-number fix mattered. **CT, however, moved further away** (−16.7% →
−17.4%) — a genuine, physically explainable side effect: the earlier (wrong)
flat chord assumption had overestimated Reynolds number by 5–12×; correcting it
revealed real low-Reynolds-number airfoil behavior (reduced lift, not just
increased drag), which VSPAERO's own simplified viscous model does not
necessarily reproduce the same way.

---

## Phase 8 — Correcting VSPAERO's Own Settings

Investigating the remaining CT gap led to inspecting VSPAERO's own case-setup
file directly, which surfaced two more real issues — this time on the
*reference* side of the comparison, not in the pipeline being validated.

- **`ReCref = 10,000,000`** — VSPAERO's reference Reynolds number for its
  viscous/stall correction was a generic default, roughly **45× higher** than
  this propeller's real operating Reynolds number (computed properly, using the
  reference chord and blade speed at 75% radius, as `≈ 220,000`).
- **`Clo2D = 0`** — the zero-angle-of-attack lift coefficient used in VSPAERO's
  viscous drag estimate was left at zero, ignoring the real cambered airfoil's
  actual `Cl(α=0) ≈ 0.08` (measured directly from the same XFOIL polars already
  generated for FLOWUnsteady).
- Initial confusion over `StallModel = 0` was resolved by checking the OpenVSP
  GUI directly: `0` means the stall model is **On** (a different enum
  convention than assumed), not disabled — a mistake corrected once the
  screenshot was reviewed, rather than left standing.

**✅ Result after correcting VSPAERO's settings**: VSPAERO's own predicted CT
dropped 27% (0.2001 → 0.1460) once given a realistic Reynolds number and real
camber data — strong confirmation the original comparison was unfair to begin
with, not that FLOWUnsteady was simply wrong.

---

## Phase 9 — Final Validated Comparison

| Quantity | VSPAERO (corrected settings) | FLOWUnsteady | % diff |
|---|---|---|---|
| CT | 0.1460 | 0.1654 | 13.2% |
| CQ | 0.0190 | 0.0209 | 10.2% |
| **Efficiency (η)** | **0.4901** | **0.5037** | **2.8%** |

Geometry: matches to numerical precision. Efficiency: within 2.8%. The
remaining 10–13% gap on CT/CQ individually is attributed to real fidelity
differences between FLOWUnsteady's blade-element + real-polar + free-wake
model and VSPAERO's VLM + empirical global stall clamp — not a remaining bug,
though the true `CLmax` (the XFOIL sweep's current +12° cutoff never actually
reached stall) remains a loose end for a tighter future check.

---

## Phase 10 — Repository Organization

Final cleanup pass before publishing the validated pipeline:

**❌ Failure points (concise, all resolved)**:
- A 401 MB VSPAERO solver cache file (`.adb`) was sitting at the repo root —
  GitHub hard-rejects any file over 100 MB; confirmed via `git ls-files` that it
  had never actually been committed, so no history rewrite was needed, just
  exclusion via `.gitignore` (which, on inspection, was already correctly
  configured to exclude it and every other large regeneratable solver file).
- Several PowerShell sessions were disrupted by pasting old terminal
  **scrollback output** back in as if it were a new command (producing a string
  of `"X" is not recognized` errors) — not a real tooling problem, just
  terminal-history confusion, resolved by re-running single clean commands.
- A `git push` was rejected (`fetch first`) because the remote had an unrelated
  README change — resolved with a normal `git pull --no-rebase` merge, no
  conflict.

**✅ Final structure**: the repository was reorganized into `tools/` (the four
reusable comparison/validation scripts), `examples/propellor02/{openvsp,
flowunsteady, paraview}` (cleanly separating each tool's native outputs), and
the heavy, fully-regeneratable raw simulation run folder renamed and left local
(already gitignored). The README was rewritten to document the full
methodology end-to-end, including an interactive `<model-viewer>`-based 3D
comparison page (hosted via GitHub Pages) overlaying the OpenVSP and
FLOWUnsteady meshes directly in-browser.

---

## Summary: Bugs Fixed vs. Dead Ends Explored

| # | Issue | Category |
|---|---|---|
| 1 | Feather-axis reference mismatch in sweep/height | **Real bug, fixed** |
| 2 | Sign error in the Bug-1 fix itself | **Real bug, fixed** (caught before shipping) |
| 3 | Julia Int64/Float64 formatting bug | **Real bug, fixed** |
| 4 | Chord-curve parsing fallback (geometry table) | **Real bug, fixed** |
| 5 | Chord-curve parsing fallback (XFOIL Re estimate, separate instance) | **Real bug, fixed** |
| 6 | Unwanted spline smoothing distorting tip chord | **Real bug, fixed** |
| 7 | Per-station rebase destroying real LE curvature (comparison tool, not pipeline) | **Real bug, fixed** |
| 8 | Off-by-one in VLM VTK parser (comparison tool) | **Real bug, fixed** |
| 9 | Wrong OpenVSP API for DegenGeom filename | **Real bug, fixed** |
| 10 | VSPAERO `ReCref`/`Clo2D` defaults unrealistic for this propeller | **Real bug (in the reference tool), fixed** |
| — | Element clustering (`blade_r`) as tip-error cause | Dead end — ruled out empirically |
| — | Literal zero-chord tip point as cause | Dead end — based on a misread file |
| — | `StallModel=0` assumed to mean "disabled" | Dead end — wrong enum assumption, corrected from GUI |
| — | A separate `_loft_*.vtk` "full-thickness" file assumed to exist | Dead end — never verified, corrected by listing real simulation output (only `_vlm_` files exist) |
