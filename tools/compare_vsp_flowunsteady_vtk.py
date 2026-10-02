#!/usr/bin/env python3
"""
compare_vsp_flowunsteady_vtk.py
---------------------------------
Compare a propeller blade geometry from OpenVSP (camber-surface VTK, as
produced by openvsp_prop_to_vtk.py) against FLOWUnsteady's own blade VLM
VTK output, without assuming they share coordinate convention, chord
reference point ("pitch axis"), or span-station discretization.

REQUIRED INPUT FILES
----------------------
1. --vsp   A camber-surface VTK STRUCTURED_GRID, e.g. "camber_surface_blade0.vtk"
           (produced by openvsp_prop_to_vtk.py from an OpenVSP DegenGeom export).
2. --fu    A FLOWUnsteady rotor VLM VTK, e.g.
           "<run_name>_Rotor_Blade1_vlm_<step>.vtk", written automatically
           whenever you run/generate a rotor with save_path set in
           FLOWUnsteady (uns.generate_rotor + uns.run_simulation, or any
           postprocessing step that calls FLOWVLM's VTK writer).

WHAT THIS SCRIPT DOES
------------------------
1. Parses both files generically (auto-detects grid sizes; does not assume
   a fixed number of span or chord stations).
2. Auto-detects which coordinate axis is the "radius" (span) axis in each
   file, as the axis with the largest coordinate range.
3. Reduces each blade cross-section to a chord vector (leading-edge point
   to trailing-edge point) at every span station, from which it computes
   chord length, chordwise angle, and reference-point ("pitch axis")
   location.
4. Detects whether each file's stored coordinates are referenced to the
   leading edge or to the mid-chord (or something in between), since the
   two tools are not guaranteed to agree on this.
5. Interpolates the finer-resolution dataset onto the coarser dataset's
   radial stations, and reports chord-length agreement.
6. Solves for the best-fit rotation + reflection relating the two files'
   local chordwise coordinate planes (since VSP and FLOWUnsteady each pick
   their own local in-plane axes), using a simple least-squares fit.
7. Applies that fit to re-express the FULL OpenVSP camber surface (all
   chordwise points, not just the two edges) in FLOWUnsteady's frame, and
   writes it out as an aligned VTK you can overlay directly in ParaView.

USAGE
-----
    python3 compare_vsp_flowunsteady_vtk.py \
        --vsp camber_surface_blade0.vtk \
        --fu propeller_Rotor_Blade1_vlm_0.vtk \
        --outdir ./out
"""

import argparse
import os
import numpy as np


# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------

def parse_vsp_structured_grid(path):
    """Parse a VTK legacy STRUCTURED_GRID (as written by
    openvsp_prop_to_vtk.py): DIMENSIONS npn nxsecs 1, points chordwise-fastest.
    Returns array of shape (nxsecs, npn, 3).
    """
    with open(path) as f:
        lines = f.readlines()

    dims_line = next(l for l in lines if l.startswith("DIMENSIONS"))
    npn, nxsecs, _ = (int(v) for v in dims_line.split()[1:4])

    pts_idx = next(i for i, l in enumerate(lines) if l.startswith("POINTS"))
    npts = int(lines[pts_idx].split()[1])
    pts = np.array([[float(v) for v in lines[pts_idx + 1 + i].split()]
                     for i in range(npts)])

    if npts != nxsecs * npn:
        raise ValueError(f"Point count {npts} != nxsecs*npn ({nxsecs}*{npn})")

    return pts.reshape(nxsecs, npn, 3)


def parse_fu_vlm_vtk(path):
    """Parse a FLOWUnsteady rotor VLM VTK (UNSTRUCTURED_GRID of quad panels
    forming a single chordwise strip: two parallel edge lines connected by
    quads). Returns (edgeA, edgeB), each shape (n_span, 3), running from
    hub to tip. Ignores any extra VTK_VERTEX cells / CELL_DATA (aero results).
    """
    with open(path) as f:
        lines = f.readlines()

    pts_idx = next(i for i, l in enumerate(lines) if l.startswith("POINTS"))
    npts = int(lines[pts_idx].split()[1])
    pts = np.array([[float(v) for v in lines[pts_idx + 1 + i].split()]
                     for i in range(npts)])

    cells_idx = next(i for i, l in enumerate(lines) if l.startswith("CELLS"))
    n_cells = int(lines[cells_idx].split()[1])
    quads = []
    i = cells_idx + 1
    count = 0
    while count < n_cells:
        parts = [int(v) for v in lines[i].split()]
        if parts[0] == 4:
            quads.append(parts[1:5])
        i += 1
        count += 1

    if not quads:
        raise ValueError(
            "No quad cells found in FU VTK - is this a VLM lattice file "
            "with a bound-vortex/trailing-edge panel strip?"
        )

    # Expected pattern from FLOWVLM's writer: quad k = [i, i+n, i+n+1, i+1]
    n_line = quads[0][1] - quads[0][0]
    for q in quads:
        if q[1] - q[0] != n_line or q[2] - q[3] != n_line:
            raise ValueError(
                "Unexpected quad connectivity pattern - this parser assumes "
                "a single chordwise strip (actuator-line VLM output). "
                "Multi-row actuator-surface VTKs need a different parser."
            )

    n_span = n_line
    edgeA = pts[0:n_span]
    edgeB = pts[n_span:2 * n_span]
    return edgeA, edgeB


# --------------------------------------------------------------------------
# Geometry reduction helpers
# --------------------------------------------------------------------------

def detect_radius_axis(points_flat):
    """points_flat: (N,3) array. Returns axis index (0,1,2) with largest range."""
    ranges = points_flat.max(axis=0) - points_flat.min(axis=0)
    return int(np.argmax(ranges))


def chord_vectors_vsp(grid, radius_axis):
    """grid: (nxsecs, npn, 3). Returns radius (nxsecs,), start pts, end pts,
    both (nxsecs,3), using the first and last chordwise points as the
    two blade edges."""
    radius = grid[:, 0, radius_axis]
    start = grid[:, 0, :]
    end = grid[:, -1, :]
    return radius, start, end


def in_plane_components(points, radius_axis):
    """Drop the radius axis, return the remaining two coordinates in their
    original column order."""
    axes = [a for a in range(3) if a != radius_axis]
    return points[:, axes]


def anchor_offset_pattern(radius_axis_start, mid, start, radius_axis):
    """Return the mean magnitude of `start` and of `mid`, projected onto the
    in-plane axes, relative to the radius axis line (x=y=0 in-plane) -
    used to detect whether a file is LE-anchored (start ~ 0) or
    mid-chord-anchored (mid ~ 0)."""
    axes = [a for a in range(3) if a != radius_axis]
    start_mag = np.linalg.norm(start[:, axes], axis=1).mean()
    mid_mag = np.linalg.norm(mid[:, axes], axis=1).mean()
    return start_mag, mid_mag


def fit_rotation_reflection(c_from, c_to):
    """Given two arrays of complex numbers (in-plane chord vectors) sampled
    at matching stations, find the best-fit complex multiplier m (rotation
    + scale) such that c_to ~ m * c_from, trying both direct and
    conjugated (reflected) versions of c_from. Returns (m, reflect, resid)."""
    def fit(c_src):
        denom = np.sum(np.abs(c_src) ** 2)
        m = np.sum(c_to * np.conj(c_src)) / denom
        resid = np.sum(np.abs(c_to - m * c_src) ** 2)
        return m, resid

    m_direct, r_direct = fit(c_from)
    m_reflect, r_reflect = fit(np.conj(c_from))

    if r_reflect < r_direct:
        return m_reflect, True, r_reflect
    return m_direct, False, r_direct


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vsp", required=True, help="OpenVSP camber-surface VTK "
                                                    "(structured grid)")
    ap.add_argument("--fu", required=True, help="FLOWUnsteady rotor VLM VTK")
    ap.add_argument("--outdir", default=".", help="Where to write the aligned "
                                                     "VTK output")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ---- parse ----
    vsp_grid = parse_vsp_structured_grid(args.vsp)          # (nxsecs, npn, 3)
    fu_edgeA, fu_edgeB = parse_fu_vlm_vtk(args.fu)           # each (n_span, 3)

    # ---- detect radius axis for each file independently ----
    vsp_axis = detect_radius_axis(vsp_grid.reshape(-1, 3))
    fu_axis = detect_radius_axis(np.vstack([fu_edgeA, fu_edgeB]))
    print(f"Detected radius axis: VSP -> {'xyz'[vsp_axis]}, "
          f"FU -> {'xyz'[fu_axis]}")

    vsp_radius, vsp_start, vsp_end = chord_vectors_vsp(vsp_grid, vsp_axis)
    fu_radius = fu_edgeA[:, fu_axis]
    fu_start, fu_end = fu_edgeA, fu_edgeB

    vsp_mid = 0.5 * (vsp_start + vsp_end)
    fu_mid = 0.5 * (fu_start + fu_end)

    print(f"Radius range: VSP [{vsp_radius.min():.4f}, {vsp_radius.max():.4f}]  "
          f"FU [{fu_radius.min():.4f}, {fu_radius.max():.4f}]")

    # ---- detect pitch-axis anchor convention ----
    vsp_start_mag, vsp_mid_mag = anchor_offset_pattern(None, vsp_mid, vsp_start, vsp_axis)
    fu_start_mag, fu_mid_mag = anchor_offset_pattern(None, fu_mid, fu_start, fu_axis)

    def describe_anchor(start_mag, mid_mag):
        if mid_mag < 0.1 * start_mag:
            return "mid-chord (50%)"
        elif start_mag < 0.1 * mid_mag:
            return "leading edge (0%)"
        else:
            return f"ambiguous (start_offset={start_mag:.5f}, mid_offset={mid_mag:.5f})"

    print(f"VSP pitch-axis convention: {describe_anchor(vsp_start_mag, vsp_mid_mag)}")
    print(f"FU  pitch-axis convention: {describe_anchor(fu_start_mag, fu_mid_mag)}")

    # ---- interpolate the finer dataset onto the coarser one's radii ----
    if len(vsp_radius) <= len(fu_radius):
        common_radius = vsp_radius
        vsp_start_c, vsp_end_c = vsp_start, vsp_end
        fu_start_c = np.stack([np.interp(common_radius, fu_radius, fu_start[:, a])
                                for a in range(3)], axis=1)
        fu_end_c = np.stack([np.interp(common_radius, fu_radius, fu_end[:, a])
                              for a in range(3)], axis=1)
    else:
        common_radius = fu_radius
        fu_start_c, fu_end_c = fu_start, fu_end
        vsp_start_c = np.stack([np.interp(common_radius, vsp_radius, vsp_start[:, a])
                                 for a in range(3)], axis=1)
        vsp_end_c = np.stack([np.interp(common_radius, vsp_radius, vsp_end[:, a])
                               for a in range(3)], axis=1)

    vsp_chord_c = vsp_end_c - vsp_start_c
    fu_chord_c = fu_end_c - fu_start_c
    vsp_len = np.linalg.norm(in_plane_components(vsp_chord_c, vsp_axis), axis=1)
    fu_len = np.linalg.norm(in_plane_components(fu_chord_c, fu_axis), axis=1)

    pct_diff = 100 * (fu_len - vsp_len) / vsp_len

    print("\n--- Chord length comparison (common radii) ---")
    print(f"{'radius':>10} {'VSP chord':>12} {'FU chord':>12} {'% diff':>8}")
    for r, a, b, d in zip(common_radius, vsp_len, fu_len, pct_diff):
        print(f"{r:10.4f} {a:12.5f} {b:12.5f} {d:8.2f}")
    print(f"\nmax |%diff| = {np.max(np.abs(pct_diff)):.2f}   "
          f"mean |%diff| = {np.mean(np.abs(pct_diff)):.2f}   "
          f"rms %diff = {np.sqrt(np.mean(pct_diff**2)):.2f}")

    # ---- solve rotation/reflection between the two in-plane conventions ----
    vsp_axes = [a for a in range(3) if a != vsp_axis]
    fu_axes = [a for a in range(3) if a != fu_axis]
    c_vsp = vsp_chord_c[:, vsp_axes[0]] + 1j * vsp_chord_c[:, vsp_axes[1]]
    c_fu = fu_chord_c[:, fu_axes[0]] + 1j * fu_chord_c[:, fu_axes[1]]

    m, reflected, resid = fit_rotation_reflection(c_vsp, c_fu)
    angle_deg = np.degrees(np.angle(m))
    scale = np.abs(m)
    print(f"\nBest-fit transform (VSP in-plane -> FU in-plane): "
          f"{'reflection + ' if reflected else ''}rotation {angle_deg:.2f} deg, "
          f"scale {scale:.4f}  (fit residual {resid:.3e})")

    # ---- apply the fit to the FULL vsp camber grid, writing it out in FU's
    #      frame. IMPORTANT: we do NOT rebase each station to its own
    #      leading edge here. OpenVSP's raw in-plane coordinates are already
    #      expressed relative to the true rotation axis (that's exactly why
    #      its mid-chord line came out at (0,0) for every station - a real
    #      geometric fact, not an artifact). A per-station rebase-to-LE
    #      would force the LE to (0,0) at every radius too, erasing the real
    #      twist-induced curvature of the LE locus. So we only rotate/
    #      reflect/scale the raw points - no translation is removed.
    nx, npn, _ = vsp_grid.shape
    aligned = np.zeros((nx, npn, 3))
    for j in range(nx):
        for k in range(npn):
            p = vsp_grid[j, k, :]
            c = p[vsp_axes[0]] + 1j * p[vsp_axes[1]]
            c_t = np.conj(c) if reflected else c
            c_new = m * c_t
            new_pt = np.zeros(3)
            new_pt[fu_axes[0]] = c_new.real
            new_pt[fu_axes[1]] = c_new.imag
            new_pt[fu_axis] = vsp_radius[j]
            aligned[j, k, :] = new_pt

    out_path = os.path.join(args.outdir, "aligned_vsp_to_fu.vtk")
    with open(out_path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("OpenVSP camber surface aligned into FLOWUnsteady frame\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {npn} {nx} 1\n")
        f.write(f"POINTS {nx*npn} float\n")
        for j in range(nx):
            for k in range(npn):
                x, y, z = aligned[j, k, :]
                f.write(f"{x:.9e} {y:.9e} {z:.9e}\n")
    print(f"\nWrote aligned overlay mesh: {out_path}")
    print("Load this alongside your original FU VLM VTK in ParaView to "
          "visually compare.")


if __name__ == "__main__":
    main()
