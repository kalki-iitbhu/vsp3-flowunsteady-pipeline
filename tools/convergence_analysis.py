#!/usr/bin/env python3
"""
convergence_analysis.py
--------------------------
Read the DegenGeom CSVs produced by openvsp_convergence_sweep.py (one per
spanwise-tessellation level) and check whether the extracted camber-surface
geometry has converged - i.e. whether increasing tessellation further stops
changing the answer.

Does NOT require OpenVSP - just re-uses the PLATE-section parsing logic
from openvsp_prop_to_vtk.py, so keep that file in the same directory (or
adjust the import below to point at it).

METRIC
------
By default, tracks chord length at a chosen radius fraction (r/R) of the
blade, interpolated from each export's own station grid, plus its percent
change from the previous (next-coarser) level - this is the standard
grid-convergence-index style check: watch the answer stop moving.

Default r/R = 0.95 (near the tip, where convergence usually matters most
for a tapered/twisted blade) - override with --rR if you want a different
station, or --metric le-curve-length for a whole-blade shape metric instead
of a single-station one.

USAGE
-----
    python3 convergence_analysis.py --manifest ./sweep/manifest.csv --rR 0.95
    python3 convergence_analysis.py --manifest ./sweep/manifest.csv --metric le-curve-length
"""

import argparse
import csv
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from openvsp_prop_to_vtk import parse_plate_sections  # noqa: E402


def chord_length_at_rR(blade, target_rR, R_tip):
    """Interpolate chord length (LE-to-TE distance) at a given r/R for one
    blade's PLATE data."""
    nx, npn = blade["nxsecs"], blade["npnts"]
    pts = blade["points"]
    radii = np.array([pts[j * npn][col] for j in range(nx)
                       for col in ["xyCamber"]])  # assumes y = radius, as in
                                                    # our earlier exports; if
                                                    # your model's radius axis
                                                    # differs, adjust here
    chordlens = []
    for j in range(nx):
        le = np.array([pts[j * npn][c] for c in ("xxCamber", "xyCamber", "xzCamber")])
        te = np.array([pts[j * npn + npn - 1][c] for c in ("xxCamber", "xyCamber", "xzCamber")])
        chordlens.append(np.linalg.norm(te - le))
    chordlens = np.array(chordlens)

    order = np.argsort(radii)
    radii_sorted = radii[order]
    chordlens_sorted = chordlens[order]

    target_r = target_rR * R_tip
    return float(np.interp(target_r, radii_sorted, chordlens_sorted))


def le_curve_length(blade):
    """Whole-blade shape metric: total arc length of the leading-edge locus
    across span. Sensitive to both discretization AND real geometric
    fidelity (sweep/twist wobble), unlike a single-station chord length."""
    nx, npn = blade["nxsecs"], blade["npnts"]
    pts = blade["points"]
    le_pts = np.array([[pts[j * npn][c] for c in ("xxCamber", "xyCamber", "xzCamber")]
                        for j in range(nx)])
    radii = le_pts[:, 1]  # assumes y = radius; adjust if needed
    order = np.argsort(radii)
    le_pts = le_pts[order]
    seg_lengths = np.linalg.norm(np.diff(le_pts, axis=0), axis=1)
    return float(np.sum(seg_lengths))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True,
                     help="manifest.csv written by openvsp_convergence_sweep.py "
                          "(columns: tess_u, csv_path)")
    ap.add_argument("--rR", type=float, default=0.95,
                     help="r/R station to track chord length at (used unless "
                          "--metric le-curve-length is given)")
    ap.add_argument("--metric", choices=["chord-at-rR", "le-curve-length"],
                     default="chord-at-rR")
    ap.add_argument("--tol-pct", type=float, default=0.5,
                     help="Convergence tolerance in percent change between "
                          "consecutive levels")
    ap.add_argument("--blade", type=int, default=0,
                     help="Which blade/PLATE section to use if the export "
                          "has multiple blades")
    args = ap.parse_args()

    levels = []
    with open(args.manifest) as f:
        for row in csv.DictReader(f):
            levels.append((int(row["tess_u"]), row["csv_path"]))
    levels.sort()

    results = []
    for tess, path in levels:
        blades = parse_plate_sections(path)
        blade = blades[args.blade]

        nx, npn = blade["nxsecs"], blade["npnts"]
        pts = blade["points"]
        radii = np.array([pts[j * npn]["xyCamber"] for j in range(nx)])
        R_tip = radii.max()

        if args.metric == "chord-at-rR":
            value = chord_length_at_rR(blade, args.rR, R_tip)
        else:
            value = le_curve_length(blade)

        results.append((tess, nx, value))
        print(f"Tess_U={tess:4d}  n_stations={nx:3d}  "
              f"{args.metric}={value:.6f}")

    print("\n--- convergence check (percent change from previous level) ---")
    converged_at = None
    for i in range(1, len(results)):
        tess, nx, value = results[i]
        _, _, prev_value = results[i - 1]
        pct_change = 100 * abs(value - prev_value) / abs(prev_value) if prev_value else float("nan")
        flag = ""
        if pct_change < args.tol_pct and converged_at is None:
            converged_at = tess
            flag = "  <- convergence tolerance first met here"
        print(f"Tess_U={tess:4d}  n_stations={nx:3d}  "
              f"pct_change={pct_change:.3f}%{flag}")

    if converged_at is not None:
        print(f"\nRECOMMENDATION: Tess_U >= {converged_at} looks sufficient "
              f"(consecutive-level change stayed under {args.tol_pct}% from "
              f"there on). Consider adding a safety margin (e.g. 1.5x) "
              f"before locking this in for production runs.")
    else:
        print(f"\nDid not reach {args.tol_pct}% convergence within the tested "
              f"range - try increasing --tess-max in the sweep script, or "
              f"loosen --tol-pct if this precision is sufficient for your use.")

    # ---- plot, if matplotlib is available ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        tess_vals = [r[0] for r in results]
        metric_vals = [r[2] for r in results]

        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(tess_vals, metric_vals, marker="o")
        ax.set_xlabel("Tess_U (spanwise tessellation)")
        ax.set_ylabel(args.metric)
        ax.set_title(f"OpenVSP DegenGeom convergence: {args.metric}")
        ax.grid(True, alpha=0.3)
        if converged_at is not None:
            ax.axvline(converged_at, color="r", linestyle="--",
                        label=f"converged @ Tess_U={converged_at}")
            ax.legend()
        fig.tight_layout()
        out_png = os.path.join(os.path.dirname(args.manifest) or ".",
                                "convergence_plot.png")
        fig.savefig(out_png, dpi=150)
        print(f"\nWrote plot: {out_png}")
    except ImportError:
        print("\n(matplotlib not available - skipping plot; the printed "
              "table above has everything you need)")


if __name__ == "__main__":
    main()
