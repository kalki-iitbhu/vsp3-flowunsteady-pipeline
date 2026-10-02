#!/usr/bin/env python3
"""
openvsp_prop_to_vtk.py
-----------------------
Extract the mean-camber-surface geometry of a propeller (or any lifting
surface) from an OpenVSP DegenGeom CSV export, and write it out as VTK
mesh(es) - one per blade/copy found in the file - ready to load in ParaView
or compare against another tool's geometry (e.g. FLOWUnsteady).

REQUIRED INPUT FILE
--------------------
An OpenVSP "DegenGeom" CSV export of your propeller model:
    1. Build/open your propeller as a PROP component in OpenVSP.
    2. Analysis -> DegenGeom
    3. Set the output file type to CSV, give it a filename, and run it.
This script only needs that ..._DegenGeom.csv file (the companion .m file,
if OpenVSP also writes one, is not required here).

WHAT GETS EXTRACTED
--------------------
DegenGeom reduces the full 3D surface through several representations
(full surface, plate, stick, point). The "PLATE" section is the mean-camber
surface: at every span/chord station it stores both:
  - x, y, z                     : the flattened plate reference point
  - xxCamber, xyCamber, xzCamber: the ACTUAL 3D camber-surface location
                                   (the true midpoint between the upper and
                                   lower surface at that station)
This script extracts the xxCamber/xyCamber/xzCamber columns, since those are
the literal "mid of upper and lower surface" points.

A propeller with multiple blades will contain one PLATE section per blade
copy (SurfNdx 0, 1, 2, ...) - this script writes one VTK file per blade,
plus a combined file with all blades together, plus a flat CSV of all
extracted points.

OUTPUT
------
For an input like "myprop_DegenGeom.csv", this writes into --outdir:
    camber_surface_blade0.vtk, camber_surface_blade1.vtk, ...   (per blade)
    camber_surface_all_blades.vtk                                (combined)
    camber_surface_points.csv                                    (flat table)

USAGE
-----
    python3 openvsp_prop_to_vtk.py myprop_DegenGeom.csv --outdir ./out
    python3 openvsp_prop_to_vtk.py myprop_DegenGeom.csv --blade 0
"""

import argparse
import csv
import os


def parse_plate_sections(path):
    """Read all PLATE sections (one per blade/surface copy) from a
    DegenGeom CSV file. Returns a list of dicts, each with:
        nxsecs, npnts, points (list of dict keyed by column name)
    """
    with open(path, "r", newline="") as f:
        lines = f.readlines()

    blades = []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i].strip()
        if line.startswith("PLATE,"):
            parts = line.split(",")
            nxsecs, npnts = int(parts[1]), int(parts[2])

            i += 1
            if not lines[i].strip().startswith("# nx,ny,nz"):
                raise ValueError(f"Unexpected DegenGeom format near line {i}")
            i += 1
            for _ in range(nxsecs):
                i += 1  # per-cross-section normal vectors, not needed here

            header_line = lines[i].strip()
            if not header_line.startswith("# x,y,z"):
                raise ValueError(f"Expected plate-points header near line {i}, "
                                  f"got: {header_line}")
            cols = [c.strip() for c in header_line.lstrip("#").split(",")]
            i += 1

            pts = []
            for _ in range(nxsecs * npnts):
                vals = [float(x) for x in lines[i].split(",")]
                pts.append(dict(zip(cols, vals)))
                i += 1

            blades.append({"nxsecs": nxsecs, "npnts": npnts, "points": pts})
            continue
        i += 1

    if not blades:
        raise ValueError(
            "No PLATE sections found. Make sure the DegenGeom export "
            "actually ran on a surface (wing/prop) component, not just "
            "a point-mass or duct-only geometry."
        )
    return blades


def write_structured_grid_vtk(path, nxsecs, npnts, points):
    """One blade's camber surface as a VTK legacy STRUCTURED_GRID.
    Point ordering: chordwise fastest-varying, spanwise slowest-varying."""
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("OpenVSP DegenGeom camber surface (mid of upper/lower)\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {npnts} {nxsecs} 1\n")
        f.write(f"POINTS {nxsecs * npnts} float\n")
        for p in points:
            f.write(f"{p['xxCamber']:.9e} {p['xyCamber']:.9e} {p['xzCamber']:.9e}\n")


def write_combined_polydata_vtk(path, blades):
    """All blades together as independent quad patches in one POLYDATA file."""
    all_points, quads, offset = [], [], 0
    for b in blades:
        nx, npn = b["nxsecs"], b["npnts"]
        for p in b["points"]:
            all_points.append((p["xxCamber"], p["xyCamber"], p["xzCamber"]))
        for j in range(nx - 1):
            for k in range(npn - 1):
                p0 = offset + j * npn + k
                p1 = offset + j * npn + (k + 1)
                p2 = offset + (j + 1) * npn + (k + 1)
                p3 = offset + (j + 1) * npn + k
                quads.append((p0, p1, p2, p3))
        offset += nx * npn

    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("OpenVSP DegenGeom camber surface - all blades\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write(f"POINTS {len(all_points)} float\n")
        for x, y, z in all_points:
            f.write(f"{x:.9e} {y:.9e} {z:.9e}\n")
        f.write(f"POLYGONS {len(quads)} {len(quads) * 5}\n")
        for q in quads:
            f.write(f"4 {q[0]} {q[1]} {q[2]} {q[3]}\n")


def write_csv(path, blades):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["blade", "xsec_idx", "chord_idx", "x", "y", "z"])
        for bi, b in enumerate(blades):
            nx, npn = b["nxsecs"], b["npnts"]
            for j in range(nx):
                for k in range(npn):
                    p = b["points"][j * npn + k]
                    w.writerow([bi, j, k, p["xxCamber"], p["xyCamber"], p["xzCamber"]])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("degengeom_csv", help="Path to the OpenVSP DegenGeom CSV export")
    ap.add_argument("--outdir", default=".", help="Directory to write output files into")
    ap.add_argument("--blade", type=int, default=None,
                     help="Only export this blade index (default: export all found)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    blades = parse_plate_sections(args.degengeom_csv)
    print(f"Found {len(blades)} PLATE section(s) (blade/surface copies) "
          f"in {args.degengeom_csv}")

    selected = blades if args.blade is None else [blades[args.blade]]
    indices = range(len(blades)) if args.blade is None else [args.blade]

    for bi, b in zip(indices, selected):
        out_path = os.path.join(args.outdir, f"camber_surface_blade{bi}.vtk")
        write_structured_grid_vtk(out_path, b["nxsecs"], b["npnts"], b["points"])
        print(f"  wrote {out_path}  ({b['nxsecs']} span stations x "
              f"{b['npnts']} chordwise points)")

    if args.blade is None and len(blades) > 1:
        combined_path = os.path.join(args.outdir, "camber_surface_all_blades.vtk")
        write_combined_polydata_vtk(combined_path, blades)
        print(f"  wrote {combined_path}")

    csv_path = os.path.join(args.outdir, "camber_surface_points.csv")
    write_csv(csv_path, blades if args.blade is None else selected)
    print(f"  wrote {csv_path}")


if __name__ == "__main__":
    main()
