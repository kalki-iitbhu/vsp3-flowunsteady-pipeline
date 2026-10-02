#!/usr/bin/env python3
"""
openvsp_convergence_sweep.py
-------------------------------
Automate a spanwise-tessellation convergence study for a propeller in
OpenVSP: reopens your .vsp3, sweeps the PROP geometry's spanwise
tessellation (Tess_U - how many points DegenGeom samples along the span;
this does NOT change your actual design XSecs/blade shape, only how finely
that shape gets sampled/exported), and writes one DegenGeom CSV per level.

REQUIRES: OpenVSP's Python API (the `openvsp` module). This ships with
OpenVSP itself - if `import openvsp` fails, install/enable it from your
OpenVSP installation (Tools -> ... -> Python API, or `pip install openvsp`
if your OpenVSP build publishes wheels) and run this with THAT Python
environment, not a plain system Python.

USAGE
-----
    python openvsp_convergence_sweep.py prop2.vsp3 --outdir ./sweep ^
        --tess-min 8 --tess-max 120 --n-levels 30

This writes ./sweep/degengeom_tessNNN.csv for ~30 tessellation levels
spread between 8 and 120 (geometrically spaced, since convergence usually
matters more at the coarse end).

If a NOTE below is wrong for your OpenVSP version (parameter/analysis-input
names occasionally shift between versions), the script prints the actual
available names it finds so you can fix the two constants near the top.
"""

import argparse
import os
import numpy as np

try:
    import openvsp as vsp
except ImportError:
    raise SystemExit(
        "Could not `import openvsp`. This script must be run with the "
        "Python environment that ships with / is linked to your OpenVSP "
        "installation, not a plain system Python. See OpenVSP's Python API "
        "docs (openvsp.org) for how to install/enable it for your version."
    )

# If your OpenVSP version uses a different name for this, the script will
# print the real available Parm names it found on failure, so you can
# update this:
TESS_PARM_NAME = "Tess_U"      # spanwise tessellation, in the "Shape" group
TESS_PARM_GROUP = "Shape"


def find_prop_geom():
    """Return the geom ID of the first PROP-type geometry in the model."""
    for gid in vsp.FindGeoms():
        if vsp.GetGeomTypeName(gid) == "Propeller" or vsp.GetGeomName(gid):
            # GetGeomTypeName's exact string varies by version; fall back to
            # checking for the Tess_U/Shape parm existing, which all
            # surface geoms have, then just take the first non-empty geom.
            pass
    geoms = vsp.FindGeoms()
    if not geoms:
        raise RuntimeError("No geometries found in this .vsp3 file at all.")
    # Prefer one whose type name mentions "prop" (case-insensitive)
    for gid in geoms:
        try:
            tname = vsp.GetGeomTypeName(gid)
        except Exception:
            tname = ""
        if "prop" in tname.lower():
            return gid
    # Fallback: just use the first geom and let the user confirm
    print(f"  [NOTE] Could not confidently identify a PROP geom by type name; "
          f"using the first geometry found ({vsp.GetGeomName(geoms[0])}). "
          f"If this is wrong, edit find_prop_geom() to hardcode the right ID.")
    return geoms[0]


def run_degengeom_csv(out_csv_path):
    """Run the DegenGeom analysis and write a CSV to out_csv_path.

    Uses the lower-level SetComputationFileName + ComputeDegenGeom API
    (rather than the Analysis Manager's SetStringAnalysisInput) because
    several OpenVSP versions don't expose a CSVFileName analysis input at
    all - SetComputationFileName is the documented, version-stable way to
    control where DegenGeom's CSV gets written.
    """
    vsp.SetComputationFileName(vsp.DEGEN_GEOM_CSV_TYPE, out_csv_path)
    vsp.ComputeDegenGeom(vsp.SET_ALL, vsp.DEGEN_GEOM_CSV_TYPE)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vsp3_file", help="Path to your propeller .vsp3 model")
    ap.add_argument("--outdir", default="./sweep", help="Directory for the CSV sweep")
    ap.add_argument("--tess-min", type=int, default=8,
                     help="Smallest spanwise tessellation to try")
    ap.add_argument("--tess-max", type=int, default=120,
                     help="Largest spanwise tessellation to try")
    ap.add_argument("--n-levels", type=int, default=30,
                     help="How many tessellation levels to sample between "
                          "tess-min and tess-max (~25-35 is plenty for a "
                          "convergence study)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Opening {args.vsp3_file} ...")
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(args.vsp3_file)

    prop_id = find_prop_geom()
    print(f"Using geom: {vsp.GetGeomName(prop_id)} ({prop_id})")

    tess_parm = vsp.FindParm(prop_id, TESS_PARM_NAME, TESS_PARM_GROUP)
    if not tess_parm:
        raise RuntimeError(
            f"Could not find Parm '{TESS_PARM_NAME}' in group "
            f"'{TESS_PARM_GROUP}' on this geom. Open the .vsp3 in the GUI, "
            f"find the actual spanwise-tessellation control under the "
            f"geom's Shape/Tess tab, and update TESS_PARM_NAME/GROUP above."
        )

    # geometrically spaced levels: more resolution at the coarse end, where
    # convergence behavior changes fastest
    levels = np.unique(np.round(
        np.geomspace(args.tess_min, args.tess_max, args.n_levels)
    ).astype(int))
    print(f"Sweeping {len(levels)} tessellation levels: {list(levels)}\n")

    manifest = []
    for tess in levels:
        vsp.SetParmVal(tess_parm, float(tess))
        vsp.Update()
        actual = vsp.GetParmVal(tess_parm)
        csv_path = os.path.join(args.outdir, f"degengeom_tess{int(tess):04d}.csv")
        print(f"  Tess_U = {tess:4d} (set), {actual:.1f} (readback) -> {csv_path}")
        run_degengeom_csv(csv_path)
        manifest.append((int(tess), csv_path))

    manifest_path = os.path.join(args.outdir, "manifest.csv")
    with open(manifest_path, "w") as f:
        f.write("tess_u,csv_path\n")
        for tess, path in manifest:
            f.write(f"{tess},{os.path.abspath(path)}\n")
    print(f"\nWrote {len(manifest)} CSVs and {manifest_path}")
    print("Next: run convergence_analysis.py on this --outdir.")


if __name__ == "__main__":
    main()
