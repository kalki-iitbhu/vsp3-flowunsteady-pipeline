#!/usr/bin/env python3
"""
bem_vsp3_to_flowunsteady.py
=============================================================================
Merged, fully-automated OpenVSP -> XFOIL -> FLOWUnsteady pipeline.

Replaces the old two-step "run auto_generate_polars.py, then manually copy
chorddist/twistdist/etc. out of the .bem into the Julia template" workflow
(see old README Step 2 / Step 4). Now you give it BOTH source files exported
from the SAME OpenVSP propeller and it does everything:

    1. .bem file   (File -> Export -> Blade Element)
       -> real, exported chorddist / twistdist / sweepdist(Skew) /
          heightdist(Rake) distribution tables (r/R = 0 hub row auto-added)

    2. .vsp3 file  (the project file itself)
       -> per-station Four-Series airfoil params (Camber/CamberLoc/ThickChord)
       -> exact NACA contour generation
       -> local XFOIL run per station -> REAL polar CSVs
          (with the narrowing-AoA fallback chain + nearest-station polar
          substitution for stations that never converge)

OUTPUT (all under -o/--outdir, default "./flowunsteady_output"):
    airfoils/contours/*.dat        per-vsp3-station NACA coordinates
    airfoils/*_polar.csv           per-vsp3-station XFOIL polar
                                    (Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr --
                                    matches AirfoilPrep.read_polar's expected
                                    7-column layout)
    <prefix>_bladetable.csv        full .bem blade table (reference only)
    <prefix>_chorddist.csv / _pitchdist.csv / _heightdist.csv / _skewdist.csv
    <prefix>.jl                    <-- the final, ready-to-run FLOWUnsteady script
                                        (chorddist/twistdist/sweepdist/heightdist
                                        come from the .bem; airfoil_contours come
                                        from the vsp3+XFOIL run, with REAL polar
                                        filenames -- no "polar_placeholder.csv")

USAGE:
    python bem_vsp3_to_flowunsteady.py my_prop.bem my_prop.vsp3 [-o OUTDIR] [-p PREFIX]

REQUIREMENTS:
    - Python 3 + numpy
    - xfoil.exe -> set XFOIL_PATH below
    - Same limitation as before: only Four-Series XSecs in the .vsp3 are
      auto-handled for the airfoil/polar side.
"""

import argparse
import csv
import os
import re
import subprocess
import sys

import numpy as np

# ============================== CONFIG ======================================

XFOIL_PATH = r"D:\Cranefield\XFOIL6.99\xfoil.exe"  # <-- EDIT THIS to your local xfoil.exe

# Rotor operating conditions -- used only to estimate local Reynolds number
# per vsp3 station (does not affect geometry or the .bem-derived distributions).
RPM = 9200
J   = 0.4
RHO = 1.225
MU  = 1.81e-5

NCRIT      = 9
ALPHA_STEP = 0.5
ITER       = 100
N_PANELS   = 160
XFOIL_TIMEOUT_S = 90

ALPHA_SWEEP_CHAIN = [
    (-8, 12),   # full attempt
    (-6, 8),    # first fallback
    (-4, 6),    # second fallback, for genuinely stubborn thin/high-Re sections
]

XS_FOUR_SERIES = 4


# =============================================================================
# ============================  .BEM  PARSING  ===============================
# =============================================================================

def _clean_lines(filepath):
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        raw = f.readlines()
    return [ln.rstrip("\n").rstrip("\r") for ln in raw]


def parse_bem(filepath):
    """Parse a .bem file -> (header dict, table_header, table_rows, sections, section_order)."""
    lines = _clean_lines(filepath)
    n = len(lines)
    i = 0

    header = {}
    while i < n:
        line = lines[i].strip()
        if line == "":
            i += 1
            continue
        if line.startswith("Section") or line.lower().startswith("radius/r"):
            break
        if ":" in line:
            key, val = line.split(":", 1)
            header[key.strip()] = val.strip()
        i += 1

    while i < n and lines[i].strip() == "":
        i += 1

    table_header = None
    table_rows = []
    if i < n and lines[i].strip().lower().startswith("radius/r"):
        table_header = [c.strip() for c in lines[i].split(",")]
        i += 1
        while i < n and lines[i].strip() != "" and not lines[i].strip().startswith("Section"):
            parts = [p.strip() for p in lines[i].split(",") if p.strip() != ""]
            if parts:
                table_rows.append([float(v) for v in parts])
            i += 1

    while i < n and lines[i].strip() == "":
        i += 1

    sections = {}
    section_order = []
    while i < n:
        line = lines[i].strip()
        if line.startswith("Section"):
            m = re.match(r"Section\s+(\d+)", line)
            sec_num = int(m.group(1)) if m else len(sections)
            i += 1
            pts = []
            while i < n and lines[i].strip() != "" and not lines[i].strip().startswith("Section"):
                parts = [p.strip() for p in lines[i].split(",") if p.strip() != ""]
                if len(parts) >= 2:
                    pts.append((float(parts[0]), float(parts[1])))
                i += 1
            sections[sec_num] = pts
            section_order.append(sec_num)
        else:
            i += 1

    return header, table_header, table_rows, sections, section_order


def _find_col(table_header, keyword):
    keyword = keyword.lower()
    for idx, name in enumerate(table_header):
        if name.lower().startswith(keyword):
            return idx
    return None


def add_hub_at_zero(rows, r_col=0):
    """Duplicate the hub (first) row's values at r/R=0 and prepend it, unless
    the table already starts at r/R=0."""
    if not rows:
        return rows
    rows_sorted = sorted(rows, key=lambda r: r[r_col])
    hub_row = rows_sorted[0]
    if abs(hub_row[r_col]) < 1e-9:
        return rows_sorted
    zero_row = list(hub_row)
    zero_row[r_col] = 0.0
    return [zero_row] + rows_sorted


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in rows:
            w.writerow(row)


def _fmt(v):
    return f"{v:.10g}"


def format_julia_matrix(rows, per_line=4, indent="    "):
    entries = [f"{_fmt(a)} {_fmt(b)}" for a, b in rows]
    lines = []
    for i in range(0, len(entries), per_line):
        chunk = entries[i:i + per_line]
        suffix = ";" if i + per_line < len(entries) else ""
        lines.append(indent + "; ".join(chunk) + suffix)
    return "[\n" + "\n".join(lines) + "\n]"


# =============================================================================
# ============================  .VSP3  PARSING  ==============================
# =============================================================================

def parse_vsp3(vsp3_path):
    with open(vsp3_path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    def top_val(name):
        m = re.search(rf"<{name} Value=\"([^\"]+)\"", text)
        return float(m.group(1)) if m else None

    diameter = top_val("Diameter")
    num_blade = top_val("NumBlade")
    if diameter is None:
        raise RuntimeError("Could not find <Diameter> in .vsp3 -- is this a PROP geom?")

    R_tip = diameter / 2.0
    B = int(num_blade) if num_blade else None

    chord_block_m = re.search(r'<ParmContainer[^>]*Name="Chord"[^>]*>(.*?)</ParmContainer>', text, re.S)
    chorddist = []
    if chord_block_m:
        cblock = chord_block_m.group(1)
        r_vals = {int(i): float(v) for i, v in re.findall(r'<r_(\d+) Value="([^"]+)"', cblock)}
        c_vals = {int(i): float(v) for i, v in re.findall(r'<crd_(\d+) Value="([^"]+)"', cblock)}
        for idx in sorted(r_vals):
            if idx in c_vals:
                chorddist.append((r_vals[idx], c_vals[idx]))
    if not chorddist:
        print("  [WARNING] Could not parse Chord curve control points -- "
              "falling back to a flat chord/R=0.15 assumption for Re estimates only.")
        chorddist = [(0.0, 0.15), (1.0, 0.15)]

    xsecsurf_match = re.search(r"<XSecSurf>(.*?)</XSecSurf>", text, re.S)
    if not xsecsurf_match:
        raise RuntimeError("No <XSecSurf> block found -- is this a PROP geom .vsp3?")
    body = xsecsurf_match.group(1)

    station_blocks = re.findall(r"<XSec>\s*<ParmContainer>.*?</XSec>\s*</XSec>", body, re.S)
    if not station_blocks:
        raise RuntimeError("No XSec station blocks found inside XSecSurf")

    stations = []
    skipped = []
    for blk in station_blocks:
        def getval(name):
            m = re.search(rf"<{name} Value=\"([^\"]+)\"", blk)
            return float(m.group(1)) if m else None

        type_m = re.search(r"<Type>(\d+)</Type>", blk)
        xsec_type = int(type_m.group(1)) if type_m else None

        r_frac = getval("RadiusFrac")
        camber = getval("Camber")
        camber_loc = getval("CamberLoc")
        thick_chord = getval("ThickChord")

        if xsec_type == XS_FOUR_SERIES and None not in (r_frac, camber, camber_loc, thick_chord):
            stations.append({"r_R": r_frac, "camber": camber,
                              "camber_loc": camber_loc, "thickness": thick_chord})
        else:
            print(f"  [skip] station r/R={r_frac}: unsupported XSec type ({xsec_type}) "
                  f"-- only Four-Series is auto-handled currently")
            skipped.append((r_frac, xsec_type))

    stations.sort(key=lambda s: s["r_R"])
    return {"R_tip": R_tip, "B": B, "chorddist": chorddist, "stations": stations,
            "skipped_stations": skipped}


# --------------------------- pre-flight sanity checks ------------------------

# Thresholds are heuristics, not hard physical limits -- they flag the
# combinations most likely to (a) make XFOIL hang/fail to converge, or
# (b) indicate a genuinely degenerate/mistaken station definition.
THICK_ROOT_WARN   = 0.20   # t/c above this is a common non-convergence cause
THIN_TIP_WARN     = 0.04   # t/c below this, combined with camber, is the classic hang case
HIGH_CAMBER_WARN  = 0.03   # camber (as a fraction) considered "high" for the thin-tip check
DTHICK_PER_SPAN_WARN = 0.5 # |delta t/c| / |delta r/R| between adjacent stations
DCAMBER_PER_SPAN_WARN = 0.3


def preflight_check(vsp3_parsed):
    """Scan every parsed vsp3 station BEFORE any XFOIL run and print warnings
    for geometry likely to be slow, non-convergent, or degenerate. Purely
    informational -- never blocks the run, just tells you what to expect and
    gives you a chance to Ctrl+C and fix the .vsp3 first if you want to."""
    stations = vsp3_parsed["stations"]
    skipped = vsp3_parsed.get("skipped_stations", [])
    warnings = []

    if not stations:
        warnings.append("No Four-Series stations found at all -- nothing to run XFOIL on. "
                         "Check that every XSec in the PROP geom uses the Four-Series airfoil type.")

    if skipped:
        detail = ", ".join(f"r/R={r} (type {t})" for r, t in skipped)
        warnings.append(f"{len(skipped)} station(s) skipped because they aren't Four-Series "
                         f"and will have NO geometry/polar/entry in airfoil_contours: {detail}")

    for s in stations:
        r_R, camber, camber_loc, t = s["r_R"], s["camber"], s["camber_loc"], s["thickness"]
        label = f"r/R={r_R:.4f}"

        if camber != 0 and camber_loc == 0:
            warnings.append(f"{label}: Camber={camber:.4f} but CamberLoc=0 -- degenerate "
                             f"camber line (will produce a flat/garbage camber). Set "
                             f"CamberLoc or zero out Camber for this station.")

        if t > THICK_ROOT_WARN:
            warnings.append(f"{label}: thickness={t*100:.1f}% t/c is quite thick -- a "
                             f"common cause of XFOIL non-convergence, especially near the root.")

        if t < THIN_TIP_WARN and camber > HIGH_CAMBER_WARN:
            warnings.append(f"{label}: thin ({t*100:.1f}% t/c) + high camber "
                             f"({camber*100:.1f}%) -- the classic XFOIL hang/blow-up "
                             f"combo, especially if this is near the tip.")

    for a, b in zip(stations, stations[1:]):
        dr = b["r_R"] - a["r_R"]
        if dr <= 0:
            continue
        dt = abs(b["thickness"] - a["thickness"])
        dc = abs(b["camber"] - a["camber"])
        if dt / dr > DTHICK_PER_SPAN_WARN:
            warnings.append(f"Large thickness jump between r/R={a['r_R']:.4f} and "
                             f"r/R={b['r_R']:.4f} ({dt*100:.1f}% t/c over span {dr:.3f}) "
                             f"-- abrupt spanwise change, may cause convergence/aero jumps.")
        if dc / dr > DCAMBER_PER_SPAN_WARN:
            warnings.append(f"Large camber jump between r/R={a['r_R']:.4f} and "
                             f"r/R={b['r_R']:.4f} ({dc*100:.1f}% camber over span {dr:.3f}).")

    print("=" * 60)
    if warnings:
        print(f"PRE-FLIGHT WARNINGS ({len(warnings)}) -- not fatal, flagged before "
              f"burning any XFOIL time:")
        for w in warnings:
            print(f"  [!] {w}")
    else:
        print("Pre-flight check: no obvious geometry red flags found.")
    print("=" * 60)


def check_bem_vsp3_consistency(bem_header, vsp3_parsed):
    """Cross-check tip radius and blade count between the .bem and .vsp3.
    Informational only -- the .bem value always wins downstream since it's
    what generate_rotor() actually uses."""
    R_bem = float(bem_header.get("Diameter", 0.0)) / 2.0
    B_bem = int(float(bem_header.get("Num_Blade", 0))) if bem_header.get("Num_Blade") else None
    R_vsp3 = vsp3_parsed["R_tip"]
    B_vsp3 = vsp3_parsed["B"]

    if abs(R_vsp3 - R_bem) > 0.02 * max(R_bem, 1e-9):
        print(f"  [NOTE] .bem tip radius ({R_bem:.4f} m) and .vsp3 tip radius "
              f"({R_vsp3:.4f} m) differ by more than 2% -- double check both files "
              f"came from the same, up-to-date geometry. The .bem value will be used.")
    if B_vsp3 is not None and B_bem is not None and B_vsp3 != B_bem:
        print(f"  [NOTE] .bem Num_Blade ({B_bem}) and .vsp3 NumBlade ({B_vsp3}) "
              f"disagree -- using the .bem value.")


# --------------------------- NACA geometry / XFOIL --------------------------

def naca4_exact(camber, camber_loc, thickness, n=100):
    m, p, t = camber, camber_loc, thickness
    beta = np.linspace(0, np.pi, n + 1)
    x = (1 - np.cos(beta)) / 2

    yt = 5 * t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
                  + 0.2843*x**3 - 0.1015*x**4)

    if m == 0 or p == 0:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    else:
        yc = np.where(x < p, m/p**2 * (2*p*x - x**2),
                       m/(1-p)**2 * ((1-2*p) + 2*p*x - x**2))
        dyc = np.where(x < p, 2*m/p**2 * (p - x),
                        2*m/(1-p)**2 * (p - x))

    theta = np.arctan(dyc)
    xu = x - yt*np.sin(theta); yu = yc + yt*np.cos(theta)
    xl = x + yt*np.sin(theta); yl = yc - yt*np.cos(theta)

    X = np.concatenate([xu[::-1], xl[1:]])
    Y = np.concatenate([yu[::-1], yl[1:]])
    X[0] = 1.0
    X[-1] = 1.0
    return X, Y


def station_label(s):
    return f"r{s['r_R']:.4f}".replace(".", "")


def write_dat(path, X, Y, label="AIRFOIL"):
    with open(path, "w") as f:
        f.write(f"{label}\n")
        for x, y in zip(X, Y):
            f.write(f"{x:.6f} {y:.6f}\n")


def local_chord(r_R, chorddist):
    rs = [p[0] for p in chorddist]
    cs = [p[1] for p in chorddist]
    return float(np.interp(r_R, rs, cs))


def estimate_reynolds(r_R, R_tip, chord_over_R):
    n_rev_s = RPM / 60.0
    D = 2 * R_tip
    v_axial = J * n_rev_s * D
    v_tang = 2*np.pi*n_rev_s * r_R * R_tip
    v_total = np.hypot(v_axial, v_tang)
    chord = chord_over_R * R_tip
    Re = RHO * v_total * chord / MU
    return max(Re, 1000.0)


def run_xfoil(dat_path, re_number, raw_polar_path, alpha_min, alpha_max):
    if not os.path.isfile(XFOIL_PATH):
        raise FileNotFoundError(
            f"xfoil.exe not found at {XFOIL_PATH} -- download it from "
            f"https://web.mit.edu/drela/Public/web/xfoil/ and update XFOIL_PATH"
        )
    if os.path.exists(raw_polar_path):
        os.remove(raw_polar_path)

    commands = f"""PLOP
G

LOAD {dat_path}

PANE
PPAR
N {N_PANELS}


OPER
VISC {re_number:.0f}
ITER {ITER}
PACC
{raw_polar_path}

ASEQ {alpha_min} {alpha_max} {ALPHA_STEP}

PACC

QUIT
"""
    return subprocess.run([XFOIL_PATH], input=commands, text=True,
                           capture_output=True, timeout=XFOIL_TIMEOUT_S)


def parse_xfoil_polar(raw_polar_path):
    """Parse XFOIL's raw PACC dump. Keeps the full 7-column layout
    (Alpha, Cl, Cd, Cdp, Cm, Top_Xtr, Bot_Xtr) since AirfoilPrep.read_polar
    (Julia side) expects that exact column count/order -- dropping Cdp or
    Top_Xtr/Bot_Xtr shifts every later column and corrupts Cm."""
    rows = []
    if not os.path.isfile(raw_polar_path):
        return rows
    with open(raw_polar_path, "r", errors="replace") as f:
        lines = f.readlines()

    data_start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("alpha") and "CL" in line:
            data_start = i + 2
            break
    if data_start is None:
        return rows

    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            vals = [float(x) for x in parts[:7]]
        except ValueError:
            continue
        while len(vals) < 7:      # Top_Xtr/Bot_Xtr can be absent at some Re/Ncrit combos
            vals.append(0.0)
        alpha, cl, cd, cdp, cm, top_xtr, bot_xtr = vals[:7]
        rows.append((alpha, cl, cd, cdp, cm, top_xtr, bot_xtr))
    return rows


def run_xfoil_chain(dat_path, re_number, raw_polar_path, label):
    for i, (amin, amax) in enumerate(ALPHA_SWEEP_CHAIN):
        tag = "full sweep" if i == 0 else f"fallback {i} ({amin} to {amax} deg)"
        try:
            run_xfoil(dat_path, re_number, raw_polar_path, amin, amax)
            rows = parse_xfoil_polar(raw_polar_path)
            if len(rows) >= 8:
                if i > 0:
                    print(f"  [{tag}] succeeded: {len(rows)} points")
                return rows
            print(f"  [{tag}] only {len(rows)} points -- trying next fallback...")
        except subprocess.TimeoutExpired:
            print(f"  [{tag}] timed out -- trying next fallback...")
    print(f"  [FAILED] {label}: no attempt in the fallback chain produced a "
          f"usable polar -- this station needs manual attention")
    return []


def write_polar_csv(path, rows):
    """Write the 7-column format AirfoilPrep.read_polar expects:
    Alpha, Cl, Cd, Cdp, Cm, Top_Xtr, Bot_Xtr (same layout as XFOIL's own
    PACC dump). A 4-column Alpha/Cl/Cd/Cm-only file misaligns every column
    from Cdp onward on the Julia side and corrupts Cm with `missing`."""
    with open(path, "w") as f:
        f.write("Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr\n")
        for alpha, cl, cd, cdp, cm, top_xtr, bot_xtr in rows:
            f.write(f"{alpha:.2f},{cl:.5f},{cd:.5f},{cdp:.5f},{cm:.5f},"
                    f"{top_xtr:.4f},{bot_xtr:.4f}\n")


def generate_airfoil_contours(vsp3_parsed, outdir):
    """Run the NACA-generation + XFOIL chain for every Four-Series vsp3 station.
    Returns a list of (r/R, contour_relpath, polar_filename) tuples, with the
    root duplicated at r/R=0.0 (required by FLOWVLM), and with failed stations'
    POLAR (never geometry) substituted from the nearest converged station."""
    R_tip = vsp3_parsed["R_tip"]
    chorddist = vsp3_parsed["chorddist"]
    stations = vsp3_parsed["stations"]

    airfoils_dir = os.path.join(outdir, "airfoils")
    contour_dir = os.path.join(airfoils_dir, "contours")
    os.makedirs(contour_dir, exist_ok=True)

    entries = []          # (r_R, contour_relpath, polar_name_or_None)
    failed_stations = []

    for s in stations:
        label = station_label(s)
        print(f"--- Station r/R={s['r_R']:.4f}  (camber={s['camber']:.4f}, "
              f"camber_loc={s['camber_loc']:.3f}, t/c={s['thickness']:.4f}) ---")

        X, Y = naca4_exact(s["camber"], s["camber_loc"], s["thickness"])
        dat_path = os.path.join(contour_dir, f"{label}.dat")
        write_dat(dat_path, X, Y, label=label)
        print(f"  Wrote contour: {dat_path}  ({len(X)} points)")

        chord_over_R = local_chord(s["r_R"], chorddist)
        Re = estimate_reynolds(s["r_R"], R_tip, chord_over_R)
        print(f"  Estimated Re = {Re:.0f}")

        raw_polar_path = os.path.join(airfoils_dir, f"{label}_raw.txt")
        csv_path = os.path.join(airfoils_dir, f"{label}_polar.csv")

        rows = run_xfoil_chain(dat_path, Re, raw_polar_path, label)
        if rows:
            write_polar_csv(csv_path, rows)
            print(f"  -> {len(rows)} AoA points written to {csv_path}\n")
        else:
            csv_path = None
            failed_stations.append(s)
            print()

        contour_rel = os.path.relpath(dat_path, outdir).replace("\\", "/")
        polar_name = os.path.basename(csv_path) if csv_path else None
        entries.append((s["r_R"], contour_rel, polar_name))

    # --- substitute failed stations' polar with the nearest converged one ---
    ok_entries = [(r, c, p) for (r, c, p) in entries if p is not None]
    if failed_stations:
        print("=" * 60)
        print(f"{len(failed_stations)} station(s) never converged. Substituting "
              f"the nearest successfully-converged station's POLAR ONLY "
              f"(contour geometry stays station-specific and accurate):")
        resolved = []
        for r, c, p in entries:
            if p is not None:
                resolved.append((r, c, p))
                continue
            if not ok_entries:
                raise RuntimeError("No stations converged at all -- cannot substitute. "
                                    "Check XFOIL_PATH and try again.")
            nearest = min(ok_entries, key=lambda t: abs(t[0] - r))
            print(f"  r/R={r:.4f}: using polar from r/R={nearest[0]:.4f} ({nearest[2]})")
            resolved.append((r, c, nearest[2]))
        entries = resolved

    # --- duplicate root at r/R=0.0, required by FLOWVLM ---
    if entries:
        r0, c0, p0 = entries[0]
        entries = [(0.0, c0, p0)] + entries

    return entries


# =============================================================================
# =========================  JULIA SCRIPT ASSEMBLY  ==========================
# =============================================================================

def build_julia_script(bem_header, extended_rows, r_col, chord_col, twist_col,
                        rake_col, skew_col, airfoil_entries, prefix):
    diameter = float(bem_header.get("Diameter", 0.0))
    R = diameter / 2.0
    num_blades = int(float(bem_header.get("Num_Blade", 2)))

    hub_rR = extended_rows[1][r_col] if len(extended_rows) > 1 and extended_rows[0][r_col] == 0.0 \
        else extended_rows[0][r_col]
    Rhub = hub_rR * R

    chorddist_rows = [[row[r_col], row[chord_col]] for row in extended_rows]
    twistdist_rows = [[row[r_col], row[twist_col]] for row in extended_rows]
    sweepdist_rows = [[row[r_col], row[skew_col]] for row in extended_rows]   # Skew/R -> sweepdist
    heightdist_rows = [[row[r_col], row[rake_col]] for row in extended_rows]  # Rake/R -> heightdist

    chorddist_str = format_julia_matrix(chorddist_rows)
    twistdist_str = format_julia_matrix(twistdist_rows)
    sweepdist_str = format_julia_matrix(sweepdist_rows)
    heightdist_str = format_julia_matrix(heightdist_rows)

    contour_entries = []
    for r, c, p in airfoil_entries:
        contour_entries.append(
            f'    ({_fmt(r)}, readdlm(joinpath(@__DIR__, "{c}"), skipstart=1), "{p}"),'
        )
    airfoil_contours_str = "[\n" + "\n".join(contour_entries) + "\n]"

    run_name = f"{prefix}-propeller"

    script = f'''#=##############################################################################
# DESCRIPTION
    Simulation of Custom {num_blades}-Bladed Propeller
    Auto-generated by bem_vsp3_to_flowunsteady.py from {prefix}.bem + {prefix}.vsp3.
    chorddist/twistdist/sweepdist/heightdist come from the REAL .bem export.
    airfoil_contours come from exact NACA geometry (per-station Camber/CamberLoc/
    ThickChord parsed from the .vsp3) with REAL locally-computed XFOIL polars --
    no manual copy-paste, no placeholder polar files.
=###############################################################################

import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm
import DelimitedFiles: readdlm

run_name        = "{run_name}"        # Name of this simulation
save_path       = run_name                  # Where to save this simulation
paraview        = true                      # Whether to visualize with Paraview

# ----------------- SIMULATION PARAMETERS --------------------------------------
RPM             = {RPM}                     # RPM
J               = {J}                       # Advance ratio Vinf/(nD)
AOA             = 0                         # (deg) Angle of attack
rho             = {RHO}                     # (kg/m^3) air density
mu              = {MU}                      # (kg/ms) air dynamic viscosity
speedofsound    = 342.35                    # (m/s) speed of sound

# ----------------- GEOMETRY PARAMETERS ----------------------------------------
pitch           = 0.0                       # (deg) collective pitch of blades
CW              = false                     # Clock-wise rotation
xfoil           = false                     # Real polars already generated below -- don't re-run XFOIL internally
ncrit           = {NCRIT}                   # Turbulence criterion for XFOIL
n               = 20                        # Number of blade elements per blade
r               = 1/5                       # Geometric expansion of elements

# Main rotor parameters (from {prefix}.bem)
R               = {_fmt(R)}                     # (m) Radius of blade tip
Rhub            = {_fmt(Rhub)}                     # (m) Radius of hub
B               = {num_blades}                         # Number of blades

# --- Geometric Distributions (r/R = 0 hub row duplicated, from the .bem) ------
chorddist = {chorddist_str}

twistdist = {twistdist_str}

sweepdist = {sweepdist_str}

heightdist = {heightdist_str}

# ----------------- AIRFOIL CONTOURS + REAL XFOIL POLARS -----------------------
# (r/R, contour matrix loaded from airfoils/contours/*.dat, polar CSV filename)
# generated from the .vsp3's per-station Four-Series parameters -- see
# airfoils/ next to this script.
airfoil_contours = {airfoil_contours_str}

# ----------------- 1) VEHICLE DEFINITION --------------------------------------
println("Generating geometry...")

# Polar filenames in airfoil_contours (e.g. "r02000_polar.csv") are looked up
# by AirfoilPrep.read_polar as joinpath(data_path, "airfoils", filename) --
# it appends "airfoils" itself, so data_path must be THIS script's directory,
# NOT .../airfoils (that would double up to .../airfoils/airfoils/...).
data_path = @__DIR__

rotor = uns.generate_rotor(R, Rhub, B, chorddist, twistdist, sweepdist, heightdist, airfoil_contours;
                                        pitch=pitch,
                                        n=n, CW=CW, blade_r=r,
                                        altReD=[RPM, J, mu/rho],
                                        xfoil=xfoil,
                                        ncrit=ncrit,
                                        data_path=data_path,
                                        verbose=true,
                                        verbose_xfoil=false,
                                        plot_disc=true
                                        );

println("Generating vehicle...")

# --- SOLVER CONFIGURATION ---
VehicleType     = uns.UVLMVehicle           # Unsteady solver
const_solution  = VehicleType==uns.QVLMVehicle
nrevs           = 4
nsteps_per_rev  = 36
nsteps          = const_solution ? 2 : nrevs*nsteps_per_rev
ttot            = nsteps/nsteps_per_rev / (RPM/60)
p_per_step      = 2
shed_starting   = true
shed_unsteady   = true
max_particles   = ((2*n+1)*B)*nsteps*p_per_step + 1
sigma_rotor_surf= R/40
lambda_vpm      = 2.125
sigma_vpm_overwrite = lambda_vpm * 2*pi*R/(nsteps_per_rev*p_per_step)
vlm_rlx         = 0.7
hubtiploss_correction = vlm.hubtiploss_nocorrection
vpm_viscous     = vpm.Inviscid()

if VehicleType == uns.QVLMVehicle
    uns.vlm.VLMSolver._mute_warning(true)
end

# Generate vehicle systems
system = vlm.WingSystem()
vlm.addwing(system, "Rotor", rotor)
rotors = [rotor];
rotor_systems = (rotors, );

wake_system = vlm.WingSystem()
if VehicleType != uns.QVLMVehicle
    vlm.addwing(wake_system, "Rotor", rotor)
end

vehicle = VehicleType(   system;
                            rotor_systems=rotor_systems,
                            wake_system=wake_system
                         );

# ------------- 2) MANEUVER DEFINITION -----------------------------------------
Vvehicle(t) = zeros(3)
anglevehicle(t) = zeros(3)
RPMcontrol(t) = 1.0

angles = ()
RPMs = (RPMcontrol, )

maneuver = uns.KinematicManeuver(angles, RPMs, Vvehicle, anglevehicle)

# ------------- 3) SIMULATION DEFINITION ---------------------------------------
magVinf         = J*RPM/60*(2*R)
Vinf(X, t)      = magVinf*[cosd(AOA), sind(AOA), 0]

Vref = 0.0
RPMref = RPM
Vinit = Vref*Vvehicle(0)
Winit = pi/180*(anglevehicle(1e-6) - anglevehicle(0))/(1e-6*ttot)

simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                    Vinit=Vinit, Winit=Winit);

# ------------- 4) MONITORS DEFINITIONS ----------------------------------------
figs, figaxs = [], []

monitor_rotor = uns.generate_monitor_rotors(rotors, J, rho, RPM, nsteps;
                                            t_scale=RPM/60,
                                            t_lbl="Revolutions",
                                            out_figs=figs,
                                            out_figaxs=figaxs,
                                            save_path=save_path,
                                            run_name=run_name,
                                            figname="rotor monitor",
                                            )

# ------------- 5) RUN SIMULATION ----------------------------------------------
println("Running simulation...")

uns.run_simulation(simulation, nsteps;
                    Vinf=Vinf,
                    rho=rho, mu=mu, sound_spd=speedofsound,
                    p_per_step=p_per_step,
                    max_particles=max_particles,
                    vpm_viscous=vpm_viscous,
                    sigma_vlm_surf=sigma_rotor_surf,
                    sigma_rotor_surf=sigma_rotor_surf,
                    sigma_vpm_overwrite=sigma_vpm_overwrite,
                    vlm_rlx=vlm_rlx,
                    hubtiploss_correction=hubtiploss_correction,
                    shed_unsteady=shed_unsteady,
                    shed_starting=shed_starting,
                    extra_runtime_function=monitor_rotor,
                    save_path=save_path,
                    run_name=run_name,
                    );

# ----------------- 6) VISUALIZATION -------------------------------------------
if paraview
    println("Calling Paraview...")
    files = joinpath(save_path, run_name*"_pfield...xmf;")
    for bi in 1:B
        global files
        files *= run_name*"_Rotor_Blade$(bi)_loft...vtk;"
        files *= run_name*"_Rotor_Blade$(bi)_vlm...vtk;"
    end
    run(`paraview --data=$(files)`)
end
'''
    return script


# =============================================================================
# =================================  MAIN  ====================================
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Merge a .bem + .vsp3 export of the same OpenVSP propeller "
                    "into one self-contained FLOWUnsteady .jl script with real XFOIL polars.")
    ap.add_argument("bem_file", help="Path to the .bem file")
    ap.add_argument("vsp3_file", help="Path to the .vsp3 file (same propeller)")
    ap.add_argument("-o", "--outdir", default="flowunsteady_output", help="Output directory")
    ap.add_argument("-p", "--prefix", default=None,
                     help="Filename prefix for outputs (default: .bem filename stem)")
    args = ap.parse_args()

    if not os.path.isfile(args.bem_file):
        sys.exit(f"Error: file not found: {args.bem_file}")
    if not os.path.isfile(args.vsp3_file):
        sys.exit(f"Error: file not found: {args.vsp3_file}")

    prefix = args.prefix or os.path.splitext(os.path.basename(args.bem_file))[0]
    os.makedirs(args.outdir, exist_ok=True)

    # ---------------- .bem: distributions ----------------
    print("=" * 60)
    print(f"Parsing {args.bem_file} ...")
    bem_header, table_header, table_rows, sections, section_order = parse_bem(args.bem_file)
    print("=== .bem header ===")
    for k, v in bem_header.items():
        print(f"  {k}: {v}")

    if not (table_header and table_rows):
        sys.exit("Error: no blade geometry table (Radius/R, Chord/R, Twist, Rake/R, "
                  "Skew/R, ...) found in this .bem file -- cannot build chorddist/"
                  "twistdist/sweepdist/heightdist.")

    r_col = _find_col(table_header, "radius/r") or 0
    chord_col = _find_col(table_header, "chord/r")
    twist_col = _find_col(table_header, "twist")
    rake_col = _find_col(table_header, "rake/r")
    skew_col = _find_col(table_header, "skew/r")
    missing = [name for name, col in
               [("Chord/R", chord_col), ("Twist", twist_col),
                ("Rake/R", rake_col), ("Skew/R", skew_col)] if col is None]
    if missing:
        sys.exit(f"Error: missing column(s) {missing} in the .bem blade table.")

    extended_rows = add_hub_at_zero(table_rows, r_col=r_col)

    full_path = os.path.join(args.outdir, f"{prefix}_bladetable.csv")
    write_csv(full_path, table_header, extended_rows)
    print(f"Wrote {full_path}  ({len(extended_rows)} rows)")

    dist_specs = [
        ("chorddist", "chord/r", "c/R"),
        ("pitchdist", "twist", "twist (deg)"),
        ("heightdist", "rake/r", "z/R (Rake, out-of-plane)"),
        ("skewdist", "skew/r", "y/R (Skew, in-plane)"),
    ]
    for fname_suffix, keyword, out_label in dist_specs:
        col = _find_col(table_header, keyword)
        dist_rows = [[row[r_col], row[col]] for row in extended_rows]
        dist_path = os.path.join(args.outdir, f"{prefix}_{fname_suffix}.csv")
        write_csv(dist_path, ["r/R", out_label], dist_rows)
        print(f"Wrote {dist_path}  ({len(dist_rows)} rows)")

    # ---------------- .vsp3: airfoil contours + real XFOIL polars ----------------
    print("=" * 60)
    print(f"Parsing {args.vsp3_file} ...")
    vsp3_parsed = parse_vsp3(args.vsp3_file)
    print(f"  Rtip = {vsp3_parsed['R_tip']:.4f} m   NumBlade = {vsp3_parsed['B']}")
    print(f"Found {len(vsp3_parsed['stations'])} Four-Series stations\n")

    # ---------------- pre-flight checks (before any XFOIL time is spent) ------
    check_bem_vsp3_consistency(bem_header, vsp3_parsed)
    preflight_check(vsp3_parsed)

    airfoil_entries = generate_airfoil_contours(vsp3_parsed, args.outdir)

    # ---------------- assemble final Julia script ----------------
    print("=" * 60)
    julia_code = build_julia_script(
        bem_header, extended_rows, r_col, chord_col, twist_col, rake_col, skew_col,
        airfoil_entries, prefix,
    )
    jl_path = os.path.join(args.outdir, f"{prefix}.jl")
    with open(jl_path, "w") as f:
        f.write(julia_code)

    print(f"Done. Self-contained FLOWUnsteady script -> {jl_path}")
    print(f"Make sure the 'airfoils/' folder next to it (under {args.outdir}) travels with it.")
    print("=" * 60)


if __name__ == "__main__":
    main()
