#!/usr/bin/env python3
"""
First production run of the UMM public SPARC residual-table pipeline.

Follows notes/UMM_SPARC_Residual_Table_Design.md §3–§9 exactly.
Does not modify frozen paper/UMM_Paper_Draft_v13.tex.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

# Allow running as script: python analysis/sparc_residuals/run_pipeline.py
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from analysis.sparc_residuals.fit import (
    epsilon_for_models,
    fit_global_A,
    fit_nfw_m200,
    mean_outer_sigma_bar,
    model_v_baryon,
    model_v_umm,
)
from analysis.sparc_residuals.parse import build_galaxies, parse_mass_models, parse_table1
from analysis.sparc_residuals.physics import (
    N_G_PRIMARY,
    UPSILON_BULGE_FIDUCIAL,
    UPSILON_DISK_FIDUCIAL,
    v_total_nfw_baryon,
)
from analysis.sparc_residuals.select import (
    apply_mandatory_cuts,
    cut_tallies,
    selected_ids,
)

SPARC_DIR = _ROOT / "Data" / "SPARC"
OUT_DIR = _ROOT / "analysis" / "sparc_residuals" / "results"
INPUT_FILES = {
    "table1": SPARC_DIR / "SPARC_Lelli2016c.mrt",
    "mass_models": SPARC_DIR / "MassModels_Lelli2016c.mrt",
    "rotmod_zip": SPARC_DIR / "Rotmod_LTG.zip",
    "sfb_zip": SPARC_DIR / "sfb_LTG.zip",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def selection_id_hash(galaxy_list: List[str], cut_spec: Dict[str, Any], file_hashes: Dict[str, str]) -> str:
    payload = json.dumps(
        {"galaxies": galaxy_list, "cuts": cut_spec, "files": file_hashes},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def stratum_from_sigma(values: List[float], val: float) -> str:
    """Tertiles of mean outer Σ_bar over the sample."""
    if not values or math.isnan(val):
        return "unknown"
    s = sorted(values)
    n = len(s)
    t1 = s[max(0, n // 3 - 1)]
    t2 = s[max(0, (2 * n) // 3 - 1)]
    if val <= t1:
        return "low"
    if val <= t2:
        return "intermediate"
    return "high"


def run(upsilon_disk: float = UPSILON_DISK_FIDUCIAL, tag: str = "fiducial") -> Dict[str, Any]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- §2 hashes ---
    file_hashes = {}
    file_paths = {}
    for key, path in INPUT_FILES.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing SPARC input: {path}")
        file_hashes[key] = sha256_file(path)
        file_paths[key] = str(path.relative_to(_ROOT))

    # --- parse ---
    table1 = parse_table1(INPUT_FILES["table1"])
    mass_models = parse_mass_models(INPUT_FILES["mass_models"])
    galaxies = build_galaxies(table1, mass_models)

    # --- §3 selection ---
    cut_spec = {
        "Q": 1,
        "incl_min_deg": 30.0,
        "R_last_over_Rdisc_min": 2.5,
        "N_R_min": 8,
        "upsilon_disk": upsilon_disk,
        "upsilon_bulge": UPSILON_BULGE_FIDUCIAL,
        "n_g": N_G_PRIMARY,
    }
    cuts = apply_mandatory_cuts(galaxies, r_last_over_rdisc_min=2.5)
    tallies = cut_tallies(cuts)
    selected = selected_ids(cuts)

    # Widen coverage if N < 20 (design §3.3)
    coverage_widened = False
    if len(selected) < 20:
        cuts = apply_mandatory_cuts(galaxies, r_last_over_rdisc_min=2.0)
        tallies = cut_tallies(cuts)
        selected = selected_ids(cuts)
        coverage_widened = True
        cut_spec["R_last_over_Rdisc_min"] = 2.0
        cut_spec["coverage_widened"] = True

    if len(selected) < 5:
        raise RuntimeError(f"Too few galaxies after cuts: N={len(selected)}")

    sel_id = selection_id_hash(selected, cut_spec, file_hashes)

    # --- §5 global A ---
    A_best, A_err, chi2_min = fit_global_A(
        galaxies, selected, upsilon_disk=upsilon_disk, n_g=N_G_PRIMARY
    )

    # optional packing-index sensitivity (one-line)
    A_pack, _, _ = fit_global_A(
        galaxies, selected, upsilon_disk=upsilon_disk, n_g=0.45
    )

    # --- §6–§7 residuals ---
    rows: List[Dict[str, Any]] = []
    sigma_outers: List[float] = []
    for gid in selected:
        g = galaxies[gid]
        sig_o = mean_outer_sigma_bar(g, upsilon_disk=upsilon_disk)
        sigma_outers.append(sig_o)

    for gid in selected:
        g = galaxies[gid]
        meta = g.meta
        assert meta is not None

        v_umm = model_v_umm(g, A_best, upsilon_disk=upsilon_disk, n_g=N_G_PRIMARY)
        v_bar_only = model_v_baryon(g, upsilon_disk=upsilon_disk)
        m200 = fit_nfw_m200(g, upsilon_disk=upsilon_disk)
        v_nfw = [
            v_total_nfw_baryon(
                p.r_kpc,
                p.vgas,
                p.vdisk,
                p.vbul,
                m200,
                upsilon_disk,
                UPSILON_BULGE_FIDUCIAL,
            )
            for p in g.points
        ]

        eps_u, n_out, r_out = epsilon_for_models(g, v_umm)
        eps_b, _, _ = epsilon_for_models(g, v_bar_only)
        eps_n, _, _ = epsilon_for_models(g, v_nfw)
        sig_o = mean_outer_sigma_bar(g, upsilon_disk=upsilon_disk)

        rows.append(
            {
                "galaxy_id": gid,
                "morph_type": meta.morph_type,
                "Q": meta.Q,
                "incl_deg": meta.incl_deg,
                "dist_Mpc": meta.D_Mpc,
                "R_disc_kpc": meta.Rdisk_kpc,
                "R_last_kpc": g.r_last,
                "R_out_kpc": r_out,
                "N_outer": n_out,
                "mean_Sigma_bar_outer": sig_o,
                "stratum": "",  # filled below
                "epsilon_umm": eps_u,
                "epsilon_baryon": eps_b,
                "epsilon_nfw": eps_n,
                "A_global": A_best,
                "Upsilon_star": upsilon_disk,
                "gas_source": "SPARC_HI",
                "notes": f"M200_NFW={m200:.3e};coverage={cuts[gid].coverage_ratio:.2f}",
            }
        )

    sigs = [r["mean_Sigma_bar_outer"] for r in rows if not math.isnan(r["mean_Sigma_bar_outer"])]
    for r in rows:
        r["stratum"] = stratum_from_sigma(sigs, r["mean_Sigma_bar_outer"])

    # --- summaries ---
    def med(key: str) -> float:
        vals = [r[key] for r in rows if r[key] == r[key]]
        return statistics.median(vals) if vals else float("nan")

    eps_u_list = [r["epsilon_umm"] for r in rows if r["epsilon_umm"] == r["epsilon_umm"]]
    frac_lt = sum(1 for e in eps_u_list if e < 0.25) / len(eps_u_list) if eps_u_list else float("nan")

    summary = {
        "N_gal": len(rows),
        "med_epsilon_umm": med("epsilon_umm"),
        "med_epsilon_baryon": med("epsilon_baryon"),
        "med_epsilon_nfw": med("epsilon_nfw"),
        "frac_eps_lt_0.25": frac_lt,
        "A_global": A_best,
        "A_err": A_err,
        "A_unit": "m s^-2",
        "chi2_min_outer": chi2_min,
        "selection_id": sel_id,
        "upsilon_disk": upsilon_disk,
        "upsilon_bulge": UPSILON_BULGE_FIDUCIAL,
        "n_g": N_G_PRIMARY,
        "A_global_ng045_sensitivity": A_pack,
        "coverage_widened_to_2.0": coverage_widened,
        "cut_tallies": tallies,
        "tag": tag,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }

    # --- write outputs (fiducial only overwrites primary deliverable names) ---
    suffix = "" if tag == "fiducial" else f"_{tag}"

    manifest = {
        "design_note": "notes/UMM_SPARC_Residual_Table_Design.md",
        "mapping": {
            "Table1.mrt": file_paths["table1"],
            "Table2/mass_models": file_paths["mass_models"],
            "Rotmod_LTG.zip": file_paths["rotmod_zip"],
            "sfb_LTG.zip": file_paths["sfb_zip"],
        },
        "input_files": {
            k: {"path": file_paths[k], "sha256": file_hashes[k]} for k in file_hashes
        },
        "cut_spec": cut_spec,
        "cut_tallies": tallies,
        "frozen_galaxy_list": selected,
        "N_selected": len(selected),
        "selection_id": sel_id,
        "coverage_widened_to_2.0": coverage_widened,
        "notes": (
            "MassModels_Lelli2016c.mrt is the radial mass-model table (design Table2). "
            "SPARC_Lelli2016c.mrt is the galaxy sample (design Table1). "
            "Gas surface density uses SPARC Vgas thin-sheet proxy; gas_source=SPARC_HI. "
            "No post-hoc residual-based galaxy drops."
        ),
    }
    man_path = OUT_DIR / f"selection_manifest{suffix}.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")

    # CSV
    cols = [
        "galaxy_id",
        "morph_type",
        "Q",
        "incl_deg",
        "dist_Mpc",
        "R_disc_kpc",
        "R_last_kpc",
        "R_out_kpc",
        "N_outer",
        "mean_Sigma_bar_outer",
        "stratum",
        "epsilon_umm",
        "epsilon_baryon",
        "epsilon_nfw",
        "A_global",
        "Upsilon_star",
        "gas_source",
        "notes",
    ]
    csv_path = OUT_DIR / f"residuals_main{suffix}.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in cols})

    json_path = OUT_DIR / f"residuals_main{suffix}.json"
    json_path.write_text(json.dumps(rows, indent=2) + "\n")

    sum_path = OUT_DIR / f"residuals_summary{suffix}.json"
    sum_path.write_text(json.dumps(summary, indent=2) + "\n")

    return {
        "manifest_path": str(man_path),
        "summary": summary,
        "rows": rows,
        "selected": selected,
        "tallies": tallies,
        "A_best": A_best,
        "A_err": A_err,
    }


def write_report(fid: Dict[str, Any], ups07: Dict[str, Any] | None) -> None:
    s = fid["summary"]
    lines = [
        "# UMM SPARC Residual-Table Production Report",
        "",
        f"**Timestamp (UTC):** {s['timestamp_utc']}",
        f"**Selection ID:** `{s['selection_id']}`",
        f"**Design contract:** `notes/UMM_SPARC_Residual_Table_Design.md` §3–§9",
        "",
        "## Status of numbers",
        "",
        "These residuals are **measured catalog residuals** from a documented run of the",
        "public SPARC mass models (Lelli, McGaugh & Schombert 2016) with pre-registered",
        "selection cuts. They **replace** the illustrative schematic sample-summaries of",
        "frozen `paper/UMM_Paper_Draft_v13.tex` for use in a future versioned successor",
        "(v14+). The frozen v13 body was **not** edited in this run.",
        "",
        "## Sample",
        "",
        f"- **N_gal:** {s['N_gal']}",
        f"- **Coverage cut:** R_last ≥ {2.0 if s.get('coverage_widened_to_2.0') else 2.5} R_disc"
        + (" (widened from 2.5 because N would have been <20)" if s.get("coverage_widened_to_2.0") else ""),
        f"- **Cut tallies:** `{json.dumps(s['cut_tallies'])}`",
        f"- **Optional robustness cut (§3.2):** not applied as a rejection cut (no galaxies dropped by ΔV quality filter).",
        "",
        "## Fitted global amplitude",
        "",
        f"- **A = {s['A_global']:.4e} ± {s['A_err']:.4e} m s^-2** (n_g = {s['n_g']}, Υ_* = {s['upsilon_disk']})",
        f"- Outer-disk χ²_min = {s['chi2_min_outer']:.2f}",
        f"- Sensitivity (n_g = 0.45 packing index): A ≈ {s['A_global_ng045_sensitivity']:.4e} m s^-2",
        "",
        "## Median outer-disk residuals",
        "",
        f"| Model | med(ε) |",
        f"|-------|--------|",
        f"| UMM (single global A) | **{s['med_epsilon_umm']:.4f}** |",
        f"| Pure baryons (A=0) | **{s['med_epsilon_baryon']:.4f}** |",
        f"| Simple NFW (1 mass/galaxy) | **{s['med_epsilon_nfw']:.4f}** |",
        "",
        f"- Fraction with ε_umm < 0.25: **{s['frac_eps_lt_0.25']:.3f}**",
        "",
    ]
    if ups07 is not None:
        s7 = ups07["summary"]
        de = s7["med_epsilon_umm"] - s["med_epsilon_umm"]
        lines += [
            "## Optional robustness: Υ_* = 0.7",
            "",
            f"- A(Υ=0.7) = {s7['A_global']:.4e} ± {s7['A_err']:.4e} m s^-2",
            f"- med(ε_umm, Υ=0.7) = {s7['med_epsilon_umm']:.4f}",
            f"- **Δε = med(ε_umm,0.7) − med(ε_umm,0.5) = {de:+.4f}**",
            "",
        ]
    lines += [
        "## Deliverables",
        "",
        "- `analysis/sparc_residuals/results/selection_manifest.json`",
        "- `analysis/sparc_residuals/results/residuals_main.csv`",
        "- `analysis/sparc_residuals/results/residuals_main.json`",
        "- `analysis/sparc_residuals/results/residuals_summary.json`",
        "- `analysis/sparc_residuals/run_pipeline.py` (reproducible entrypoint)",
        "",
        "## Skeptic notes",
        "",
        "- Selection frozen before residual ranking; no high-ε post-hoc drops.",
        "- UMM column uses one global A for all galaxies.",
        "- NFW baseline is one-parameter (M200) with Dutton–Macciò c(M), not full ΛCDM.",
        "- Σ_gas uses Vgas thin-sheet proxy (SPARC_HI); not THINGS map re-reduction.",
        "",
    ]
    (OUT_DIR / "production_report.md").write_text("\n".join(lines) + "\n")


def write_v14_status() -> None:
    text = """# Status for v14 — honesty sentences that may now be updated

After this production run, a future versioned successor to frozen v13 may replace the SPARC honesty language in §Galactic specialization (`\\label{sec:SPARC}`) and the Limitations paragraph “SPARC residual numbers are sample-summaries, not a catalog re-pipeline” with text that cites `analysis/sparc_residuals/results/` (selection_id, measured med(ε), fitted A±σ_A) as a documented public-catalog residual table. Do **not** edit `paper/UMM_Paper_Draft_v13.tex` in place; ship changes only in v14+. The collaboration note and Appendix A residual-geometry dictionary remain untouched. Illustrative wording in figure captions for fig2/table residual panels should be revised only when those figures are regenerated from the measured CSV.
"""
    (OUT_DIR / "status_for_v14.md").write_text(text)


def main() -> int:
    print("Running fiducial pipeline (Υ_*=0.5)...")
    fid = run(upsilon_disk=0.5, tag="fiducial")
    print(
        f"N={fid['summary']['N_gal']}  A={fid['A_best']:.4e}±{fid['A_err']:.4e}  "
        f"med_eps_umm={fid['summary']['med_epsilon_umm']:.4f}"
    )
    print("Running Υ_*=0.7 robustness...")
    ups07 = run(upsilon_disk=0.7, tag="Upsilon0p7")
    write_report(fid, ups07)
    write_v14_status()
    # also copy primary names are already without suffix for fiducial
    print("Wrote results to", OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
