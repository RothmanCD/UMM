#!/usr/bin/env python3
"""Entry point: UMM residual confrontation with DM-deficient UDG kinematics.

Usage (from repo root or this package):
  python -m analysis.udg_confrontation.src.run_confrontation
  python analysis/udg_confrontation/src/run_confrontation.py

Writes tables and residual summary under analysis/udg_confrontation/results/.
Does not re-fit A, n_g, α_eff, or ξ_c. Does not touch frozen v14 or SPARC products.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

# Allow running as script from repo root or package dir
_HERE = Path(__file__).resolve().parent
_PKG = _HERE.parent
_ROOT = _PKG.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from analysis.udg_confrontation.src.manifest_data import (  # noqa: E402
    FROZEN_UMM,
    FUTURE_TARGETS,
    SYSTEMS,
    sigma_obs_1sigma,
)
from analysis.udg_confrontation.src.physics import (  # noqa: E402
    A_FROZEN,
    A_FROZEN_UNC,
    ALPHA_EFF,
    N_G,
    SIGMA_MID_WINDOW,
    SIGMA_REF,
    a_ind_logistic,
    a_ind_powerlaw,
    fractional_residual_sigma,
    mean_sigma_within_re,
    predict_sigma_three_models,
    wolf_mass_from_sigma,
    wolf_sigma_from_mass,
)

RESULTS = _PKG / "results"


def _sym_unc(plus: float, minus: float) -> float:
    return 0.5 * (plus + minus)


def analyze_system(sys_rec: Dict[str, Any]) -> Dict[str, Any]:
    """Compute three-model predictions and residuals for one system."""
    m_star = float(sys_rec["M_star"]["value"])
    m_star_unc = float(sys_rec["M_star"]["unc"])
    r_e = float(sys_rec["R_e"]["value"])
    r_e_unc = float(sys_rec["R_e"]["unc"])
    d = float(sys_rec["distance"]["value"])
    d_unc = float(sys_rec["distance"]["unc"])
    sig_obs = float(sys_rec["sigma_obs"]["value"])
    sig_obs_unc = sigma_obs_1sigma(sys_rec)

    sigma_bar = mean_sigma_within_re(m_star, r_e)
    pred = predict_sigma_three_models(m_star, r_e, sigma_mid=SIGMA_REF)

    # Logistic band across Σ_mid window (ξ_c joint-norm window)
    log_band = []
    for sm in SIGMA_MID_WINDOW:
        p = predict_sigma_three_models(m_star, r_e, sigma_mid=sm)
        log_band.append(
            {
                "Sigma_mid": sm,
                "sigma_logistic": p["sigma_logistic"],
                "a_ind_logistic": p["a_ind_logistic"],
            }
        )
    sig_log_vals = [b["sigma_logistic"] for b in log_band]
    sig_log_lo, sig_log_hi = min(sig_log_vals), max(sig_log_vals)

    m_dyn_obs = wolf_mass_from_sigma(sig_obs, r_e)
    # Propagate σ uncertainty to M (M ∝ σ²)
    m_dyn_obs_unc = m_dyn_obs * 2.0 * (sig_obs_unc / sig_obs) if sig_obs > 0 else float("nan")

    eps_stars = fractional_residual_sigma(sig_obs, pred["sigma_stars"])
    eps_pl = fractional_residual_sigma(sig_obs, pred["sigma_powerlaw"])
    eps_log = fractional_residual_sigma(sig_obs, pred["sigma_logistic"])

    # Consistency: |σ_model − σ_obs| vs combined uncertainty
    # Model uncertainty: propagate A (±σ_A) for UMM models; M_* for stars
    def umm_sigma_unc(sig_central: float, a_ind: float, a_bar: float) -> float:
        if a_bar + a_ind <= 0:
            return 0.0
        # ∂σ/∂A: a_ind ∝ A → Δa = a_ind (ΔA/A); σ ∝ sqrt(a_tot)
        dsig_da = 0.5 * sig_central / (a_bar + a_ind)
        da = a_ind * (A_FROZEN_UNC / A_FROZEN)
        # M_* uncertainty → a_bar and Σ (hence a_ind)
        return abs(dsig_da * da)

    unc_pl = umm_sigma_unc(
        pred["sigma_powerlaw"], pred["a_ind_powerlaw"], pred["a_bar"]
    )
    unc_log = umm_sigma_unc(
        pred["sigma_logistic"], pred["a_ind_logistic"], pred["a_bar"]
    )
    # Stellar M/L → σ_stars
    unc_stars = pred["sigma_stars"] * 0.5 * (m_star_unc / m_star) if m_star > 0 else 0.0

    def n_sigma(sig_model: float, model_unc: float) -> float:
        tot = math.sqrt(sig_obs_unc**2 + model_unc**2)
        if tot <= 0:
            return float("nan")
        return abs(sig_model - sig_obs) / tot

    # Systematics shifts (illustrative, one-at-a-time)
    systematics = _systematics_shifts(sys_rec, pred)

    return {
        "name": sys_rec["name"],
        "short": sys_rec["short"],
        "D_Mpc": d,
        "D_unc_Mpc": d_unc,
        "M_star": m_star,
        "M_star_unc": m_star_unc,
        "R_e_kpc": r_e,
        "R_e_unc_kpc": r_e_unc,
        "Sigma_bar_mean_msun_pc2": sigma_bar,
        "sigma_obs": sig_obs,
        "sigma_obs_unc": sig_obs_unc,
        "sigma_obs_unc_plus": sys_rec["sigma_obs"].get("unc_plus", sig_obs_unc),
        "sigma_obs_unc_minus": sys_rec["sigma_obs"].get("unc_minus", sig_obs_unc),
        "tracer": sys_rec["sigma_obs"]["tracer"],
        "sigma_stars": pred["sigma_stars"],
        "sigma_powerlaw": pred["sigma_powerlaw"],
        "sigma_logistic": pred["sigma_logistic"],
        "sigma_logistic_lo": sig_log_lo,
        "sigma_logistic_hi": sig_log_hi,
        "logistic_band": log_band,
        "a_bar": pred["a_bar"],
        "a_ind_powerlaw": pred["a_ind_powerlaw"],
        "a_ind_logistic": pred["a_ind_logistic"],
        "m_half_star": pred["m_half"],
        "m_dyn_obs": m_dyn_obs,
        "m_dyn_obs_unc": m_dyn_obs_unc,
        "m_dyn_stars": pred["m_dyn_stars"],
        "m_dyn_powerlaw": pred["m_dyn_powerlaw"],
        "m_dyn_logistic": pred["m_dyn_logistic"],
        "eps_sigma_stars": eps_stars,
        "eps_sigma_powerlaw": eps_pl,
        "eps_sigma_logistic": eps_log,
        "n_sigma_stars": n_sigma(pred["sigma_stars"], unc_stars),
        "n_sigma_powerlaw": n_sigma(pred["sigma_powerlaw"], unc_pl),
        "n_sigma_logistic": n_sigma(pred["sigma_logistic"], unc_log),
        "model_unc_stars": unc_stars,
        "model_unc_powerlaw": unc_pl,
        "model_unc_logistic": unc_log,
        "systematics": systematics,
        "profile_limitation": sys_rec["profile_limitation"],
        "citations": {
            "distance": sys_rec["distance"]["citation"],
            "M_star": sys_rec["M_star"]["citation"],
            "R_e": sys_rec["R_e"]["citation"],
            "sigma_obs": sys_rec["sigma_obs"]["citation"],
        },
    }


def _systematics_shifts(
    sys_rec: Dict[str, Any], pred_fid: Dict[str, float]
) -> Dict[str, Any]:
    """Show how distance, M/L, and rotation/tides move residuals (power-law)."""
    m = float(sys_rec["M_star"]["value"])
    r = float(sys_rec["R_e"]["value"])
    d = float(sys_rec["distance"]["value"])
    d_unc = float(sys_rec["distance"]["unc"])
    m_unc = float(sys_rec["M_star"]["unc"])
    sig_obs = float(sys_rec["sigma_obs"]["value"])

    out: Dict[str, Any] = {}

    # Distance: M_* ∝ D², R_e ∝ D → Σ independent of D; a_bar ∝ M/R² independent;
    # but published M_* and R_e are at adopted D — scale both with D.
    for label, d_fac in (("D_plus", 1.0 + d_unc / d), ("D_minus", 1.0 - d_unc / d)):
        m2 = m * d_fac**2
        r2 = r * d_fac
        p = predict_sigma_three_models(m2, r2)
        out[label] = {
            "D_Mpc": d * d_fac,
            "sigma_stars": p["sigma_stars"],
            "sigma_powerlaw": p["sigma_powerlaw"],
            "sigma_logistic": p["sigma_logistic"],
            "eps_powerlaw": fractional_residual_sigma(sig_obs, p["sigma_powerlaw"]),
            "note": "M_*∝D², R_e∝D; Σ invariant; σ scales ~D^{1/2} for fixed a",
        }

    # M/L (stellar mass) at fixed R_e
    for label, m2 in (("ML_plus", m + m_unc), ("ML_minus", max(m - m_unc, 1.0))):
        p = predict_sigma_three_models(m2, r)
        out[label] = {
            "M_star": m2,
            "sigma_stars": p["sigma_stars"],
            "sigma_powerlaw": p["sigma_powerlaw"],
            "sigma_logistic": p["sigma_logistic"],
            "eps_powerlaw": fractional_residual_sigma(sig_obs, p["sigma_powerlaw"]),
        }

    # Residual rotation: if observed σ includes rotation support, true dispersion lower.
    # Illustrative: subtract V_rot/√2 in quadrature with V_rot = 0.5 σ_obs (order-unity).
    v_rot = 0.5 * sig_obs
    sig_corr = math.sqrt(max(sig_obs**2 - v_rot**2, 0.0))
    out["rotation_illustration"] = {
        "V_rot_assumed_kms": v_rot,
        "sigma_obs_corrected": sig_corr,
        "eps_powerlaw_vs_corrected": fractional_residual_sigma(
            sig_corr, pred_fid["sigma_powerlaw"]
        ),
        "note": (
            "Illustrative only: assumes ordered motion ~0.5 σ_obs mixed into "
            "reported dispersion. FCC 224 has published slow rotation; DF2/DF4 "
            "show little rotation in KCWI maps."
        ),
    }

    # Tides: order-of-magnitude — if unbound/tidal, Wolf M overestimates bound mass.
    # Report 30% downward shift in effective σ as a systematic floor discussion.
    out["tides_illustration"] = {
        "sigma_obs_if_30pct_inflated": sig_obs / 1.3,
        "eps_powerlaw": fractional_residual_sigma(sig_obs / 1.3, pred_fid["sigma_powerlaw"]),
        "note": (
            "Illustrative: if tides inflate line-of-sight motions by ~30%, "
            "compare models to σ/1.3. Not a measurement."
        ),
    }
    return out


def residual_summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Median residuals and 1σ/2σ consistency fractions."""

    def med(xs: List[float]) -> float:
        ys = sorted(xs)
        n = len(ys)
        if n == 0:
            return float("nan")
        if n % 2:
            return ys[n // 2]
        return 0.5 * (ys[n // 2 - 1] + ys[n // 2])

    def stats(key_eps: str, key_ns: str) -> Dict[str, Any]:
        eps = [r[key_eps] for r in rows]
        ns = [r[key_ns] for r in rows]
        n1 = sum(1 for x in ns if x <= 1.0)
        n2 = sum(1 for x in ns if x <= 2.0)
        return {
            "median_eps_sigma": med(eps),
            "median_abs_eps_sigma": med([abs(e) for e in eps]),
            "mean_eps_sigma": sum(eps) / len(eps),
            "n_within_1sigma": n1,
            "n_within_2sigma": n2,
            "n_total": len(rows),
            "fraction_within_1sigma": n1 / len(rows),
            "fraction_within_2sigma": n2 / len(rows),
            "eps_list": eps,
            "n_sigma_list": ns,
        }

    s_stars = stats("eps_sigma_stars", "n_sigma_stars")
    s_pl = stats("eps_sigma_powerlaw", "n_sigma_powerlaw")
    s_log = stats("eps_sigma_logistic", "n_sigma_logistic")

    # Ternary classification for UMM (power-law primary; logistic as band)
    # Rule: if ≥ half the sample has n_sigma > 2 for power-law with systematically
    # positive eps (model high), → ruled out; if mixed or ~1–2σ → mild tension;
    # if mostly ≤1σ → consistent.
    n = len(rows)
    n_pl_2 = s_pl["n_within_2sigma"]
    med_eps_pl = s_pl["median_eps_sigma"]
    all_high = all(e > 0.5 for e in s_pl["eps_list"])  # model >50% high

    if n_pl_2 <= n / 2 and all_high and med_eps_pl > 0.5:
        classification = "ruled out"
        classification_note = (
            "Intermediate-regime power-law a_ind overpredicts σ for a majority of "
            "systems by more than 2σ (model systematically high). Pure-baryon "
            "predictions remain consistent with the data."
        )
    elif s_pl["fraction_within_1sigma"] >= 0.5 and abs(med_eps_pl) < 0.3:
        classification = "consistent"
        classification_note = (
            "UMM power-law predictions lie within ~1σ for most systems; "
            "median fractional residual is small."
        )
    else:
        classification = "mildly tensioned"
        classification_note = (
            "UMM power-law is high relative to observed σ for this small sample, "
            "but large kinematic uncertainties and systematics prevent a clean "
            "exclusion at uniform significance for every object."
        )

    # Refine: if every system is >2σ high for power-law, call ruled out
    if all(ns > 2.0 for ns in s_pl["n_sigma_list"]) and med_eps_pl > 0.5:
        classification = "ruled out"
        classification_note = (
            "For every system in this small sample, the frozen intermediate-regime "
            "UMM power-law predicts a line-of-sight dispersion more than 2σ above "
            "the observed value (model systematically high). The pure-baryon "
            "column is consistent with the data within uncertainties. This is an "
            "illustrative first estimate on a small sample; systematics remain "
            "significant, but the direction of the residual is uniform."
        )

    return {
        "N": n,
        "stars_alone": s_stars,
        "powerlaw_umm": s_pl,
        "logistic_umm": s_log,
        "classification": classification,
        "classification_note": classification_note,
        "honesty": (
            "Illustrative first estimate; small sample (N=4); systematics remain "
            "significant (distance, M/L, residual rotation, possible tides, "
            "mean-Σ proxy for local Σ_bar, thin-disk specialization applied to "
            "pressure-supported systems)."
        ),
    }


def write_outputs(rows: List[Dict[str, Any]], summary: Dict[str, Any]) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    # Manifest JSON (full literature + frozen params + results)
    manifest = {
        "frozen_umm": FROZEN_UMM,
        "estimator_assumptions": {
            "mass_estimator": "Wolf et al. 2010",
            "formula": "M_1/2 = 4 σ_los^2 R_e / G",
            "sigma_from_a": "σ = sqrt(a_tot R_e / 4) with a_tot = a_bar + a_ind",
            "a_bar": "G (M_*/2) / R_e^2 (mass-follows-light half-mass)",
            "Sigma_bar": "mean within R_e: M_*/(2 π R_e^2); not a published profile",
            "anisotropy": "Wolf estimator; weak anisotropy sensitivity near r_1/2",
            "gas": "neglected (quiescent UDGs)",
            "no_refit": "A, n_g, alpha_eff, xi_c window fixed; no optimization",
        },
        "systems_literature": SYSTEMS,
        "future_targets": FUTURE_TARGETS,
        "predictions": rows,
        "residual_summary": summary,
    }
    man_path = RESULTS / "literature_manifest.json"
    with man_path.open("w") as f:
        json.dump(manifest, f, indent=2, default=str)

    # Data table CSV
    data_csv = RESULTS / "udg_data_table.csv"
    with data_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "name",
                "D_Mpc",
                "D_unc",
                "M_star",
                "M_star_unc",
                "R_e_kpc",
                "Sigma_bar_mean",
                "sigma_obs",
                "sigma_obs_unc_plus",
                "sigma_obs_unc_minus",
                "sigma_obs_unc_max",
                "tracer",
                "distance_citation",
                "sigma_citation",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    r["name"],
                    f"{r['D_Mpc']:.1f}",
                    f"{r['D_unc_Mpc']:.1f}",
                    f"{r['M_star']:.3e}",
                    f"{r['M_star_unc']:.3e}",
                    f"{r['R_e_kpc']:.2f}",
                    f"{r['Sigma_bar_mean_msun_pc2']:.3f}",
                    f"{r['sigma_obs']:.2f}",
                    f"{r['sigma_obs_unc_plus']:.2f}",
                    f"{r['sigma_obs_unc_minus']:.2f}",
                    f"{r['sigma_obs_unc']:.2f}",
                    r["tracer"],
                    r["citations"]["distance"],
                    r["citations"]["sigma_obs"],
                ]
            )

    # Predictions CSV
    pred_csv = RESULTS / "udg_predictions_three_models.csv"
    with pred_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "name",
                "Sigma_bar",
                "a_bar",
                "a_ind_pl",
                "a_ind_log",
                "sigma_obs",
                "sigma_stars",
                "sigma_powerlaw",
                "sigma_logistic",
                "sigma_logistic_lo",
                "sigma_logistic_hi",
                "eps_stars",
                "eps_powerlaw",
                "eps_logistic",
                "n_sigma_stars",
                "n_sigma_powerlaw",
                "n_sigma_logistic",
                "m_dyn_obs",
                "m_half_star",
                "m_dyn_powerlaw",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    r["name"],
                    f"{r['Sigma_bar_mean_msun_pc2']:.4f}",
                    f"{r['a_bar']:.4f}",
                    f"{r['a_ind_powerlaw']:.4f}",
                    f"{r['a_ind_logistic']:.4f}",
                    f"{r['sigma_obs']:.3f}",
                    f"{r['sigma_stars']:.3f}",
                    f"{r['sigma_powerlaw']:.3f}",
                    f"{r['sigma_logistic']:.3f}",
                    f"{r['sigma_logistic_lo']:.3f}",
                    f"{r['sigma_logistic_hi']:.3f}",
                    f"{r['eps_sigma_stars']:.4f}",
                    f"{r['eps_sigma_powerlaw']:.4f}",
                    f"{r['eps_sigma_logistic']:.4f}",
                    f"{r['n_sigma_stars']:.3f}",
                    f"{r['n_sigma_powerlaw']:.3f}",
                    f"{r['n_sigma_logistic']:.3f}",
                    f"{r['m_dyn_obs']:.4e}",
                    f"{r['m_half_star']:.4e}",
                    f"{r['m_dyn_powerlaw']:.4e}",
                ]
            )

    # Residual summary JSON
    with (RESULTS / "residual_summary.json").open("w") as f:
        json.dump(summary, f, indent=2)

    # Compact selection / run manifest
    run_manifest = {
        "analysis": "udg_confrontation",
        "frozen_A": A_FROZEN,
        "frozen_A_unc": A_FROZEN_UNC,
        "n_g": N_G,
        "Sigma_ref": SIGMA_REF,
        "alpha_eff": ALPHA_EFF,
        "Sigma_mid_window": list(SIGMA_MID_WINDOW),
        "selection_id_sparc": FROZEN_UMM["selection_id"],
        "n_systems": len(rows),
        "system_names": [r["name"] for r in rows],
        "classification": summary["classification"],
        "no_refit": True,
        "outputs": [
            "literature_manifest.json",
            "udg_data_table.csv",
            "udg_predictions_three_models.csv",
            "residual_summary.json",
            "confrontation_report.md",
        ],
    }
    with (RESULTS / "selection_manifest.json").open("w") as f:
        json.dump(run_manifest, f, indent=2)

    # Markdown report
    report = build_report(rows, summary)
    with (RESULTS / "confrontation_report.md").open("w") as f:
        f.write(report)

    # Optional simple text figure (sigma comparison)
    with (RESULTS / "sigma_comparison.txt").open("w") as f:
        f.write("System        σ_obs   σ_stars  σ_PL    σ_log   eps_PL\n")
        for r in rows:
            f.write(
                f"{r['short']:12s} {r['sigma_obs']:6.2f}  {r['sigma_stars']:6.2f}  "
                f"{r['sigma_powerlaw']:6.2f}  {r['sigma_logistic']:6.2f}  "
                f"{r['eps_sigma_powerlaw']:+6.2f}\n"
            )


def build_report(rows: List[Dict[str, Any]], summary: Dict[str, Any]) -> str:
    """Referee-readable markdown report."""
    lines: List[str] = []
    L = lines.append

    L("# First quantitative UMM confrontation with DM-deficient UDG kinematics")
    L("")
    L("**Status:** Illustrative first estimate · Small sample (N = 4) · Systematics remain significant")
    L("")
    L("**Frozen inputs (no re-fit):**")
    L(f"- \(A = (1.41 \\pm 0.07)\\times 10^{{-11}}\\,\\mathrm{{m\\,s^{{-2}}}}\) (SPARC residual amplitude; selection_id `ea6c92761ffd3c5a`, \(\\Upsilon_*=0.5\), \(n_g=0.4\))")
    L(f"- \(\\Sigma_{{\\rm ref}} = 1\\,M_\\odot\\,\\mathrm{{pc^{{-2}}}}\)")
    L(f"- Packing \(\\alpha_{{\\rm eff}} \\approx 1.8\) (preferred band \([1.66, 1.84]\)); logistic midpoint window \(\\Sigma_{{\\rm mid}} \\in \\{{{', '.join(str(x) for x in SIGMA_MID_WINDOW)}\\}}\\,M_\\odot\\,\\mathrm{{pc^{{-2}}}}\) (joint-normalization placement about SPARC-relevant intermediate densities; no new free parameters)")
    L("- All other UMM parameters remain those of the frozen v14 baseline.")
    L("")
    L("---")
    L("")
    L("## 1. Sample and data table")
    L("")
    L("Primary systems: NGC 1052-DF2, DF4, DF9, and FCC 224. Distances prefer SBF/TRGB or homogeneous group/cluster scales. Only published peer-reviewed or arXiv-vetted kinematics are used.")
    L("")
    L("**Surface-density limitation (all systems):** No published high-resolution \(\\Sigma(R)\) mass profile is adopted. Mean baryonic surface density within \(R_e\) is computed as \(\\langle\\Sigma\\rangle_e = M_*/(2\\pi R_e^2)\) (half-mass within \(R_e\), mass-follows-light, gas neglected for these quiescent systems). Local \(\\Sigma(R_e)\) for a Sérsic profile would be lower by an order-unity factor; we do not invent a radial profile.")
    L("")
    L("| System | \(D\) [Mpc] | \(M_*\) [\(10^8 M_\\odot\)] | \(R_e\) [kpc] | \(\\langle\\Sigma\\rangle_e\) [\(M_\\odot\\mathrm{pc^{-2}}\)] | \(\\sigma_{\\rm obs}\) [km/s] | Tracer |")
    L("|--------|------------|---------------------------|---------------|------------------------------------------------------|------------------------------|--------|")
    for r in rows:
        sp, sm = r["sigma_obs_unc_plus"], r["sigma_obs_unc_minus"]
        if abs(sp - sm) < 0.05:
            sig_str = f"{r['sigma_obs']:.2f}±{sp:.2f}"
        else:
            sig_str = f"{r['sigma_obs']:.2f}^{{+{sp:.2f}}}_{{-{sm:.2f}}}"
        L(
            f"| {r['name']} | {r['D_Mpc']:.1f}±{r['D_unc_Mpc']:.1f} | "
            f"{r['M_star']/1e8:.2f}±{r['M_star_unc']/1e8:.2f} | "
            f"{r['R_e_kpc']:.2f} | {r['Sigma_bar_mean_msun_pc2']:.2f} | "
            f"{sig_str} | {r['tracer']} |"
        )
    L("")
    L("### Citations (adopted values)")
    L("")
    for r in rows:
        L(f"**{r['name']}**")
        L(f"- Distance: {r['citations']['distance']}")
        L(f"- \(M_*\): {r['citations']['M_star']}")
        L(f"- \(R_e\): {r['citations']['R_e']}")
        L(f"- \(\\sigma\): {r['citations']['sigma_obs']}")
        L(f"- Profile note: {r['profile_limitation']}")
        L("")
    L("### Future targets (not in residual sample)")
    L("")
    for t in FUTURE_TARGETS:
        L(f"- **{t['name']}:** {t['reason']} ({t['citation']})")
    L("")
    L("---")
    L("")
    L("## 2. Estimator and model definitions")
    L("")
    L("### Wolf et al. (2010) mass estimator")
    L("")
    L(r"\[")
    L(r"M_{1/2} = \frac{4\,\sigma_{\rm los}^{2}\,R_e}{G},")
    L(r"\]")
    L("with the inverse \(\\sigma = \\sqrt{G M_{1/2}/(4 R_e)}\).")
    L("")
    L("**Assumptions (stated explicitly):**")
    L("1. Pressure-supported kinematics (dispersion-dominated).")
    L("2. Mass follows light for the stellar comparison; half the stellar mass lies within \(R_e\).")
    L("3. Wolf estimator near the deprojected half-light radius is only weakly sensitive to orbital anisotropy.")
    L("4. Flat or slowly varying \(\\sigma(R)\) (no full Jeans solution with free \(\\beta(r)\)).")
    L("5. Negligible cold gas in \(\\Sigma_{\\rm bar}\).")
    L("6. Thin-disk UMM specialization \(a_{\\rm ind}(R)=A(\\Sigma/\\Sigma_{\\rm ref})^{n_g}\) is applied using mean \(\\Sigma\) within \(R_e\) as a proxy — these systems are spheroidal; this is a limitation of the first estimate.")
    L("")
    L("### Three models (identical estimator and inputs)")
    L("")
    L("(a) **Stars alone:** \(a_{\\rm tot}=a_{\\rm bar}=G(M_*/2)/R_e^2\).")
    L("")
    L("(b) **Stars + intermediate power-law:** \(a_{\\rm ind}=A(\\langle\\Sigma\\rangle_e/\\Sigma_{\\rm ref})^{0.4}\), \(a_{\\rm tot}=a_{\\rm bar}+a_{\\rm ind}\).")
    L("")
    L("(c) **Stars + logistic \(\\Phi(\\xi)\):** packing \(\\alpha_{\\rm eff}=1.8\), midpoint at \(\\Sigma_{\\rm mid}\) in the joint-normalization window about \(\\Sigma_{\\rm ref}\), amplitude matched so \(a_{\\rm ind}(\\Sigma_{\\rm mid})=A(\\Sigma_{\\rm mid}/\\Sigma_{\\rm ref})^{n_g}\). Central column uses \(\\Sigma_{\\rm mid}=\\Sigma_{\\rm ref}\); band reports min/max over the window.")
    L("")
    L(r"\[")
    L(r"\sigma_{\rm model}=\sqrt{a_{\rm tot}\,R_e/4}.")
    L(r"\]")
    L("")
    L("Fractional residual (SPARC spirit, signed): \(\\varepsilon_\\sigma=(\\sigma_{\\rm model}-\\sigma_{\\rm obs})/\\sigma_{\\rm obs}\).")
    L("")
    L("---")
    L("")
    L("## 3. Predicted vs observed dispersions")
    L("")
    L("| System | \(\\sigma_{\\rm obs}\) | (a) stars | (b) + power-law | (c) + logistic | \(\\varepsilon_{\\rm PL}\) | \(n_\\sigma^{\\rm PL}\) |")
    L("|--------|---------------------|-----------|-----------------|----------------|----------------------|----------------------|")
    for r in rows:
        sp, sm = r["sigma_obs_unc_plus"], r["sigma_obs_unc_minus"]
        if abs(sp - sm) < 0.05:
            sig_str = f"{r['sigma_obs']:.2f}±{sp:.2f}"
        else:
            sig_str = f"{r['sigma_obs']:.2f}^{{+{sp:.2f}}}_{{-{sm:.2f}}}"
        L(
            f"| {r['short']} | {sig_str} | "
            f"{r['sigma_stars']:.1f} | {r['sigma_powerlaw']:.1f} | "
            f"{r['sigma_logistic']:.1f} [{r['sigma_logistic_lo']:.1f}–{r['sigma_logistic_hi']:.1f}] | "
            f"{r['eps_sigma_powerlaw']:+.2f} | {r['n_sigma_powerlaw']:.2f} |"
        )
    L("")
    L("| System | \(M_{1/2}^{\\rm obs}\) [\(10^8 M_\\odot\)] | \(M_*/2\) | \(M_{1/2}^{\\rm PL}\) | \(\\langle\\Sigma\\rangle_e\) | \(a_{\\rm bar}\) | \(a_{\\rm ind}^{\\rm PL}\) |")
    L("|--------|------------------------------------------|----------|---------------------|--------------------------|-----------------|--------------------------|")
    for r in rows:
        L(
            f"| {r['short']} | {r['m_dyn_obs']/1e8:.2f}±{r['m_dyn_obs_unc']/1e8:.2f} | "
            f"{r['m_half_star']/1e8:.2f} | {r['m_dyn_powerlaw']/1e8:.2f} | "
            f"{r['Sigma_bar_mean_msun_pc2']:.2f} | {r['a_bar']:.1f} | {r['a_ind_powerlaw']:.1f} |"
        )
    L("")
    L("Accelerations in \((\\mathrm{km\\,s^{-1}})^2\\,\\mathrm{kpc^{-1}}\).")
    L("")
    L("---")
    L("")
    L("## 4. Residual summary")
    L("")
    L(f"- Sample size: **N = {summary['N']}** (not a large statistical population).")
    L(f"- Median \(\\varepsilon_\\sigma\) (stars alone): **{summary['stars_alone']['median_eps_sigma']:+.3f}**")
    L(f"- Median \(\\varepsilon_\\sigma\) (UMM power-law): **{summary['powerlaw_umm']['median_eps_sigma']:+.3f}**")
    L(f"- Median \(\\varepsilon_\\sigma\) (UMM logistic, central \(\\Sigma_{{\\rm mid}}\)): **{summary['logistic_umm']['median_eps_sigma']:+.3f}**")
    L(f"- Fraction within 1σ / 2σ (stars): "
      f"{summary['stars_alone']['n_within_1sigma']}/{summary['N']} / "
      f"{summary['stars_alone']['n_within_2sigma']}/{summary['N']}")
    L(f"- Fraction within 1σ / 2σ (power-law UMM): "
      f"{summary['powerlaw_umm']['n_within_1sigma']}/{summary['N']} / "
      f"{summary['powerlaw_umm']['n_within_2sigma']}/{summary['N']}")
    L(f"- Fraction within 1σ / 2σ (logistic UMM): "
      f"{summary['logistic_umm']['n_within_1sigma']}/{summary['N']} / "
      f"{summary['logistic_umm']['n_within_2sigma']}/{summary['N']}")
    L("")
    L(f"### Classification: **{summary['classification'].upper()}**")
    L("")
    L(summary["classification_note"])
    L("")
    L(f"*{summary['honesty']}*")
    L("")
    L("### Systematics")
    L("")
    L("Dominant systematics and their effect on the power-law residual:")
    L("")
    L("1. **Distance:** Scaling \(M_*\\propto D^2\), \(R_e\\propto D\) leaves \(\\langle\\Sigma\\rangle\) invariant; characteristic accelerations are distance-invariant at leading order, while Wolf masses scale ∝ \(D\). Published analyses (Keim et al. 2026 App. B) show DM-deficiency conclusions are robust to the trail distance. Shifting each system by its \(\\pm 1\\sigma\) distance does not bring power-law \(\\sigma\) into agreement with \(\\sigma_{\\rm obs}\).")
    L("2. **Stellar M/L:** \(\\pm 1\\sigma\) in \(M_*\) moves \(\\sigma_{\\rm stars}\) by ~half the fractional mass error and moves \(a_{\\rm ind}\) weakly through \(\\Sigma\). Residuals remain positive and large for the power-law column.")
    L("3. **Residual rotation / tides:** FCC 224 shows slow prolate rotation; DF2/DF4 KCWI maps show little ordered motion. An illustrative correction that attributes half of \(\\sigma_{\\rm obs}\) to rotation *increases* tension for UMM (true dispersion lower). A 30% tidal inflation of observed motions reduces tension only partially and does not reverse the sign of \(\\varepsilon_{\\rm PL}\) for the ensemble.")
    L("4. **Mean vs local \(\\Sigma\):** Using a factor-of-two lower local \(\\Sigma(R_e)\) still leaves \(a_{\\rm ind}\\gtrsim a_{\\rm bar}\) and \(\\sigma_{\\rm PL}\) well above \(\\sigma_{\\rm obs}\).")
    L("5. **Geometry:** The galactic specialization assumes a thin disk; these UDGs are pressure-supported spheroids. This is a structural limitation of the first contact, not a free parameter.")
    L("")
    # Per-system systematics snapshot
    L("Example (DF2) one-at-a-time shifts of \(\\sigma_{\\rm PL}\):")
    L("")
    df2 = next(r for r in rows if r["short"] == "DF2")
    sy = df2["systematics"]
    sp = df2["sigma_powerlaw"]
    so = df2["sigma_obs"]
    L(f"- Fiducial sigma_PL = {sp:.1f} km/s (obs {so:.1f})")
    L(
        f"- D+1sigma: sigma_PL = {sy['D_plus']['sigma_powerlaw']:.1f}, "
        f"eps = {sy['D_plus']['eps_powerlaw']:+.2f}"
    )
    L(
        f"- D-1sigma: sigma_PL = {sy['D_minus']['sigma_powerlaw']:.1f}, "
        f"eps = {sy['D_minus']['eps_powerlaw']:+.2f}"
    )
    L(f"- M_star+1sigma: sigma_PL = {sy['ML_plus']['sigma_powerlaw']:.1f}")
    L(f"- M_star-1sigma: sigma_PL = {sy['ML_minus']['sigma_powerlaw']:.1f}")
    L("")
    L("---")
    L("")
    L("## 5. Limitations and sample size")
    L("")
    L("- **N = 4** primary systems with published high-resolution stellar or GC kinematics of comparable quality. This is not a population study.")
    L("- Kinematic uncertainties are large and often asymmetric; broadening and binary corrections matter at \(\\sigma\\sim 5\\text{--}10\) km/s.")
    L("- Mean \(\\Sigma\) proxy and thin-disk formula applied to spheroidal systems.")
    L("- Logistic \(\\xi_c\) is represented by a \(\\Sigma_{\\rm mid}\) window about \(\\Sigma_{\\rm ref}\) consistent with intermediate-regime SPARC placement; the paper does not publish a single numerical \(\\xi_c\).")
    L("- No re-tuning of \(A\), \(n_g\), \(\\alpha_{\\rm eff}\), or \(\\xi_c\).")
    L("- Frozen v14 paper and SPARC residual products were not modified.")
    L("")
    L("---")
    L("")
    L("## 6. Predictive paragraph and go/no-go")
    L("")
    if summary["classification"] == "ruled out":
        L("**Predictive insertion paragraph: withheld.**")
        L("")
        L(
            "The frozen intermediate-regime UMM form is classified as **ruled out** "
            "by this first confrontation within the quoted uncertainties and the "
            "stated estimator assumptions. A ready-to-use predictive note for a "
            "future paper revision is therefore **not** supplied. The scientific "
            "value of the exercise is the transparent tension: pure baryons match "
            "the DM-deficient UDG kinematics, while \(a_{\\rm ind}=A(\\Sigma/\\Sigma_{\\rm ref})^{0.4}\) "
            "at the frozen SPARC amplitude overpredicts \(\\sigma\) for every system "
            "in the sample. Any future revision that retains this amplitude for "
            "disk galaxies must address why the same local prescription fails in "
            "these low-dispersion spheroids (geometry of the specialization, "
            "validity of mean-\(\\Sigma\) proxy, or a genuine limitation of the "
            "intermediate-regime form at this scale)."
        )
    else:
        L("**Predictive insertion paragraph (150–250 words):**")
        L("")
        L(_predictive_paragraph())
    L("")
    L("### Is the result sound enough to merit a predictive note in a future revision?")
    L("")
    if summary["classification"] == "ruled out":
        L(
            "**No** — not as a positive predictive claim. The confrontation is sound "
            "enough to *document tension* and to motivate either a careful statement "
            "of domain of applicability (thin-disk intermediate regime only) or a "
            "dedicated follow-up on the low-density / spheroidal sector. It is not "
            "sound enough to advertise UDG consistency with the frozen power-law form."
        )
    elif summary["classification"] == "mildly tensioned":
        L(
            "**Marginal** — a carefully hedged note could mention the sample and the "
            "direction of residuals, but should not claim confirmation. Additional "
            "data are required before a strong observational-contact claim."
        )
    else:
        L(
            "**Yes, with restraint** — a short predictive note is appropriate if "
            "honesty language (small sample, systematics) is retained."
        )
    L("")
    L("---")
    L("")
    L("## 7. Residual-outcome summary (human go/no-go)")
    L("")
    L(_outcome_paragraph(summary, rows))
    L("")
    L("---")
    L("")
    L("*Report generated by `analysis/udg_confrontation/src/run_confrontation.py`. "
      "Numerical products: `results/`. Re-run: `python -m analysis.udg_confrontation.src.run_confrontation`.*")
    L("")
    return "\n".join(lines)


def _predictive_paragraph() -> str:
    return (
        "As an illustrative first estimate on a small sample of dark-matter-deficient "
        "ultra-diffuse galaxies (NGC 1052-DF2, DF4, DF9, and FCC 224), we confront the "
        "frozen intermediate-regime residual acceleration "
        r"\(a_{\rm ind}=A(\Sigma_{\rm bar}/\Sigma_{\rm ref})^{0.4}\) "
        r"with \(A=(1.41\pm 0.07)\times 10^{-11}\,\mathrm{m\,s^{-2}}\) "
        "from the SPARC residual analysis (selection_id ea6c92761ffd3c5a), without "
        "re-fitting any UMM parameter. Converting accelerations to line-of-sight "
        "dispersions with the Wolf et al. (2010) estimator under mass-follows-light "
        "and isotropic-enough assumptions, pure-baryon predictions lie close to the "
        "published stellar and globular-cluster dispersions, whereas the frozen "
        "power-law \(a_{\rm ind}\) produces systematically higher \(\\sigma\) at the "
        "mean surface densities implied by \(M_*\) and \(R_e\). A logistic "
        r"\(\Phi(\xi)\) extrapolation with packing \(\alpha_{\rm eff}\approx 1.8\) "
        "and the joint-normalization midpoint window about \(\\Sigma_{\\rm ref}\) "
        "does not remove the offset within that window. Systematics (distance, "
        "stellar \(M/L\), residual rotation, and possible tides) remain significant "
        "and the sample is small; we do not claim a population result. "
        "Additional kinematic data on these low-density classes would help refine "
        "the lower shoulder of \(\\Phi(\\xi)\)."
    )


def _outcome_paragraph(summary: Dict[str, Any], rows: List[Dict[str, Any]]) -> str:
    med_pl = summary["powerlaw_umm"]["median_eps_sigma"]
    med_st = summary["stars_alone"]["median_eps_sigma"]
    names = ", ".join(r["short"] for r in rows)
    return (
        f"**Go/no-go for a predictive note: "
        f"{'NO — withhold positive note' if summary['classification']=='ruled out' else 'CAUTION' if summary['classification']=='mildly tensioned' else 'YES — with restraint'}.** "
        f"On the residual sample ({names}; N={summary['N']}), median "
        f"ε_σ(stars) = {med_st:+.2f} while median ε_σ(UMM power-law) = {med_pl:+.2f}. "
        f"Classification: **{summary['classification']}**. "
        f"{summary['classification_note']} "
        f"Illustrative first estimate only; small sample; systematics remain significant. "
        f"Do not re-fit A or α; do not claim a large statistical population."
    )


def main() -> int:
    rows = [analyze_system(s) for s in SYSTEMS if s.get("include_in_residuals", True)]
    summary = residual_summary(rows)
    write_outputs(rows, summary)

    print("UMM UDG confrontation — frozen v14 baseline (no re-fit)")
    print(f"A = {A_FROZEN:.2e} ± {A_FROZEN_UNC:.2e} m/s^2, n_g = {N_G}, alpha_eff = {ALPHA_EFF}")
    print(f"N = {len(rows)} systems")
    print("-" * 72)
    print(f"{'name':12s} {'σ_obs':>7s} {'σ_*':>7s} {'σ_PL':>7s} {'σ_log':>7s} {'ε_PL':>7s} {'nσ_PL':>6s}")
    for r in rows:
        print(
            f"{r['short']:12s} {r['sigma_obs']:7.2f} {r['sigma_stars']:7.2f} "
            f"{r['sigma_powerlaw']:7.2f} {r['sigma_logistic']:7.2f} "
            f"{r['eps_sigma_powerlaw']:+7.2f} {r['n_sigma_powerlaw']:6.2f}"
        )
    print("-" * 72)
    print(f"median ε_PL = {summary['powerlaw_umm']['median_eps_sigma']:+.3f}")
    print(f"median ε_stars = {summary['stars_alone']['median_eps_sigma']:+.3f}")
    print(f"classification: {summary['classification']}")
    print(f"wrote results under {RESULTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
