"""Machine-readable literature manifest for DM-deficient UDG sample.

Every adopted number has value, 1σ uncertainty, unit, citation key, and notes.
Only published peer-reviewed / arXiv-vetted sources. No private re-reductions.
No surface-density *profiles* invented — mean Σ within R_e is computed from
M_* and R_e and flagged as a mean, not a published profile.
"""
from __future__ import annotations

from typing import Any, Dict, List

# ---------------------------------------------------------------------------
# Frozen UMM inputs (not fitted here)
# ---------------------------------------------------------------------------
FROZEN_UMM: Dict[str, Any] = {
    "A_m_s2": 1.41e-11,
    "A_unc_m_s2": 0.07e-11,
    "n_g": 0.4,
    "Sigma_ref_msun_pc2": 1.0,
    "Upsilon_star_sparc": 0.5,
    "selection_id": "ea6c92761ffd3c5a",
    "alpha_eff": 1.8,
    "alpha_eff_band": [1.66, 1.84],
    "Sigma_mid_window_msun_pc2": [0.3, 1.0, 3.0],
    "sources": {
        "A": "paper/UMM_Paper_Draft_v14.tex Sec. SPARC; analysis/sparc_residuals/results/",
        "alpha_eff": "notes/UMM_Packing_Coefficient_Tightening.md; v14 packing §",
        "Sigma_ref": "v14 galactic specialization, Σ_ref = 1 M_⊙ pc^{-2}",
    },
}

# ---------------------------------------------------------------------------
# Primary residual sample
# ---------------------------------------------------------------------------
SYSTEMS: List[Dict[str, Any]] = [
    {
        "name": "NGC1052-DF2",
        "short": "DF2",
        "include_in_residuals": True,
        "distance": {
            "value": 20.0,
            "unc": 1.7,
            "unit": "Mpc",
            "method": "SBF (homogeneous group scale)",
            "citation": "van Dokkum et al. 2018, Nature, 555, 629 (arXiv:1803.10237)",
            "notes": (
                "SBF D=19.0±1.7 Mpc; paper adopts ≈20 Mpc with NGC 1052. "
                "TRGB re-measurements exist (Danieli et al. 2020 ≈22 Mpc); "
                "we adopt the homogeneous 20 Mpc group scale and propagate ±1.7 Mpc. "
                "Shen et al. 2023 quote 21.7 Mpc in some comparisons."
            ),
        },
        "M_star": {
            "value": 2.0e8,
            "unc": 0.4e8,
            "unit": "Msun",
            "citation": "van Dokkum et al. 2018, Nature; M/L_V≈2.0",
            "notes": "M_stars ≈ 2×10^8 M_⊙ at D=20 Mpc; ~25% M/L systematic folded into unc.",
        },
        "R_e": {
            "value": 2.2,
            "unc": 0.1,
            "unit": "kpc",
            "citation": "van Dokkum et al. 2018, Nature (R_e=22.6'' → 2.2 kpc at 20 Mpc)",
            "notes": "Major-axis half-light radius; Sérsic n=0.6. Scales ∝ D.",
        },
        "sigma_obs": {
            "value": 8.5,
            "unc_plus": 2.3,
            "unc_minus": 3.1,
            "unit": "km/s",
            "tracer": "stars (KCWI integrated light)",
            "citation": "Danieli et al. 2019, ApJL, 874, L12 (arXiv:1901.03711)",
            "notes": (
                "Primary stellar σ. GC-based values: van Dokkum et al. 2018 "
                "σ_intr=3.2^{+5.5}_{-3.2} (90% UL 10.5); revised GC analyses exist. "
                "Keim et al. 2026 quote a broadening-corrected 6.1^{+3.7}_{-3.1} for cross-check; "
                "we retain Danieli 2019 as the published primary stellar measurement."
            ),
        },
        "sigma_baryon_lit": {
            "value": 7.0,
            "unc": 2.0,
            "unit": "km/s",
            "citation": "Danieli et al. 2019; van Dokkum et al. 2018 (stars alone ~7–8 km/s)",
            "notes": "Literature pure-baryon expectation; we recompute via Wolf below.",
        },
        "profile_limitation": (
            "No published high-resolution mass surface-density profile used. "
            "Mean Σ within R_e computed from M_* and R_e (mass-follows-light)."
        ),
    },
    {
        "name": "NGC1052-DF4",
        "short": "DF4",
        "include_in_residuals": True,
        "distance": {
            "value": 20.0,
            "unc": 1.6,
            "unit": "Mpc",
            "method": "SBF / TRGB (homogeneous group scale)",
            "citation": "van Dokkum et al. 2019, ApJL, 874, L5 (SBF 19.9±2.8); Danieli et al. 2020 TRGB ≈20.0",
            "notes": "Adopt D=20 Mpc on the same group scale as DF2.",
        },
        "M_star": {
            "value": 1.5e8,
            "unc": 0.4e8,
            "unit": "Msun",
            "citation": "van Dokkum et al. 2019, ApJL (M/L_V=2.0±0.5; L_V=(7.7±0.8)×10^7 L_⊙)",
            "notes": "M_stars=(1.5±0.4)×10^8 M_⊙ at D=20 Mpc.",
        },
        "R_e": {
            "value": 1.6,
            "unc": 0.1,
            "unit": "kpc",
            "citation": "van Dokkum et al. 2019, ApJL (R_e=1.6 kpc; n=0.79)",
            "notes": "Scales ∝ D. Shen et al. 2023 quote 16.5'' = 1.6 kpc.",
        },
        "sigma_obs": {
            "value": 8.0,
            "unc_plus": 2.3,
            "unc_minus": 1.9,
            "unit": "km/s",
            "tracer": "stars (KCWI integrated light)",
            "citation": "Shen et al. 2023, ApJ, 958, 3 (arXiv:2309.08592)",
            "notes": (
                "Stellar σ_stars=8.0^{+2.3}_{-1.9} km/s. Combined star+GC fiducial "
                "σ_f=6.3^{+2.5}_{-1.6} km/s also reported; we use stellar as primary. "
                "GC-only van Dokkum et al. 2019: 4.2^{+4.4}_{-2.2} km/s."
            ),
        },
        "sigma_baryon_lit": {
            "value": 7.0,
            "unc": 1.0,
            "unit": "km/s",
            "citation": "Shen et al. 2023 (stars alone 7±1 km/s)",
            "notes": "Literature pure-baryon expectation.",
        },
        "profile_limitation": (
            "No published radial Σ profile used; mean Σ within R_e from M_* and R_e."
        ),
    },
    {
        "name": "NGC1052-DF9",
        "short": "DF9",
        "include_in_residuals": True,
        "distance": {
            "value": 20.6,
            "unc": 2.0,
            "unit": "Mpc",
            "method": "Trail-position distance (group kinematic trail)",
            "citation": "Keim et al. 2026, ApJ (arXiv:2603.15860); trail scale Keim et al. 2025",
            "notes": (
                "Expected distance from position along the NGC 1052 trail. "
                "Results are weakly distance-dependent (σ_baryon ∝ D^{1/2} scaling "
                "cancels in mass comparison at leading order; see Keim et al. App. B)."
            ),
        },
        "M_star": {
            "value": 1.4e8,
            "unc": 0.3e8,
            "unit": "Msun",
            "citation": "Keim et al. 2026 (from Gannon et al. 2023 photometry; M/L_g=2.2)",
            "notes": "Total stellar mass 1.4×10^8 M_⊙; half-light stellar mass 0.71±0.16×10^8.",
        },
        "R_e": {
            "value": 1.1,
            "unc": 0.1,
            "unit": "kpc",
            "citation": "Keim et al. 2026; Gannon et al. 2023 (11.1'' at 20.6 Mpc → 1.1 kpc)",
            "notes": "Circularized effective radius used with Wolf estimator.",
        },
        "sigma_obs": {
            "value": 6.4,
            "unc_plus": 4.0,
            "unc_minus": 4.3,
            "unit": "km/s",
            "tracer": "stars (KCWI; broadening- and binary-corrected)",
            "citation": "Keim et al. 2026, ApJ (arXiv:2603.15860)",
            "notes": (
                "Gravitational σ after subtracting σ_broadening=5.4±2 and "
                "σ_binaries=2.3±0.6 from fit σ_fit=8.7^{+3.4}_{-3.7}."
            ),
        },
        "sigma_baryon_lit": {
            "value": 8.3,
            "unc_plus": 0.9,
            "unc_minus": 1.4,
            "unit": "km/s",
            "citation": "Keim et al. 2026 (Wolf from M_e,* alone)",
            "notes": "Literature pure-baryon expectation from stellar mass within R_e.",
        },
        "profile_limitation": (
            "No published Σ profile; mean Σ within R_e from M_* and R_e."
        ),
    },
    {
        "name": "FCC224",
        "short": "FCC224",
        "include_in_residuals": True,
        "distance": {
            "value": 20.0,
            "unc": 1.5,
            "unit": "Mpc",
            "method": "Fornax cluster distance (SBF/PNLF scale)",
            "citation": "Buzzo et al. 2025, A&A (arXiv:2502.05405); Tang et al. 2025a photometry at 20 Mpc",
            "notes": (
                "Fornax outskirts; systemic velocity 1405±3 km/s consistent with cluster. "
                "GCLF and stellar mass quoted at 20 Mpc. Distance unc. ~typical Fornax ladder."
            ),
        },
        "M_star": {
            "value": 1.74e8,
            "unc": 0.16e8,
            "unit": "Msun",
            "citation": "Buzzo et al. 2025; Tang et al. 2025a: log(M_*/M_⊙)=8.24±0.04",
            "notes": "10^{8.24}=1.74×10^8 M_⊙.",
        },
        "R_e": {
            "value": 1.89,
            "unc": 0.05,
            "unit": "kpc",
            "citation": "Buzzo et al. 2025; Tang et al. 2025a: R_e=1.89±0.01 kpc",
            "notes": "UDG within uncertainties; μ_g,0=23.97±0.03; b/a=0.64; n=0.75.",
        },
        "sigma_obs": {
            "value": 7.82,
            "unc_plus": 6.74,
            "unc_minus": 4.36,
            "unit": "km/s",
            "tracer": "stars (KCWI red-arm CaT; broadening-corrected)",
            "citation": "Buzzo et al. 2025, A&A (arXiv:2502.05405)",
            "notes": (
                "Published σ_stars=7.82^{+6.74}_{-4.36} km/s after subtracting "
                "σ_broadening=6.9±2.6 from raw σ_galaxy=10.43±5.76 km/s "
                "(Buzzo et al. 2025 §3.3.1). Asymmetric 1σ adopted as published; "
                "not a fabricated symmetric error. Galaxy is DM-deficient within 1 R_e "
                "per authors. Slow prolate rotation 7.5±3.0 km/s reported "
                "(not subtracted here; see systematics)."
            ),
        },
        "sigma_baryon_lit": {
            "value": 7.8,
            "unc": 2.0,
            "unit": "km/s",
            "citation": "Buzzo et al. 2025 (stars-alone Jeans/SMHM comparison value)",
            "notes": "Authors' stars-alone prediction 7.8±2.0 km/s.",
        },
        "profile_limitation": (
            "No published full Σ profile used; mean Σ within R_e from M_* and R_e. "
            "KCWI FoV covers only ~0.4 R_e,circ for population gradients."
        ),
    },
]

# ---------------------------------------------------------------------------
# Future targets (no residual entry — kinematics absent or not comparable)
# ---------------------------------------------------------------------------
FUTURE_TARGETS: List[Dict[str, Any]] = [
    {
        "name": "NGC1052 trail members (excl. DF2/DF4/DF9)",
        "reason": (
            "Keim et al. 2025/2026: remaining trail galaxies up to ~100× fainter; "
            "radial velocities partially available but internal kinematics not yet "
            "of KCWI/MUSE quality comparable to DF2/DF4/DF9."
        ),
        "citation": "van Dokkum et al. 2022a; Keim et al. 2025, 2026",
    },
    {
        "name": "FCC240",
        "reason": (
            "Buzzo et al. 2026 (arXiv:2605.24099): close companion to FCC 224; "
            "analogue discussion; published kinematics of residual-sample quality "
            "not adopted here as primary."
        ),
        "citation": "Buzzo et al. 2026, ApJ (arXiv:2605.24099)",
    },
    {
        "name": "Tidal-dwarf candidates (e.g. NGC 5291N, etc.)",
        "reason": (
            "Classic TDGs can be DM-poor but generally have rotation-supported HI "
            "kinematics and different formation channels; not included without a "
            "homogeneous pressure-supported σ of UDG quality."
        ),
        "citation": "Lelli et al. 2015; Bournaud et al. 2007 (context only)",
    },
    {
        "name": "Almost-dark HI clouds / HI-bearing UDGs",
        "reason": (
            "Interesting for low-Σ shoulder of Φ(ξ) but require published stellar or "
            "GC dispersions of comparable quality; list as future targets only."
        ),
        "citation": "—",
    },
]


def sigma_obs_1sigma(system: Dict[str, Any]) -> float:
    """Symmetric 1σ proxy = max of +/− published uncertainties."""
    s = system["sigma_obs"]
    return float(max(s.get("unc_plus", s.get("unc", 0.0)), s.get("unc_minus", s.get("unc", 0.0))))
