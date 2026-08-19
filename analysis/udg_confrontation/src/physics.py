"""Frozen-v14 UMM predictors for UDG kinematics (no re-fitting).

All amplitudes and indices are fixed to the v14 / SPARC residual baseline:
  A = (1.41 ± 0.07) × 10^{-11} m s^{-2}
  n_g = 0.4
  Σ_ref = 1 M_⊙ pc^{-2}
  selection_id ea6c92761ffd3c5a, Υ_* = 0.5

Unit conventions match analysis/sparc_residuals/physics.py.
"""
from __future__ import annotations

import math
from typing import Dict, Optional, Tuple

# --- Frozen v14 / SPARC residual baseline (do not re-fit) ---
A_FROZEN = 1.41e-11  # m s^{-2}
A_FROZEN_UNC = 0.07e-11  # m s^{-2} (1σ on global SPARC fit)
N_G = 0.4
SIGMA_REF = 1.0  # M_⊙ pc^{-2}
ALPHA_EFF = 1.8  # packing-derived ≈1.8 (preferred band 1.75±0.09)
# Joint-normalization window for logistic midpoint surface density (M_⊙ pc^{-2}).
# Midpoint of Φ(ξ) is placed near SPARC-relevant intermediate densities about Σ_ref
# without introducing new free parameters beyond the paper's joint-norm window.
SIGMA_MID_WINDOW = (0.3, 1.0, 3.0)  # low, central, high of ξ_c window via Σ_mid

G_KPC = 4.30091e-6  # kpc M_⊙^{-1} (km/s)^2
# a [m s^{-2}] → (km/s)^2 kpc^{-1}  (same factor as sparc_residuals/physics.py)
A_SI_TO_ASTRO = 3.085677581e13


def a_si_to_astro(a_si: float) -> float:
    """Convert acceleration from m s^{-2} to (km/s)^2 kpc^{-1}."""
    return a_si * A_SI_TO_ASTRO


def mean_sigma_within_re(m_star: float, r_e_kpc: float) -> float:
    """Mean stellar surface density within R_e [M_⊙ pc^{-2}].

    Half the stellar mass lies within the projected half-light radius, so
        ⟨Σ⟩_e = (M_*/2) / (π R_e^2) = M_* / (2 π R_e^2)
    with R_e in pc. Gas is neglected for these quiescent systems (published
    HI non-detections or negligible cold gas). This is a *mean* density, not a
    published radial profile — limitation is stated at call sites.
    """
    if m_star <= 0 or r_e_kpc <= 0:
        return 0.0
    r_e_pc = r_e_kpc * 1000.0
    return m_star / (2.0 * math.pi * r_e_pc * r_e_pc)


def a_ind_powerlaw(
    sigma_bar: float,
    A: float = A_FROZEN,
    n_g: float = N_G,
    sigma_ref: float = SIGMA_REF,
) -> float:
    """Intermediate-regime a_ind = A (Σ_bar/Σ_ref)^{n_g} in (km/s)^2 kpc^{-1}."""
    if sigma_bar <= 0 or A <= 0:
        return 0.0
    a_si = A * (sigma_bar / sigma_ref) ** n_g
    return a_si_to_astro(a_si)


def a_ind_logistic(
    sigma_bar: float,
    A: float = A_FROZEN,
    alpha_eff: float = ALPHA_EFF,
    sigma_mid: float = SIGMA_REF,
    sigma_ref: float = SIGMA_REF,
) -> float:
    """Logistic Φ(ξ) extrapolation with packing α_eff and midpoint at Σ_mid.

    Construction (no new free parameters beyond the joint-norm window):
      ξ − ξ_c = (1/2) ln(Σ / Σ_mid)
      f = Φ/Φ_max = 1 / (1 + exp(−α_eff (ξ − ξ_c)))
      a_ind = a_mid * f / f_mid,  f_mid = 1/2 at Σ = Σ_mid

    Amplitude is fixed so that at Σ = Σ_ref the intermediate power-law value
    A (Σ_ref/Σ_ref)^{n_g} = A is recovered when Σ_mid = Σ_ref:
      a_ind(Σ_mid) = A * (Σ_mid / Σ_ref)^{n_g}   [match local power-law at mid]
    and a_ind = 2 * a_mid * f  so a_ind(Σ_mid) = a_mid.

    When Σ_mid = Σ_ref: a_ind(Σ_ref) = A, and the local logarithmic slope at
    midpoint is α_eff/4 ≈ 0.45, consistent with packing n_g.
    """
    if sigma_bar <= 0 or A <= 0 or sigma_mid <= 0:
        return 0.0
    # Amplitude at logistic midpoint: match power-law (frozen n_g) at Σ_mid
    a_mid_si = A * (sigma_mid / sigma_ref) ** N_G
    # f = 1 / (1 + (Σ_mid/Σ)^{α_eff/2})
    ratio = (sigma_mid / sigma_bar) ** (0.5 * alpha_eff)
    f = 1.0 / (1.0 + ratio)
    # At midpoint f=1/2 → a = a_mid; scale as a = 2 a_mid f
    a_si = 2.0 * a_mid_si * f
    return a_si_to_astro(a_si)


def a_bar_half_light(m_star: float, r_e_kpc: float) -> float:
    """Characteristic Newtonian acceleration from stars within R_e.

    Uses half the stellar mass as the mass interior to R_e (by definition of
    half-light radius for mass-follows-light):
        a_bar = G (M_*/2) / R_e^2
    in (km/s)^2 kpc^{-1}.
    """
    if m_star <= 0 or r_e_kpc <= 0:
        return 0.0
    m_half = 0.5 * m_star
    return G_KPC * m_half / (r_e_kpc ** 2)


def wolf_mass_from_sigma(sigma_kms: float, r_e_kpc: float) -> float:
    """Wolf et al. (2010) half-light dynamical mass [M_⊙].

        M_{1/2} = 4 σ_los^2 R_e / G

    Equivalent form: M(r_{1/2}) = 3 σ^2 r_{1/2}/G with r_{1/2} ≈ (4/3) R_e.
    Assumptions: pressure-supported system, flat or mildly varying σ profile,
    mass-follows-light for the comparison sample, isotropic-enough orbits for
    the Wolf estimator (weak anisotropy sensitivity near r_{1/2}).
    """
    if sigma_kms <= 0 or r_e_kpc <= 0:
        return 0.0
    return 4.0 * (sigma_kms ** 2) * r_e_kpc / G_KPC


def wolf_sigma_from_mass(m_half: float, r_e_kpc: float) -> float:
    """Invert Wolf estimator: σ_los from M_{1/2} and R_e [km/s]."""
    if m_half <= 0 or r_e_kpc <= 0:
        return 0.0
    return math.sqrt(G_KPC * m_half / (4.0 * r_e_kpc))


def sigma_from_a_tot(a_tot_astro: float, r_e_kpc: float) -> float:
    """Line-of-sight dispersion from total acceleration at R_e via Wolf.

    Identify a_tot with G M_{1/2} / R_e^2 so that
        σ = sqrt(a_tot * R_e / 4).
    """
    if a_tot_astro <= 0 or r_e_kpc <= 0:
        return 0.0
    return math.sqrt(a_tot_astro * r_e_kpc / 4.0)


def predict_sigma_three_models(
    m_star: float,
    r_e_kpc: float,
    A: float = A_FROZEN,
    n_g: float = N_G,
    alpha_eff: float = ALPHA_EFF,
    sigma_mid: float = SIGMA_REF,
    sigma_bar: Optional[float] = None,
) -> Dict[str, float]:
    """Three frozen-v14 predictions for σ_los at R_e.

    (a) stars alone
    (b) stars + intermediate power-law a_ind
    (c) stars + logistic Φ(ξ) at packing α_eff and chosen Σ_mid (ξ_c window)

    Returns accelerations [astro units], Σ, and σ [km/s] for each model.
    """
    if sigma_bar is None:
        sigma_bar = mean_sigma_within_re(m_star, r_e_kpc)
    a_bar = a_bar_half_light(m_star, r_e_kpc)
    a_pl = a_ind_powerlaw(sigma_bar, A=A, n_g=n_g)
    a_log = a_ind_logistic(
        sigma_bar, A=A, alpha_eff=alpha_eff, sigma_mid=sigma_mid
    )

    sig_a = sigma_from_a_tot(a_bar, r_e_kpc)
    sig_b = sigma_from_a_tot(a_bar + a_pl, r_e_kpc)
    sig_c = sigma_from_a_tot(a_bar + a_log, r_e_kpc)

    m_half = 0.5 * m_star
    return {
        "sigma_bar_msun_pc2": sigma_bar,
        "a_bar": a_bar,
        "a_ind_powerlaw": a_pl,
        "a_ind_logistic": a_log,
        "sigma_stars": sig_a,
        "sigma_powerlaw": sig_b,
        "sigma_logistic": sig_c,
        "m_half": m_half,
        "m_dyn_stars": wolf_mass_from_sigma(sig_a, r_e_kpc),
        "m_dyn_powerlaw": wolf_mass_from_sigma(sig_b, r_e_kpc),
        "m_dyn_logistic": wolf_mass_from_sigma(sig_c, r_e_kpc),
    }


def fractional_residual_sigma(sigma_obs: float, sigma_model: float) -> float:
    """ε_σ = (σ_model − σ_obs) / σ_obs  (signed; positive = model high)."""
    if sigma_obs is None or sigma_obs <= 0:
        return float("nan")
    return (sigma_model - sigma_obs) / sigma_obs


def fractional_residual_mass(m_obs: float, m_model: float) -> float:
    """ε_M = (M_model − M_obs) / M_obs."""
    if m_obs is None or m_obs <= 0:
        return float("nan")
    return (m_model - m_obs) / m_obs
