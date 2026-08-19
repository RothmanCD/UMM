"""Physical helpers for UMM SPARC residual pipeline (design note §4–§7)."""
from __future__ import annotations

import math
from typing import Iterable, List, Optional, Sequence, Tuple

# Constants
G_KPC = 4.30091e-6  # kpc M_sun^-1 (km/s)^2
SIGMA_REF = 1.0  # M_sun pc^-2
N_G_PRIMARY = 0.4
UPSILON_DISK_FIDUCIAL = 0.5
UPSILON_BULGE_FIDUCIAL = 0.7  # SPARC-preferred bulge when disk is 0.5


def v_bar_squared(
    vgas: float,
    vdisk: float,
    vbul: float,
    upsilon_disk: float,
    upsilon_bulge: float,
) -> float:
    """Newtonian baryonic V^2 from SPARC mass-model components (Υ=1 tables).

    SPARC convention: Vgas may be negative; use Vgas*|Vgas|.
    """
    return vgas * abs(vgas) + upsilon_disk * (vdisk ** 2) + upsilon_bulge * (vbul ** 2)


def v_bar(
    vgas: float,
    vdisk: float,
    vbul: float,
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
) -> float:
    vb2 = v_bar_squared(vgas, vdisk, vbul, upsilon_disk, upsilon_bulge)
    return math.sqrt(max(vb2, 0.0))


def sigma_star_msun_pc2(
    sb_disk: float,
    sb_bul: float,
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
) -> float:
    """Stellar mass surface density from 3.6μm SB (L_sun pc^-2) and M/L."""
    return upsilon_disk * max(sb_disk, 0.0) + upsilon_bulge * max(sb_bul, 0.0)


def sigma_gas_from_vgas(vgas: float, r_kpc: float) -> float:
    """Approximate HI+He mass surface density from Vgas via thin-sheet proxy.

    Σ ≈ Vgas^2 / (2 π G R) with unit conversion to M_sun pc^-2.
    Uses |Vgas| so density is non-negative; gas_source labeled SPARC_HI.
    """
    if r_kpc <= 0:
        return 0.0
    # a = V^2/R in (km/s)^2 / kpc; Σ = a/(2πG) with R in pc
    # Σ [Msun/pc2] = V^2 / (2π G_pc * R_pc), G_pc = 4.30091e-3
    g_pc = 4.30091e-3
    r_pc = r_kpc * 1000.0
    return (vgas ** 2) / (2.0 * math.pi * g_pc * r_pc)


def sigma_bar(
    sb_disk: float,
    sb_bul: float,
    vgas: float,
    r_kpc: float,
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
) -> float:
    return sigma_star_msun_pc2(sb_disk, sb_bul, upsilon_disk, upsilon_bulge) + sigma_gas_from_vgas(
        vgas, r_kpc
    )


def a_ind(sigma_bar_val: float, A: float, n_g: float = N_G_PRIMARY) -> float:
    """Independent acceleration a_ind = A (Σ_bar/Σ_ref)^n_g  [km^2 s^-2 kpc^-1]."""
    if sigma_bar_val <= 0 or A <= 0:
        return 0.0
    # A is in m s^-2; convert to (km/s)^2 / kpc for combining with V^2/R.
    # a_si = V_si^2 / R_si with V_si = V_kms*1000, R_si = R_kpc * 3.085677581e19 m
    # ⇒ V_kms^2 / R_kpc = a_si * 3.085677581e19 / 1e6 = a_si * 3.085677581e13
    a_si = A * (sigma_bar_val / SIGMA_REF) ** n_g
    return a_si * 3.085677581e13


def v_model_from_a(r_kpc: float, a_tot: float) -> float:
    if r_kpc <= 0 or a_tot <= 0:
        return 0.0
    return math.sqrt(a_tot * r_kpc)


def a_from_v(v: float, r_kpc: float) -> float:
    if r_kpc <= 0:
        return 0.0
    return (v ** 2) / r_kpc


def r_out_threshold(r_min: float, r_last: float) -> float:
    """R_out = R_min + 0.70 (R_last - R_min)."""
    return r_min + 0.70 * (r_last - r_min)


def outer_mean_abs_frac_residual(
    r: Sequence[float],
    v_obs: Sequence[float],
    v_model: Sequence[float],
    r_out: float,
) -> Tuple[float, int]:
    """ε = mean |Vobs-Vmod|/Vobs for R≥R_out, unweighted. Drops Vobs≤0."""
    vals: List[float] = []
    for ri, vo, vm in zip(r, v_obs, v_model):
        if ri < r_out:
            continue
        if vo is None or vo <= 0:
            continue
        vals.append(abs(vo - vm) / vo)
    if not vals:
        return float("nan"), 0
    return sum(vals) / len(vals), len(vals)


# --- NFW ---

def dutton_maccio_c200(m200: float, h: float = 0.7) -> float:
    """Dutton & Macciò (2014)-style mass–concentration (Planck-ish).

    log10 c = 0.905 - 0.101 log10(M200 / (1e12 h^-1 Msun))
    """
    if m200 <= 0:
        return 10.0
    m12 = m200 / (1.0e12 / h)  # M200 in units of 1e12 h^-1 Msun
    # M200 [Msun]; M200 / (1e12 * Msun/h) = M200 * h / 1e12
    m12 = m200 * h / 1.0e12
    logc = 0.905 - 0.101 * math.log10(max(m12, 1e-6))
    return 10.0 ** logc


def r200_from_m200(m200: float, h: float = 0.7) -> float:
    """r_200 [kpc] for spherical overdensity 200 ρ_crit."""
    # ρ_crit = 3 H0^2 / (8πG); H0 = 100h km/s/Mpc
    # r_200^3 = M200 / ( (4π/3) * 200 * ρ_crit )
    H0 = 100.0 * h  # km/s/Mpc
    # G_astro = 4.30091e-6 kpc/Msun (km/s)^2
    # ρ_crit in Msun/kpc^3: 3 H0^2 / (8πG) with H0 in km/s/kpc
    H0_kpc = H0 / 1000.0  # km/s/kpc
    rho_crit = 3.0 * H0_kpc ** 2 / (8.0 * math.pi * G_KPC)
    r200 = (m200 / ((4.0 * math.pi / 3.0) * 200.0 * rho_crit)) ** (1.0 / 3.0)
    return r200


def v_nfw(r_kpc: float, m200: float, h: float = 0.7) -> float:
    """Circular velocity of NFW halo at radius r [kpc]."""
    if r_kpc <= 0 or m200 <= 0:
        return 0.0
    c = dutton_maccio_c200(m200, h=h)
    r200 = r200_from_m200(m200, h=h)
    rs = r200 / c
    x = r_kpc / rs
    gc = math.log(1.0 + c) - c / (1.0 + c)
    if gc <= 0:
        return 0.0
    # M(<r) = M200 * (ln(1+x)-x/(1+x)) / gc
    if x <= 0:
        return 0.0
    mx = math.log(1.0 + x) - x / (1.0 + x)
    m_enc = m200 * mx / gc
    return math.sqrt(G_KPC * m_enc / r_kpc)


def v_total_nfw_baryon(
    r: float,
    vgas: float,
    vdisk: float,
    vbul: float,
    m200: float,
    upsilon_disk: float,
    upsilon_bulge: float,
) -> float:
    vb2 = v_bar_squared(vgas, vdisk, vbul, upsilon_disk, upsilon_bulge)
    vn = v_nfw(r, m200)
    return math.sqrt(max(vb2, 0.0) + vn ** 2)
