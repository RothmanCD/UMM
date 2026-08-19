"""Global A fit and per-galaxy NFW mass fit."""
from __future__ import annotations

import math
from typing import Dict, List, Sequence, Tuple

from .parse import GalaxyCurve
from .physics import (
    N_G_PRIMARY,
    UPSILON_BULGE_FIDUCIAL,
    UPSILON_DISK_FIDUCIAL,
    a_from_v,
    a_ind,
    outer_mean_abs_frac_residual,
    r_out_threshold,
    sigma_bar,
    v_bar,
    v_model_from_a,
    v_total_nfw_baryon,
)


def galaxy_outer_mask(curve: GalaxyCurve) -> Tuple[float, List[int]]:
    r = curve.r
    r_out = r_out_threshold(min(r), max(r))
    idx = [i for i, ri in enumerate(r) if ri >= r_out]
    return r_out, idx


def model_v_umm(
    curve: GalaxyCurve,
    A: float,
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
    n_g: float = N_G_PRIMARY,
) -> List[float]:
    out: List[float] = []
    for p in curve.points:
        vb = v_bar(p.vgas, p.vdisk, p.vbul, upsilon_disk, upsilon_bulge)
        ab = a_from_v(vb, p.r_kpc)
        sig = sigma_bar(
            p.sb_disk, p.sb_bul, p.vgas, p.r_kpc, upsilon_disk, upsilon_bulge
        )
        ai = a_ind(sig, A, n_g=n_g)
        out.append(v_model_from_a(p.r_kpc, ab + ai))
    return out


def model_v_baryon(
    curve: GalaxyCurve,
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
) -> List[float]:
    return [
        v_bar(p.vgas, p.vdisk, p.vbul, upsilon_disk, upsilon_bulge) for p in curve.points
    ]


def chi2_outer_umm(
    galaxies: Dict[str, GalaxyCurve],
    selected: Sequence[str],
    A: float,
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
    n_g: float = N_G_PRIMARY,
) -> float:
    chi2 = 0.0
    for gid in selected:
        g = galaxies[gid]
        r_out, idx = galaxy_outer_mask(g)
        vmod = model_v_umm(g, A, upsilon_disk, upsilon_bulge, n_g=n_g)
        for i in idx:
            p = g.points[i]
            if p.v_obs <= 0:
                continue
            sig = p.e_vobs if p.e_vobs and p.e_vobs > 0 else 0.05 * p.v_obs
            if sig <= 0:
                sig = 1.0
            chi2 += ((p.v_obs - vmod[i]) / sig) ** 2
    return chi2


def fit_global_A(
    galaxies: Dict[str, GalaxyCurve],
    selected: Sequence[str],
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
    n_g: float = N_G_PRIMARY,
    A_grid: Sequence[float] | None = None,
) -> Tuple[float, float, float]:
    """Return (A_best, sigma_A, chi2_min).

    Coarse-then-fine grid search on outer-disk χ²; σ_A from Δχ²=1 interval
    on a dense log-A mesh about the minimum. A in m s^-2.
    """
    if A_grid is None:
        A_grid = [10 ** x for x in _frange(-13.0, -9.0, 0.05)]

    best_A = A_grid[0]
    best_chi = float("inf")
    for A in A_grid:
        c = chi2_outer_umm(
            galaxies, selected, A, upsilon_disk, upsilon_bulge, n_g=n_g
        )
        if c < best_chi:
            best_chi = c
            best_A = A

    # Dense log mesh about the minimum for Δχ² = 1
    log_c = math.log10(max(best_A, 1e-20))
    dense = [10 ** (log_c + d) for d in _frange(-0.8, 0.8, 0.01)]
    chi_by_A: List[Tuple[float, float]] = []
    for A in dense:
        c = chi2_outer_umm(
            galaxies, selected, A, upsilon_disk, upsilon_bulge, n_g=n_g
        )
        chi_by_A.append((A, c))
        if c < best_chi:
            best_chi = c
            best_A = A

    # Formal Δχ²=1 (often too tight when χ²≫N from underestimated σ_V)
    below = [A for A, c in chi_by_A if c <= best_chi + 1.0]
    if len(below) >= 2:
        sigma_formal = 0.5 * (max(below) - min(below))
    else:
        eps = 0.05 * best_A
        c_plus = chi2_outer_umm(
            galaxies, selected, best_A + eps, upsilon_disk, upsilon_bulge, n_g=n_g
        )
        d2 = (c_plus - best_chi) / (eps ** 2)
        sigma_formal = math.sqrt(1.0 / d2) if d2 > 0 else 0.2 * best_A

    # Subsample jackknife (≤15 leave-one-outs) for realistic σ_A when χ²≫N
    jack: List[float] = []
    sel_list = list(selected)
    jgrid = [10 ** x for x in _frange(-13.0, -9.0, 0.1)]
    step = max(1, len(sel_list) // 15)
    for k in range(0, len(sel_list), step):
        leave = sel_list[:k] + sel_list[k + 1 :]
        if not leave:
            continue
        ba, bc = jgrid[0], float("inf")
        for A in jgrid:
            c = chi2_outer_umm(
                galaxies, leave, A, upsilon_disk, upsilon_bulge, n_g=n_g
            )
            if c < bc:
                bc, ba = c, A
        jack.append(ba)
    if len(jack) >= 2:
        mean_j = sum(jack) / len(jack)
        sigma_jack = math.sqrt(
            sum((a - mean_j) ** 2 for a in jack) / max(len(jack) - 1, 1)
        )
    else:
        sigma_jack = 0.0

    sigma = max(sigma_formal, sigma_jack, 0.05 * best_A)
    if not math.isfinite(sigma) or sigma <= 0:
        sigma = 0.2 * best_A
    return best_A, sigma, best_chi


def _frange(a: float, b: float, step: float) -> List[float]:
    out = []
    x = a
    while x <= b + 1e-12:
        out.append(x)
        x += step
    return out


def fit_nfw_m200(
    curve: GalaxyCurve,
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
) -> float:
    """One-parameter NFW: minimize outer-disk χ² over M200."""
    r_out, idx = galaxy_outer_mask(curve)
    if not idx:
        return 1.0e11

    def chi2(m200: float) -> float:
        c2 = 0.0
        for i in idx:
            p = curve.points[i]
            if p.v_obs <= 0:
                continue
            vm = v_total_nfw_baryon(
                p.r_kpc,
                p.vgas,
                p.vdisk,
                p.vbul,
                m200,
                upsilon_disk,
                upsilon_bulge,
            )
            sig = p.e_vobs if p.e_vobs and p.e_vobs > 0 else 0.05 * p.v_obs
            if sig <= 0:
                sig = 1.0
            c2 += ((p.v_obs - vm) / sig) ** 2
        return c2

    # grid in log M200 from 1e9 to 1e14 (covers massive outer disks e.g. UGC11914)
    logs = [9.0 + i * 0.05 for i in range(101)]  # 9.00 .. 14.00
    best_m = 1e11
    best_c = float("inf")
    for lg in logs:
        m = 10 ** lg
        c = chi2(m)
        if c < best_c:
            best_c = c
            best_m = m
    # Refine relative to frozen grid-best only (do not cascade best_m *= fac)
    grid_best = best_m
    for fac in [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5]:
        m = grid_best * fac
        c = chi2(m)
        if c < best_c:
            best_c = c
            best_m = m
    return best_m


def epsilon_for_models(
    curve: GalaxyCurve,
    v_model: Sequence[float],
) -> Tuple[float, int, float]:
    r = curve.r
    r_out = r_out_threshold(min(r), max(r))
    v_obs = [p.v_obs for p in curve.points]
    eps, n_out = outer_mean_abs_frac_residual(r, v_obs, list(v_model), r_out)
    return eps, n_out, r_out


def mean_outer_sigma_bar(
    curve: GalaxyCurve,
    upsilon_disk: float = UPSILON_DISK_FIDUCIAL,
    upsilon_bulge: float = UPSILON_BULGE_FIDUCIAL,
) -> float:
    r = curve.r
    r_out = r_out_threshold(min(r), max(r))
    vals = []
    for p in curve.points:
        if p.r_kpc < r_out:
            continue
        vals.append(
            sigma_bar(
                p.sb_disk, p.sb_bul, p.vgas, p.r_kpc, upsilon_disk, upsilon_bulge
            )
        )
    if not vals:
        return float("nan")
    return sum(vals) / len(vals)
