"""Selection cuts per UMM_SPARC_Residual_Table_Design.md §3."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .parse import GalaxyCurve


@dataclass
class CutResult:
    galaxy_id: str
    passed: bool
    reasons: List[str] = field(default_factory=list)
    Q: Optional[int] = None
    incl_deg: Optional[float] = None
    R_disc: Optional[float] = None
    R_last: Optional[float] = None
    N_R: Optional[int] = None
    coverage_ratio: Optional[float] = None


def apply_mandatory_cuts(
    galaxies: Dict[str, GalaxyCurve],
    r_last_over_rdisc_min: float = 2.5,
) -> Dict[str, CutResult]:
    """§3.1 mandatory cuts. Pre-registered; no residual-based drops."""
    results: Dict[str, CutResult] = {}
    for gid, g in sorted(galaxies.items()):
        reasons: List[str] = []
        meta = g.meta
        if meta is None:
            results[gid] = CutResult(gid, False, ["missing_table1_metadata"])
            continue
        Q = meta.Q
        inc = meta.incl_deg
        rdisc = meta.Rdisk_kpc
        n_r = g.n_r
        r_last = g.r_last if n_r else float("nan")
        cov = (r_last / rdisc) if (rdisc and rdisc > 0 and n_r) else float("nan")

        if Q != 1:
            reasons.append(f"Q={Q}!=1")
        if inc is None or inc < 30.0:
            reasons.append(f"incl={inc}<30")
        if n_r < 8:
            reasons.append(f"N_R={n_r}<8")
        if not (meta.D_Mpc and meta.D_Mpc > 0 and math_isfinite(meta.D_Mpc)):
            reasons.append("bad_distance")
        if not (inc is not None and math_isfinite(inc)):
            reasons.append("bad_inclination")
        if not (rdisc and rdisc > 0 and math_isfinite(rdisc)):
            reasons.append("bad_Rdisk")
        elif not (math_isfinite(cov) and cov >= r_last_over_rdisc_min):
            reasons.append(f"coverage={cov:.3f}<{r_last_over_rdisc_min}")

        results[gid] = CutResult(
            galaxy_id=gid,
            passed=(len(reasons) == 0),
            reasons=reasons,
            Q=Q,
            incl_deg=inc,
            R_disc=rdisc,
            R_last=r_last,
            N_R=n_r,
            coverage_ratio=cov if math_isfinite(cov) else None,
        )
    return results


def math_isfinite(x: float) -> bool:
    return x == x and x not in (float("inf"), float("-inf"))


def selected_ids(cut_results: Dict[str, CutResult]) -> List[str]:
    return [gid for gid, c in sorted(cut_results.items()) if c.passed]


def cut_tallies(cut_results: Dict[str, CutResult]) -> Dict[str, int]:
    tallies = {
        "n_total": len(cut_results),
        "n_pass": sum(1 for c in cut_results.values() if c.passed),
        "fail_Q": 0,
        "fail_incl": 0,
        "fail_NR": 0,
        "fail_coverage": 0,
        "fail_meta": 0,
    }
    for c in cut_results.values():
        if c.passed:
            continue
        rs = " ".join(c.reasons)
        if "Q=" in rs:
            tallies["fail_Q"] += 1
        if "incl=" in rs:
            tallies["fail_incl"] += 1
        if "N_R=" in rs:
            tallies["fail_NR"] += 1
        if "coverage=" in rs:
            tallies["fail_coverage"] += 1
        if "missing" in rs or "bad_" in rs:
            tallies["fail_meta"] += 1
    return tallies
