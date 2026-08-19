"""Parse SPARC Table1 (galaxy sample) and MassModels radial tables."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


HUBBLE_TYPE = {
    0: "S0",
    1: "Sa",
    2: "Sab",
    3: "Sb",
    4: "Sbc",
    5: "Sc",
    6: "Scd",
    7: "Sd",
    8: "Sdm",
    9: "Sm",
    10: "Im",
    11: "BCD",
}


@dataclass
class GalaxyMeta:
    galaxy_id: str
    T: int
    morph_type: str
    D_Mpc: float
    e_D: float
    f_D: int
    incl_deg: float
    e_inc: float
    L36: float
    Rdisk_kpc: float
    MHI: float
    RHI_kpc: float
    Vflat: float
    Q: int
    ref: str


@dataclass
class RadialPoint:
    r_kpc: float
    v_obs: float
    e_vobs: float
    vgas: float
    vdisk: float
    vbul: float
    sb_disk: float
    sb_bul: float


@dataclass
class GalaxyCurve:
    galaxy_id: str
    meta: Optional[GalaxyMeta]
    points: List[RadialPoint] = field(default_factory=list)

    @property
    def r(self) -> List[float]:
        return [p.r_kpc for p in self.points]

    @property
    def n_r(self) -> int:
        return len(self.points)

    @property
    def r_min(self) -> float:
        return min(p.r_kpc for p in self.points)

    @property
    def r_last(self) -> float:
        return max(p.r_kpc for p in self.points)


def parse_table1(path: Path) -> Dict[str, GalaxyMeta]:
    """Parse SPARC_Lelli2016c.mrt (Table 1 galaxy sample)."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    out: Dict[str, GalaxyMeta] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith(
            ("Title", "Authors", "Table", "=", "Byte", "-", "Note", " 0 =", " 1 =")
        ):
            continue
        # Data lines: free-format whitespace after optional leading spaces
        parts = line.split()
        if len(parts) < 18:
            continue
        # Galaxy name may be first token; must not look like a pure note key
        name = parts[0]
        if not any(ch.isalpha() for ch in name):
            continue
        try:
            T = int(parts[1])
            D = float(parts[2])
            e_D = float(parts[3])
            f_D = int(parts[4])
            Inc = float(parts[5])
            e_Inc = float(parts[6])
            L36 = float(parts[7])
            # e_L, Reff, SBeff, Rdisk, SBdisk, MHI, RHI, Vflat, e_Vflat, Q
            Rdisk = float(parts[11])
            MHI = float(parts[13])
            RHI = float(parts[14])
            Vflat = float(parts[15])
            Q = int(parts[17])
            ref = parts[18] if len(parts) > 18 else ""
        except (ValueError, IndexError):
            continue
        if Q not in (1, 2, 3):
            continue
        out[name] = GalaxyMeta(
            galaxy_id=name,
            T=T,
            morph_type=HUBBLE_TYPE.get(T, str(T)),
            D_Mpc=D,
            e_D=e_D,
            f_D=f_D,
            incl_deg=Inc,
            e_inc=e_Inc,
            L36=L36,
            Rdisk_kpc=Rdisk,
            MHI=MHI,
            RHI_kpc=RHI,
            Vflat=Vflat,
            Q=Q,
            ref=ref,
        )
    return out


def parse_mass_models(path: Path) -> Dict[str, List[RadialPoint]]:
    """Parse MassModels_Lelli2016c.mrt radial mass models."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    curves: Dict[str, List[RadialPoint]] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith(
            ("Title", "Authors", "Author", "Table", "=", "Byte", "-", "Note")
        ):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        name = parts[0]
        try:
            # ID D R Vobs eV Vgas Vdisk Vbul SBdisk SBbul
            _D = float(parts[1])
            R = float(parts[2])
            Vobs = float(parts[3])
            eV = float(parts[4])
            Vgas = float(parts[5])
            Vdisk = float(parts[6])
            Vbul = float(parts[7])
            SBdisk = float(parts[8])
            SBbul = float(parts[9])
        except ValueError:
            continue
        curves.setdefault(name, []).append(
            RadialPoint(
                r_kpc=R,
                v_obs=Vobs,
                e_vobs=eV,
                vgas=Vgas,
                vdisk=Vdisk,
                vbul=Vbul,
                sb_disk=SBdisk,
                sb_bul=SBbul,
            )
        )
    # sort by radius
    for name in curves:
        curves[name].sort(key=lambda p: p.r_kpc)
    return curves


def build_galaxies(
    table1: Dict[str, GalaxyMeta], mass_models: Dict[str, List[RadialPoint]]
) -> Dict[str, GalaxyCurve]:
    out: Dict[str, GalaxyCurve] = {}
    for gid, pts in mass_models.items():
        meta = table1.get(gid)
        out[gid] = GalaxyCurve(galaxy_id=gid, meta=meta, points=pts)
    return out
