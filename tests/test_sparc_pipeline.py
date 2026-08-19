#!/usr/bin/env python3
"""Tests driving the shipped SPARC residual pipeline (not a reimplementation)."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from analysis.sparc_residuals.physics import (  # noqa: E402
    a_ind,
    outer_mean_abs_frac_residual,
    r_out_threshold,
    v_bar,
)
from analysis.sparc_residuals.parse import parse_mass_models, parse_table1  # noqa: E402
from analysis.sparc_residuals.select import apply_mandatory_cuts, selected_ids  # noqa: E402
from analysis.sparc_residuals.parse import build_galaxies  # noqa: E402

RESULTS = ROOT / "analysis" / "sparc_residuals" / "results"
SPARC = ROOT / "Data" / "SPARC"


class TestNFWFit(unittest.TestCase):
    def test_refine_uses_frozen_grid_best(self):
        """Regression: refine factors must multiply grid-best, not cascading best_m."""
        import inspect
        from analysis.sparc_residuals import fit as fit_mod

        src = inspect.getsource(fit_mod.fit_nfw_m200)
        self.assertIn("grid_best", src)
        self.assertIn("grid_best * fac", src)
        # must not assign best_m = best_m * fac pattern in refine
        self.assertNotIn("m = best_m * fac", src)

    def test_m200_grid_reaches_1e14(self):
        import inspect
        from analysis.sparc_residuals import fit as fit_mod

        src = inspect.getsource(fit_mod.fit_nfw_m200)
        self.assertIn("1e14", src)
        self.assertIn("range(101)", src)

    def test_ugc11914_nfw_prefers_high_mass_if_present(self):
        """If UGC11914 is in the mass models, best M200 should not be stuck at grid edge artifact."""
        from analysis.sparc_residuals.fit import fit_nfw_m200
        from analysis.sparc_residuals.parse import (
            build_galaxies,
            parse_mass_models,
            parse_table1,
        )

        t1 = parse_table1(SPARC / "SPARC_Lelli2016c.mrt")
        mm = parse_mass_models(SPARC / "MassModels_Lelli2016c.mrt")
        gals = build_galaxies(t1, mm)
        if "UGC11914" not in gals:
            self.skipTest("UGC11914 not in catalog")
        m200 = fit_nfw_m200(gals["UGC11914"])
        # Skeptic: true min near ~4e13; must not return the broken cascade ~2.57e13
        self.assertGreater(m200, 3.0e13)
        self.assertLess(m200, 1.0e14)


class TestPhysicsHelpers(unittest.TestCase):
    def test_epsilon_definition(self):
        r = [1.0, 2.0, 3.0, 4.0, 5.0]
        vobs = [100.0, 100.0, 100.0, 100.0, 100.0]
        vmod = [100.0, 90.0, 80.0, 100.0, 50.0]
        r_out = r_out_threshold(1.0, 5.0)  # 1 + 0.7*4 = 3.8
        eps, n = outer_mean_abs_frac_residual(r, vobs, vmod, r_out)
        # points R>=3.8: R=4 (0), R=5 (0.5) → mean 0.25
        self.assertEqual(n, 2)
        self.assertAlmostEqual(eps, 0.25, places=10)

    def test_a_ind_scaling(self):
        a1 = a_ind(1.0, 1.0e-10, n_g=0.4)
        a16 = a_ind(16.0, 1.0e-10, n_g=0.4)
        # 16^0.4 = (2^4)^0.4 = 2^1.6 ≈ 3.0314
        self.assertAlmostEqual(a16 / a1, 16 ** 0.4, places=8)

    def test_v_bar_ml_half(self):
        # Υ=1 → V=Vdisk; Υ=0.5 → V=sqrt(0.5)*Vdisk for pure disk
        self.assertAlmostEqual(v_bar(0, 100, 0, 1.0, 0.7), 100.0, places=6)
        self.assertAlmostEqual(v_bar(0, 100, 0, 0.5, 0.7), math.sqrt(0.5) * 100, places=6)

    def test_r_out_formula(self):
        self.assertAlmostEqual(r_out_threshold(1.0, 11.0), 1.0 + 0.7 * 10.0, places=10)


class TestSelectionAndOutputs(unittest.TestCase):
    def test_inputs_exist_and_hashable(self):
        for name in [
            "SPARC_Lelli2016c.mrt",
            "MassModels_Lelli2016c.mrt",
        ]:
            p = SPARC / name
            self.assertTrue(p.is_file(), p)
            h = hashlib.sha256(p.read_bytes()).hexdigest()
            self.assertEqual(len(h), 64)

    def test_selection_q1_incl(self):
        t1 = parse_table1(SPARC / "SPARC_Lelli2016c.mrt")
        mm = parse_mass_models(SPARC / "MassModels_Lelli2016c.mrt")
        gals = build_galaxies(t1, mm)
        cuts = apply_mandatory_cuts(gals, r_last_over_rdisc_min=2.5)
        sel = selected_ids(cuts)
        self.assertGreaterEqual(len(sel), 20)
        for gid in sel:
            m = gals[gid].meta
            self.assertIsNotNone(m)
            self.assertEqual(m.Q, 1)
            self.assertGreaterEqual(m.incl_deg, 30.0)
            self.assertGreaterEqual(gals[gid].n_r, 8)

    def test_results_tables_exist_and_schema(self):
        for name in [
            "selection_manifest.json",
            "residuals_main.csv",
            "residuals_main.json",
            "residuals_summary.json",
            "production_report.md",
            "status_for_v14.md",
        ]:
            self.assertTrue((RESULTS / name).is_file(), name)

        summary = json.loads((RESULTS / "residuals_summary.json").read_text())
        self.assertIn("N_gal", summary)
        self.assertIn("A_global", summary)
        self.assertIn("med_epsilon_umm", summary)
        self.assertIn("selection_id", summary)

        with open(RESULTS / "residuals_main.csv") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), summary["N_gal"])
        required = [
            "galaxy_id",
            "morph_type",
            "epsilon_umm",
            "epsilon_baryon",
            "epsilon_nfw",
            "A_global",
        ]
        for col in required:
            self.assertIn(col, rows[0])

        # single global A
        As = {r["A_global"] for r in rows}
        self.assertEqual(len(As), 1)

        # medians recompute from main table
        def med(key):
            vals = sorted(float(r[key]) for r in rows)
            n = len(vals)
            if n % 2:
                return vals[n // 2]
            return 0.5 * (vals[n // 2 - 1] + vals[n // 2])

        self.assertAlmostEqual(med("epsilon_umm"), summary["med_epsilon_umm"], places=8)
        self.assertAlmostEqual(
            med("epsilon_baryon"), summary["med_epsilon_baryon"], places=8
        )

        # manifest hashes match files
        man = json.loads((RESULTS / "selection_manifest.json").read_text())
        for key, info in man["input_files"].items():
            path = ROOT / info["path"]
            self.assertTrue(path.is_file(), path)
            h = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(h, info["sha256"], key)

        # no post-hoc: frozen list length equals N
        self.assertEqual(len(man["frozen_galaxy_list"]), summary["N_gal"])

        # report claims measured
        report = (RESULTS / "production_report.md").read_text()
        self.assertIn("measured catalog residuals", report.lower())

    def test_v13_untouched(self):
        freeze = (ROOT / "notes" / "UMM_v13_Freeze.md").read_text()
        # hash line from freeze
        self.assertIn("da045cc5d68fe95141a54366501acc2049435b2f75517021d5d8d96ad7788093", freeze)
        h = hashlib.sha256((ROOT / "paper" / "UMM_Paper_Draft_v13.tex").read_bytes()).hexdigest()
        self.assertEqual(h, "da045cc5d68fe95141a54366501acc2049435b2f75517021d5d8d96ad7788093")


if __name__ == "__main__":
    unittest.main()
