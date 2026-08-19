"""Unit tests for UDG confrontation helpers — drive real shipped functions."""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.sparc_residuals.physics import a_ind as sparc_a_ind  # noqa: E402
from analysis.udg_confrontation.src.physics import (  # noqa: E402
    A_FROZEN,
    A_SI_TO_ASTRO,
    G_KPC,
    N_G,
    SIGMA_REF,
    a_bar_half_light,
    a_ind_logistic,
    a_ind_powerlaw,
    mean_sigma_within_re,
    predict_sigma_three_models,
    sigma_from_a_tot,
    wolf_mass_from_sigma,
    wolf_sigma_from_mass,
)


class TestWolfEstimator(unittest.TestCase):
    def test_roundtrip(self):
        m = 1.0e8
        r = 2.2
        sig = wolf_sigma_from_mass(m, r)
        m2 = wolf_mass_from_sigma(sig, r)
        self.assertAlmostEqual(m2 / m, 1.0, places=10)

    def test_danieli_df2_scale(self):
        """Wolf mass at Danieli-like numbers is O(10^8 M_sun)."""
        # σ=8.5, R_e=2.2 → M ~ 1.5e8
        m = wolf_mass_from_sigma(8.5, 2.2)
        self.assertGreater(m, 1.0e8)
        self.assertLess(m, 2.0e8)

    def test_sigma_from_a_matches_wolf(self):
        m_half = 1.0e8
        r = 2.2
        a = G_KPC * m_half / (r**2)
        sig_a = sigma_from_a_tot(a, r)
        sig_w = wolf_sigma_from_mass(m_half, r)
        self.assertAlmostEqual(sig_a, sig_w, places=10)


class TestAIndMatchesSPARC(unittest.TestCase):
    def test_unit_conversion_matches_sparc_physics(self):
        """Power-law a_ind must match analysis/sparc_residuals/physics.py."""
        for sigma in (0.5, 1.0, 5.0, 20.0):
            ours = a_ind_powerlaw(sigma, A=A_FROZEN, n_g=N_G)
            theirs = sparc_a_ind(sigma, A=A_FROZEN, n_g=N_G)
            self.assertAlmostEqual(ours, theirs, places=8)

    def test_ref_amplitude(self):
        a = a_ind_powerlaw(SIGMA_REF, A=A_FROZEN)
        self.assertAlmostEqual(a, A_FROZEN * A_SI_TO_ASTRO, places=6)


class TestLogistic(unittest.TestCase):
    def test_midpoint_matches_powerlaw_at_sigma_ref(self):
        a_pl = a_ind_powerlaw(SIGMA_REF)
        a_log = a_ind_logistic(SIGMA_REF, sigma_mid=SIGMA_REF)
        self.assertAlmostEqual(a_log, a_pl, places=8)

    def test_low_density_suppression(self):
        """At very low Σ, logistic should fall below pure power-law (α/2 > n_g)."""
        sigma = 0.05
        a_pl = a_ind_powerlaw(sigma)
        a_log = a_ind_logistic(sigma, sigma_mid=SIGMA_REF)
        self.assertLess(a_log, a_pl)

    def test_no_nan(self):
        for s in (0.01, 0.1, 1.0, 10.0, 100.0):
            self.assertTrue(math.isfinite(a_ind_logistic(s)))


class TestMeanSigmaAndPredictions(unittest.TestCase):
    def test_mean_sigma_df2(self):
        # M=2e8, R_e=2.2 kpc → Σ ~ few M_sun/pc^2
        s = mean_sigma_within_re(2.0e8, 2.2)
        self.assertGreater(s, 1.0)
        self.assertLess(s, 20.0)

    def test_three_models_ordering_udg(self):
        """For UDG-like inputs, a_ind > 0 so σ_PL > σ_stars."""
        p = predict_sigma_three_models(2.0e8, 2.2)
        self.assertGreater(p["sigma_powerlaw"], p["sigma_stars"])
        self.assertGreater(p["a_ind_powerlaw"], 0.0)

    def test_a_bar_positive(self):
        self.assertGreater(a_bar_half_light(2e8, 2.2), 0.0)

    def test_no_refit_constants(self):
        """Frozen A and n_g must remain the SPARC baseline values."""
        self.assertAlmostEqual(A_FROZEN, 1.41e-11)
        self.assertAlmostEqual(N_G, 0.4)


class TestEndToEndDF2Recompute(unittest.TestCase):
    """Recompute DF2 from adopted inputs; assert agreement with runner formulas."""

    def test_df2_recompute(self):
        m_star = 2.0e8
        r_e = 2.2
        p = predict_sigma_three_models(m_star, r_e)
        # Stars alone ~7 km/s
        self.assertGreater(p["sigma_stars"], 5.0)
        self.assertLess(p["sigma_stars"], 9.0)
        # Power-law substantially higher
        self.assertGreater(p["sigma_powerlaw"], p["sigma_stars"] + 3.0)
        # Explicit manual recompute
        sig_bar = mean_sigma_within_re(m_star, r_e)
        a_bar = a_bar_half_light(m_star, r_e)
        a_pl = a_ind_powerlaw(sig_bar)
        sig_pl = math.sqrt((a_bar + a_pl) * r_e / 4.0)
        self.assertAlmostEqual(sig_pl, p["sigma_powerlaw"], places=10)


class TestFCC224PublishedUncertainties(unittest.TestCase):
    """Buzzo et al. 2025 published σ_stars=7.82^{+6.74}_{-4.36} must be in manifest."""

    def test_asymmetric_not_fabricated_symmetric(self):
        from analysis.udg_confrontation.src.manifest_data import SYSTEMS, sigma_obs_1sigma

        fcc = next(s for s in SYSTEMS if s["short"] == "FCC224")
        self.assertAlmostEqual(fcc["sigma_obs"]["value"], 7.82)
        self.assertAlmostEqual(fcc["sigma_obs"]["unc_plus"], 6.74)
        self.assertAlmostEqual(fcc["sigma_obs"]["unc_minus"], 4.36)
        # 1σ proxy for n_sigma uses max of published legs (not a fake ±5.5)
        self.assertAlmostEqual(sigma_obs_1sigma(fcc), 6.74)
        self.assertNotAlmostEqual(fcc["sigma_obs"]["unc_plus"], fcc["sigma_obs"]["unc_minus"])


if __name__ == "__main__":
    unittest.main()
