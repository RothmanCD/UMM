# First quantitative UMM confrontation with DM-deficient UDG kinematics

**Status:** Illustrative first estimate · Small sample (N = 4) · Systematics remain significant

**Frozen inputs (no re-fit):**
- \(A = (1.41 \pm 0.07)\times 10^{-11}\,\mathrm{m\,s^{-2}}\) (SPARC residual amplitude; selection_id `ea6c92761ffd3c5a`, \(\Upsilon_*=0.5\), \(n_g=0.4\))
- \(\Sigma_{\rm ref} = 1\,M_\odot\,\mathrm{pc^{-2}}\)
- Packing \(\alpha_{\rm eff} \approx 1.8\) (preferred band \([1.66, 1.84]\)); logistic midpoint window \(\Sigma_{\rm mid} \in \{0.3, 1.0, 3.0\}\,M_\odot\,\mathrm{pc^{-2}}\) (joint-normalization placement about SPARC-relevant intermediate densities; no new free parameters)
- All other UMM parameters remain those of the frozen v14 baseline.

---

## 1. Sample and data table

Primary systems: NGC 1052-DF2, DF4, DF9, and FCC 224. Distances prefer SBF/TRGB or homogeneous group/cluster scales. Only published peer-reviewed or arXiv-vetted kinematics are used.

**Surface-density limitation (all systems):** No published high-resolution \(\Sigma(R)\) mass profile is adopted. Mean baryonic surface density within \(R_e\) is computed as \(\langle\Sigma\rangle_e = M_*/(2\pi R_e^2)\) (half-mass within \(R_e\), mass-follows-light, gas neglected for these quiescent systems). Local \(\Sigma(R_e)\) for a Sérsic profile would be lower by an order-unity factor; we do not invent a radial profile.

| System | \(D\) [Mpc] | \(M_*\) [\(10^8 M_\odot\)] | \(R_e\) [kpc] | \(\langle\Sigma\rangle_e\) [\(M_\odot\mathrm{pc^{-2}}\)] | \(\sigma_{\rm obs}\) [km/s] | Tracer |
|--------|------------|---------------------------|---------------|------------------------------------------------------|------------------------------|--------|
| NGC1052-DF2 | 20.0±1.7 | 2.00±0.40 | 2.20 | 6.58 | 8.50^{+2.30}_{-3.10} | stars (KCWI integrated light) |
| NGC1052-DF4 | 20.0±1.6 | 1.50±0.40 | 1.60 | 9.33 | 8.00^{+2.30}_{-1.90} | stars (KCWI integrated light) |
| NGC1052-DF9 | 20.6±2.0 | 1.40±0.30 | 1.10 | 18.41 | 6.40^{+4.00}_{-4.30} | stars (KCWI; broadening- and binary-corrected) |
| FCC224 | 20.0±1.5 | 1.74±0.16 | 1.89 | 7.75 | 7.82^{+6.74}_{-4.36} | stars (KCWI red-arm CaT; broadening-corrected) |

### Citations (adopted values)

**NGC1052-DF2**
- Distance: van Dokkum et al. 2018, Nature, 555, 629 (arXiv:1803.10237)
- \(M_*\): van Dokkum et al. 2018, Nature; M/L_V≈2.0
- \(R_e\): van Dokkum et al. 2018, Nature (R_e=22.6'' → 2.2 kpc at 20 Mpc)
- \(\sigma\): Danieli et al. 2019, ApJL, 874, L12 (arXiv:1901.03711)
- Profile note: No published high-resolution mass surface-density profile used. Mean Σ within R_e computed from M_* and R_e (mass-follows-light).

**NGC1052-DF4**
- Distance: van Dokkum et al. 2019, ApJL, 874, L5 (SBF 19.9±2.8); Danieli et al. 2020 TRGB ≈20.0
- \(M_*\): van Dokkum et al. 2019, ApJL (M/L_V=2.0±0.5; L_V=(7.7±0.8)×10^7 L_⊙)
- \(R_e\): van Dokkum et al. 2019, ApJL (R_e=1.6 kpc; n=0.79)
- \(\sigma\): Shen et al. 2023, ApJ, 958, 3 (arXiv:2309.08592)
- Profile note: No published radial Σ profile used; mean Σ within R_e from M_* and R_e.

**NGC1052-DF9**
- Distance: Keim et al. 2026, ApJ (arXiv:2603.15860); trail scale Keim et al. 2025
- \(M_*\): Keim et al. 2026 (from Gannon et al. 2023 photometry; M/L_g=2.2)
- \(R_e\): Keim et al. 2026; Gannon et al. 2023 (11.1'' at 20.6 Mpc → 1.1 kpc)
- \(\sigma\): Keim et al. 2026, ApJ (arXiv:2603.15860)
- Profile note: No published Σ profile; mean Σ within R_e from M_* and R_e.

**FCC224**
- Distance: Buzzo et al. 2025, A&A (arXiv:2502.05405); Tang et al. 2025a photometry at 20 Mpc
- \(M_*\): Buzzo et al. 2025; Tang et al. 2025a: log(M_*/M_⊙)=8.24±0.04
- \(R_e\): Buzzo et al. 2025; Tang et al. 2025a: R_e=1.89±0.01 kpc
- \(\sigma\): Buzzo et al. 2025, A&A (arXiv:2502.05405)
- Profile note: No published full Σ profile used; mean Σ within R_e from M_* and R_e. KCWI FoV covers only ~0.4 R_e,circ for population gradients.

### Future targets (not in residual sample)

- **NGC1052 trail members (excl. DF2/DF4/DF9):** Keim et al. 2025/2026: remaining trail galaxies up to ~100× fainter; radial velocities partially available but internal kinematics not yet of KCWI/MUSE quality comparable to DF2/DF4/DF9. (van Dokkum et al. 2022a; Keim et al. 2025, 2026)
- **FCC240:** Buzzo et al. 2026 (arXiv:2605.24099): close companion to FCC 224; analogue discussion; published kinematics of residual-sample quality not adopted here as primary. (Buzzo et al. 2026, ApJ (arXiv:2605.24099))
- **Tidal-dwarf candidates (e.g. NGC 5291N, etc.):** Classic TDGs can be DM-poor but generally have rotation-supported HI kinematics and different formation channels; not included without a homogeneous pressure-supported σ of UDG quality. (Lelli et al. 2015; Bournaud et al. 2007 (context only))
- **Almost-dark HI clouds / HI-bearing UDGs:** Interesting for low-Σ shoulder of Φ(ξ) but require published stellar or GC dispersions of comparable quality; list as future targets only. (—)

---

## 2. Estimator and model definitions

### Wolf et al. (2010) mass estimator

\[
M_{1/2} = \frac{4\,\sigma_{\rm los}^{2}\,R_e}{G},
\]
with the inverse \(\sigma = \sqrt{G M_{1/2}/(4 R_e)}\).

**Assumptions (stated explicitly):**
1. Pressure-supported kinematics (dispersion-dominated).
2. Mass follows light for the stellar comparison; half the stellar mass lies within \(R_e\).
3. Wolf estimator near the deprojected half-light radius is only weakly sensitive to orbital anisotropy.
4. Flat or slowly varying \(\sigma(R)\) (no full Jeans solution with free \(\beta(r)\)).
5. Negligible cold gas in \(\Sigma_{\rm bar}\).
6. Thin-disk UMM specialization \(a_{\rm ind}(R)=A(\Sigma/\Sigma_{\rm ref})^{n_g}\) is applied using mean \(\Sigma\) within \(R_e\) as a proxy — these systems are spheroidal; this is a limitation of the first estimate.

### Three models (identical estimator and inputs)

(a) **Stars alone:** \(a_{\rm tot}=a_{\rm bar}=G(M_*/2)/R_e^2\).

(b) **Stars + intermediate power-law:** \(a_{\rm ind}=A(\langle\Sigma\rangle_e/\Sigma_{\rm ref})^{0.4}\), \(a_{\rm tot}=a_{\rm bar}+a_{\rm ind}\).

(c) **Stars + logistic \(\Phi(\xi)\):** packing \(\alpha_{\rm eff}=1.8\), midpoint at \(\Sigma_{\rm mid}\) in the joint-normalization window about \(\Sigma_{\rm ref}\), amplitude matched so \(a_{\rm ind}(\Sigma_{\rm mid})=A(\Sigma_{\rm mid}/\Sigma_{\rm ref})^{n_g}\). Central column uses \(\Sigma_{\rm mid}=\Sigma_{\rm ref}\); band reports min/max over the window.

\[
\sigma_{\rm model}=\sqrt{a_{\rm tot}\,R_e/4}.
\]

Fractional residual (SPARC spirit, signed): \(\varepsilon_\sigma=(\sigma_{\rm model}-\sigma_{\rm obs})/\sigma_{\rm obs}\).

---

## 3. Predicted vs observed dispersions

| System | \(\sigma_{\rm obs}\) | (a) stars | (b) + power-law | (c) + logistic | \(\varepsilon_{\rm PL}\) | \(n_\sigma^{\rm PL}\) |
|--------|---------------------|-----------|-----------------|----------------|----------------------|----------------------|
| DF2 | 8.50^{+2.30}_{-3.10} | 7.0 | 23.6 | 21.3 [18.1–23.4] | +1.78 | 4.80 |
| DF4 | 8.00^{+2.30}_{-1.90} | 7.1 | 21.8 | 18.9 [16.0–21.2] | +1.73 | 5.87 |
| DF9 | 6.40^{+4.00}_{-4.30} | 8.3 | 21.3 | 17.1 [14.6–19.5] | +2.32 | 3.44 |
| FCC224 | 7.82^{+6.74}_{-4.36} | 7.0 | 22.7 | 20.1 [17.0–22.3] | +1.90 | 2.20 |

| System | \(M_{1/2}^{\rm obs}\) [\(10^8 M_\odot\)] | \(M_*/2\) | \(M_{1/2}^{\rm PL}\) | \(\langle\Sigma\rangle_e\) | \(a_{\rm bar}\) | \(a_{\rm ind}^{\rm PL}\) |
|--------|------------------------------------------|----------|---------------------|--------------------------|-----------------|--------------------------|
| DF2 | 1.48±1.08 | 1.00 | 11.40 | 6.58 | 88.9 | 924.2 |
| DF4 | 0.95±0.55 | 0.75 | 7.08 | 9.33 | 126.0 | 1062.8 |
| DF9 | 0.42±0.56 | 0.70 | 4.63 | 18.41 | 248.8 | 1395.2 |
| FCC224 | 1.07±1.85 | 0.87 | 9.07 | 7.75 | 104.8 | 987.1 |

Accelerations in \((\mathrm{km\,s^{-1}})^2\,\mathrm{kpc^{-1}}\).

---

## 4. Residual summary

- Sample size: **N = 4** (not a large statistical population).
- Median \(\varepsilon_\sigma\) (stars alone): **-0.106**
- Median \(\varepsilon_\sigma\) (UMM power-law): **+1.841**
- Median \(\varepsilon_\sigma\) (UMM logistic, central \(\Sigma_{\rm mid}\)): **+1.538**
- Fraction within 1σ / 2σ (stars): 4/4 / 4/4
- Fraction within 1σ / 2σ (power-law UMM): 0/4 / 0/4
- Fraction within 1σ / 2σ (logistic UMM): 0/4 / 1/4

### Classification: **RULED OUT**

For every system in this small sample, the frozen intermediate-regime UMM power-law predicts a line-of-sight dispersion more than 2σ above the observed value (model systematically high). The pure-baryon column is consistent with the data within uncertainties. This is an illustrative first estimate on a small sample; systematics remain significant, but the direction of the residual is uniform.

*Illustrative first estimate; small sample (N=4); systematics remain significant (distance, M/L, residual rotation, possible tides, mean-Σ proxy for local Σ_bar, thin-disk specialization applied to pressure-supported systems).*

### Systematics

Dominant systematics and their effect on the power-law residual:

1. **Distance:** Scaling \(M_*\propto D^2\), \(R_e\propto D\) leaves \(\langle\Sigma\rangle\) invariant; characteristic accelerations are distance-invariant at leading order, while Wolf masses scale ∝ \(D\). Published analyses (Keim et al. 2026 App. B) show DM-deficiency conclusions are robust to the trail distance. Shifting each system by its \(\pm 1\sigma\) distance does not bring power-law \(\sigma\) into agreement with \(\sigma_{\rm obs}\).
2. **Stellar M/L:** \(\pm 1\sigma\) in \(M_*\) moves \(\sigma_{\rm stars}\) by ~half the fractional mass error and moves \(a_{\rm ind}\) weakly through \(\Sigma\). Residuals remain positive and large for the power-law column.
3. **Residual rotation / tides:** FCC 224 shows slow prolate rotation; DF2/DF4 KCWI maps show little ordered motion. An illustrative correction that attributes half of \(\sigma_{\rm obs}\) to rotation *increases* tension for UMM (true dispersion lower). A 30% tidal inflation of observed motions reduces tension only partially and does not reverse the sign of \(\varepsilon_{\rm PL}\) for the ensemble.
4. **Mean vs local \(\Sigma\):** Using a factor-of-two lower local \(\Sigma(R_e)\) still leaves \(a_{\rm ind}\gtrsim a_{\rm bar}\) and \(\sigma_{\rm PL}\) well above \(\sigma_{\rm obs}\).
5. **Geometry:** The galactic specialization assumes a thin disk; these UDGs are pressure-supported spheroids. This is a structural limitation of the first contact, not a free parameter.

Example (DF2) one-at-a-time shifts of \(\sigma_{\rm PL}\):

- Fiducial sigma_PL = 23.6 km/s (obs 8.5)
- D+1sigma: sigma_PL = 24.6, eps = +1.89
- D-1sigma: sigma_PL = 22.6, eps = +1.66
- M_star+1sigma: sigma_PL = 24.6
- M_star-1sigma: sigma_PL = 22.5

---

## 5. Limitations and sample size

- **N = 4** primary systems with published high-resolution stellar or GC kinematics of comparable quality. This is not a population study.
- Kinematic uncertainties are large and often asymmetric; broadening and binary corrections matter at \(\sigma\sim 5\text{--}10\) km/s.
- Mean \(\Sigma\) proxy and thin-disk formula applied to spheroidal systems.
- Logistic \(\xi_c\) is represented by a \(\Sigma_{\rm mid}\) window about \(\Sigma_{\rm ref}\) consistent with intermediate-regime SPARC placement; the paper does not publish a single numerical \(\xi_c\).
- No re-tuning of \(A\), \(n_g\), \(\alpha_{\rm eff}\), or \(\xi_c\).
- Frozen v14 paper and SPARC residual products were not modified.

---

## 6. Predictive paragraph and go/no-go

**Predictive insertion paragraph: withheld.**

The frozen intermediate-regime UMM form is classified as **ruled out** by this first confrontation within the quoted uncertainties and the stated estimator assumptions. A ready-to-use predictive note for a future paper revision is therefore **not** supplied. The scientific value of the exercise is the transparent tension: pure baryons match the DM-deficient UDG kinematics, while \(a_{\rm ind}=A(\Sigma/\Sigma_{\rm ref})^{0.4}\) at the frozen SPARC amplitude overpredicts \(\sigma\) for every system in the sample. Any future revision that retains this amplitude for disk galaxies must address why the same local prescription fails in these low-dispersion spheroids (geometry of the specialization, validity of mean-\(\Sigma\) proxy, or a genuine limitation of the intermediate-regime form at this scale).

### Is the result sound enough to merit a predictive note in a future revision?

**No** — not as a positive predictive claim. The confrontation is sound enough to *document tension* and to motivate either a careful statement of domain of applicability (thin-disk intermediate regime only) or a dedicated follow-up on the low-density / spheroidal sector. It is not sound enough to advertise UDG consistency with the frozen power-law form.

---

## 7. Residual-outcome summary (human go/no-go)

**Go/no-go for a predictive note: NO — withhold positive note.** On the residual sample (DF2, DF4, DF9, FCC224; N=4), median ε_σ(stars) = -0.11 while median ε_σ(UMM power-law) = +1.84. Classification: **ruled out**. For every system in this small sample, the frozen intermediate-regime UMM power-law predicts a line-of-sight dispersion more than 2σ above the observed value (model systematically high). The pure-baryon column is consistent with the data within uncertainties. This is an illustrative first estimate on a small sample; systematics remain significant, but the direction of the residual is uniform. Illustrative first estimate only; small sample; systematics remain significant. Do not re-fit A or α; do not claim a large statistical population.

---

*Report generated by `analysis/udg_confrontation/src/run_confrontation.py`. Numerical products: `results/`. Re-run: `python -m analysis.udg_confrontation.src.run_confrontation`.*
