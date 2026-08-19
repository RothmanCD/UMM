# UMM SPARC Residual-Table Design  
## Documented public-catalog analysis specification (methods supplement)

**Status:** Design / ready-to-implement specification  
**Parent baseline:** frozen `paper/UMM_Paper_Draft_v13.tex` (and byte-identical v13.1)  
**Date:** 2026-08-07  
**Language:** American English  

### Honesty (mandatory)

This note is **preparation only**. It does **not** claim that a public SPARC re-reduction has been executed. Residual numbers in v13 remain *illustrative schematic sample-summaries*, not measured catalog products. Until the pipeline below is run, versioned, and released with open code and data products, **no table entry may be presented as an empirical residual**.

Terminology matches v13: independent relational flux amplitude \(\Phi(\xi)\), intermediate-regime acceleration \(a_{\rm ind}\), single global amplitude \(A\), baryonic surface density \(\Sigma_{\rm bar}=\Sigma_*+\Sigma_{\rm gas}\).

---

## 1. Purpose

Convert the v13 illustrative SPARC contact into a **fully documented, reproducible** residual table for a high-quality public-catalog subsample, with:

- explicit selection cuts;
- fixed mass-to-light (M/L) and gas protocols;
- single global \(A\) fit (no per-galaxy dark halo in the UMM column);
- exact residual definition and baselines (pure baryons; simple NFW);
- release of code and intermediate data products.

**Forced (by this design):** the schema, cuts, and formulas below are the analysis contract.  
**Requires execution (extra data / work):** actual SPARC/THINGS files, compute run, and public release.

---

## 2. Public data products (inputs)

| Product | Source | Use |
|---------|--------|-----|
| SPARC mass models & rotation curves | Lelli, McGaugh & Schombert (2016), public SPARC database | \(V_{\rm obs}(R)\), stellar/gas disk models, quality flags, inclination, distance |
| SPARC photometry / \(3.6\,\mu\mathrm{m}\) disks | Same | Stellar surface density \(\Sigma_*\) |
| HI maps (THINGS-quality where available) | Walter et al. (2008) THINGS; SPARC-native HI otherwise | \(\Sigma_{\rm gas}\) for outer disk |
| NFW concentration prior | Dutton & Macciò (2014) mass–concentration (or equivalent published \(c(M)\) at fixed cosmology) | Simple NFW baseline only |

Primary archive path (to be pinned at run time): the official SPARC data release URL/DOI cited in Lelli et al. (2016), with a **frozen download SHA-256 manifest** stored in the analysis repository.

---

## 3. Selection cuts (high-quality subsample)

Target size: **\(N \approx 28\)** galaxies (v13 notional size), or a **clearly justified larger set** if cuts below admit more objects while remaining non-cherry-picked.

### 3.1 Mandatory cuts (apply in order)

1. **Quality flag:** SPARC \(Q = 1\) only.  
2. **Inclination:** \(i \geq 30^\circ\) (reject near face-on systems with large \(V\) deprojection uncertainty).  
3. **Outer-disk coverage:** last measured point \(R_{\rm last} \geq 2.5\,R_{\rm disc}\), where \(R_{\rm disc}\) is the SPARC exponential disk scale length (or equivalent published disk scale).  
4. **Rotation-curve sampling:** at least \(N_R \geq 8\) independent radial points with published \(V_{\rm obs}\) and uncertainties.  
5. **Distance / inclination metadata:** distance and inclination present and finite in the SPARC table.  
6. **Non-interaction / morphology note:** record morphological type (Hubble or SPARC type string); do **not** exclude by type except if the SPARC notes flag the curve as unusable for circular-speed modeling (document any such exclusion with SPARC note ID).

### 3.2 Optional robustness cut (report both with and without)

7. **\(\Delta V\) quality:** mean fractional uncertainty \(\langle \sigma_V / V_{\rm obs}\rangle_{R\geq R_{\rm out}} \leq 0.15\).

### 3.3 Sample construction rules

- Cuts are **pre-registered** by this note; no post-hoc removal of high-residual galaxies.  
- If \(N < 20\) after mandatory cuts, **widen only cut 3** to \(R_{\rm last} \geq 2.0\,R_{\rm disc}\) and report both samples.  
- If \(N > 40\), keep all objects that pass cuts (larger set is preferred to cherry-picking down to 28); still report a nested “strict-28” subsample as the first 28 by SPARC name sort **only if** a fixed size is required for table layout—prefer the full clean set.

**Forced:** cut list and non-cherry-pick rule.  
**Requires listed extra data:** the live SPARC table version/DOI.

---

## 4. Mass models

### 4.1 Stellar mass-to-light protocol

- Use SPARC \(3.6\,\mu\mathrm{m}\) stellar disk (and bulge if present) mass models.  
- **Primary M/L:** single global stellar mass-to-light ratio \(\Upsilon_*^{[3.6]}\) in the SPARC-preferred range  
  \[
  \Upsilon_*^{[3.6]} \in [0.5,\,0.7]\,M_\odot/L_\odot
  \]
  with **fiducial** \(\Upsilon_*^{[3.6]} = 0.5\) (disk) unless a published SPARC “preferred” column is adopted uniformly for the whole sample.  
- **Bulge:** if SPARC supplies a separate bulge \(\Upsilon_b\), use the SPARC default bulge value consistently for all galaxies that have bulges; do not fit per-galaxy \(\Upsilon\).  
- **Robustness:** re-run with \(\Upsilon_*^{[3.6]} = 0.7\) and quote \(\Delta\epsilon\) (shift in residuals) in the notes column / appendix.

**Forced (design):** one global \(\Upsilon_*\) policy for the sample, not per-galaxy free \(\Upsilon\).  
**Requires listed extra data:** choice of fiducial \(\Upsilon_*\) inside the published range.

### 4.2 Gas surface density

- \(\Sigma_{\rm gas} = 1.33\,\Sigma_{\rm HI}\) (helium correction factor **1.33**, standard SPARC convention) unless SPARC already provides He-corrected gas.  
- Prefer **THINGS-quality** HI surface-density maps for the outer disk when the galaxy is in THINGS; otherwise use SPARC-native gas disk models.  
- Always form  
  \[
  \Sigma_{\rm bar}(R) = \Sigma_*(R) + \Sigma_{\rm gas}(R).
  \]
- Report a **stellar-only** control residual \(\epsilon_i^{(*)}\) with \(\Sigma_{\rm bar}=\Sigma_*\) in an appendix table (not the main UMM column).

---

## 5. Model accelerations and single global \(A\)

### 5.1 Baryonic acceleration

From SPARC mass models, compute the Newtonian baryonic circular contribution \(V_{\rm bar}(R)\) (stars + gas ± bulge as published), and

\[
a_{\rm bar}(R) = \frac{V_{\rm bar}(R)^2}{R}.
\]

### 5.2 UMM independent term (intermediate regime)

Following v13,

\[
a_{\rm ind}(R) = A \left(\frac{\Sigma_{\rm bar}(R)}{\Sigma_{\rm ref}}\right)^{n_g},
\qquad
\Sigma_{\rm ref} = 1\,M_\odot\,\mathrm{pc^{-2}},
\qquad
n_g = 0.4
\]

(use \(n_g=0.4\) for direct comparison with the published radial acceleration relation; optionally report \(n_g=0.45\) from packing as a one-line sensitivity).

Total model:

\[
a_{\rm tot}(R) = a_{\rm bar}(R) + a_{\rm ind}(R),
\qquad
V_{\rm model}(R) = \sqrt{R\,a_{\rm tot}(R)}.
\]

### 5.3 Single global \(A\) fit

- Fit **one** scalar \(A\) to the full selected sample by minimizing

\[
\chi^2(A)
=
\sum_{i\in\mathrm{sample}}
\sum_{R_j \geq R_{{\rm out},i}}
\left(
\frac{V_{{\rm obs},i}(R_j) - V_{{\rm model},i}(R_j;A)}{\sigma_{V,i}(R_j)}
\right)^2
\]

(outer-disk points only for the fit that defines the published \(A\); report an all-radii fit as robustness).  

- No per-galaxy dark-matter halo parameters in the UMM model.  
- Quote \(A \pm \sigma_A\) from \(\Delta\chi^2 = 1\) (or bootstrap over galaxies).  
- **Target scale (not a measured result of this design note):** v13 joint-normalization window  
  \(A = (1.2\pm 0.2)\times 10^{-10}\,\mathrm{m\,s^{-2}}\). The pipeline will **replace** this with the fitted value and uncertainty.

**Forced:** single global \(A\); formula for \(V_{\rm model}\).  
**Requires execution:** the numerical \(A\) from data.

---

## 6. Outer-disk residual definition (exact)

Match v13:

\[
\epsilon_i
=
\left\langle
\left|
\frac{V_{\rm obs}(R) - V_{\rm model}(R)}{V_{\rm obs}(R)}
\right|
\right\rangle_{R \geq R_{\rm out}},
\]

with

\[
R_{\rm out}
=
\text{radius enclosing 70\% of the radial extent of measured points}
=
R_{\min} + 0.70\,(R_{\rm last}-R_{\min}),
\]

where \(R_{\min}\) is the innermost measured point used in the SPARC curve for that galaxy.  

- Average is the **unweighted arithmetic mean** over discrete SPARC radial points with \(R\geq R_{\rm out}\) (primary); report inverse-variance weighted mean as robustness.  
- Points with \(V_{\rm obs}\leq 0\) or missing \(\sigma_V\) are dropped and counted in `notes`.

### 6.1 Summary statistics (to be computed after the run)

- \(\mathrm{med}(\epsilon_i)\) over the sample  
- fraction with \(\epsilon_i < 0.25\)  
- stratified medians by mean outer-disk \(\Sigma_{\rm bar}\) (high / intermediate / low tertiles)

These statistics **must not** be filled with v13 illustrative numbers in any “results” table until the pipeline runs.

---

## 7. Comparison baselines

### 7.1 Pure baryons

\[
V_{\rm bar\mbox{-}only}(R) = V_{\rm bar}(R),
\qquad
\epsilon_i^{\rm bar}
=
\left\langle
\left|
\frac{V_{\rm obs}-V_{\rm bar}}{V_{\rm obs}}
\right|
\right\rangle_{R\geq R_{\rm out}}.
\]

Equivalent to \(A=0\).

### 7.2 Simple NFW halo

For each galaxy, add a spherical NFW halo with:

- **one free parameter:** virial mass \(M_{200}\) (or \(V_{200}\));  
- **concentration** \(c_{200} = c(M_{200})\) from the fixed Dutton–Macciò (or equivalent) mass–concentration relation at a stated cosmology (\(H_0\), \(\Omega_m\) pinned in the release);  
- fit \(M_{200}\) by minimizing outer-disk \(\chi^2\) for that galaxy alone (or full curve; state choice and use consistently);  
- residual \(\epsilon_i^{\rm NFW}\) with the same \(\epsilon\) definition as UMM.

**Forced (design):** NFW has one mass parameter per galaxy; UMM has one global \(A\).  
**Requires listed extra data:** choice of \(c(M)\) relation and cosmology.

---

## 8. Table schema (main residual table)

One row per galaxy. Column names are fixed for machine readability.

| Column | Type | Description |
|--------|------|-------------|
| `galaxy_id` | string | SPARC galaxy name (primary key) |
| `morph_type` | string | Morphological type (SPARC / RC3-style) |
| `Q` | int | SPARC quality flag (must be 1) |
| `incl_deg` | float | Inclination in degrees |
| `dist_Mpc` | float | Distance adopted |
| `R_disc_kpc` | float | Disk scale length |
| `R_last_kpc` | float | Last measured radius |
| `R_out_kpc` | float | Outer-disk threshold used for \(\epsilon\) |
| `N_outer` | int | Number of points with \(R\geq R_{\rm out}\) |
| `mean_Sigma_bar_outer` | float | Mean \(\Sigma_{\rm bar}\) on outer points \([M_\odot\,\mathrm{pc^{-2}}]\) |
| `stratum` | string | `high` / `intermediate` / `low` (tertile of `mean_Sigma_bar_outer`) |
| `epsilon_umm` | float | Outer-disk \(\epsilon_i\) for UMM (global \(A\)) |
| `epsilon_baryon` | float | Outer-disk \(\epsilon_i^{\rm bar}\) |
| `epsilon_nfw` | float | Outer-disk \(\epsilon_i^{\rm NFW}\) |
| `A_global` | float | Fitted global \(A\) (same for all rows) \([\mathrm{m\,s^{-2}}]\) |
| `Upsilon_star` | float | Stellar M/L used |
| `gas_source` | string | `THINGS` / `SPARC_HI` / `none` |
| `notes` | string | Free text: exclusions, missing gas, bulge flags, etc. |

### 8.1 Companion summary table (sample-level)

| Column | Description |
|--------|-------------|
| `N_gal` | Sample size |
| `med_epsilon_umm` | Median \(\epsilon_{\rm umm}\) |
| `med_epsilon_baryon` | Median \(\epsilon_{\rm baryon}\) |
| `med_epsilon_nfw` | Median \(\epsilon_{\rm nfw}\) |
| `frac_eps_lt_0.25` | Fraction with \(\epsilon_{\rm umm}<0.25\) |
| `A_global` / `A_err` | Fit and uncertainty |
| `selection_id` | Hash of cut list + data manifest |

---

## 9. Code and data release plan

When the pipeline is executed, release under an open license (e.g. MIT / BSD-3 for code; CC-BY for derived tables):

1. **Code repository** (public):  
   - environment lockfile (`environment.yml` or `requirements.txt` + Python version);  
   - script to download SPARC (and THINGS cross-match list) with checksum verification;  
   - selection module implementing §3;  
   - mass-model / \(a_{\rm bar}\) builder;  
   - UMM \(a_{\rm ind}\) + global-\(A\) fitter;  
   - NFW baseline fitter;  
   - table writer (CSV + machine-readable JSON);  
   - one-command `make residuals` target.
2. **Data products:**  
   - `selection_manifest.json` (galaxy list + cut booleans + SPARC file hashes);  
   - `residuals_main.csv` (schema §8);  
   - `residuals_summary.json` (§8.1);  
   - optional per-galaxy model curves (`V_obs`, `V_bar`, `V_umm`, `V_nfw`).
3. **Paper interface:**  
   - Replace v13 illustrative Table/Fig residual panels only in a **versioned successor manuscript** (v14+), with Methods pointing to this note and the release DOI.  
   - Retain a short sentence that pre-pipeline drafts used illustrative sample-summaries.

**Until release:** all prose must keep the v13 honesty language (illustrative schematic sample-summaries; not a catalog re-pipeline).

---

## 10. Implementation checklist (for the first production run)

- [ ] Pin SPARC + THINGS download URLs and SHA-256 manifests  
- [ ] Implement cuts §3; freeze `galaxy_id` list  
- [ ] Implement M/L + gas §4  
- [ ] Fit global \(A\) §5; freeze \(A\pm\sigma_A\)  
- [ ] Compute \(\epsilon_{\rm umm}\), \(\epsilon_{\rm baryon}\), \(\epsilon_{\rm nfw}\) §6–7  
- [ ] Write tables §8  
- [ ] Public tag + DOI §9  
- [ ] Only then: promote numbers from “illustrative” to “measured” in a successor paper  

---

## 11. Skeptic notes (scope control)

| Claim | Status |
|-------|--------|
| This note produces measured \(\epsilon_i\) | **False** — design only |
| Selection is non-cherry-picked | **True by contract** once cuts are frozen before looking at residuals |
| UMM uses one global \(A\) | **True by design** |
| NFW is a full ΛCDM halo analysis | **False** — simple one-parameter NFW baseline only |
| \(n_g=0.4\) is derived in this note | **False** — adopted for RAR comparison; packing index \(0.45\) is sensitivity |

---

## 12. References (inputs)

- Lelli, McGaugh & Schombert, AJ 152, 157 (2016) — SPARC  
- McGaugh, Lelli & Schombert, PRL 117, 201101 (2016) — radial acceleration relation  
- Walter et al., AJ 136, 2563 (2008) — THINGS  
- Dutton & Macciò, MNRAS 441, 3359 (2014) — mass–concentration (or successor relation pinned at run time)  
- UMM Paper Draft v13 §Galactic specialization / Limitations (honesty flags)

---

*End of UMM SPARC Residual-Table Design. Additive methods note; does not modify frozen v13.*
