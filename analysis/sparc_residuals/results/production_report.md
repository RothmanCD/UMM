# UMM SPARC Residual-Table Production Report

**Timestamp (UTC):** 2026-08-08T04:01:07.878276+00:00
**Selection ID:** `ea6c92761ffd3c5a`
**Design contract:** `notes/UMM_SPARC_Residual_Table_Design.md` §3–§9

## Status of numbers

These residuals are **measured catalog residuals** from a documented run of the
public SPARC mass models (Lelli, McGaugh & Schombert 2016) with pre-registered
selection cuts. They **replace** the illustrative schematic sample-summaries of
frozen `paper/UMM_Paper_Draft_v13.tex` for use in a future versioned successor
(v14+). The frozen v13 body was **not** edited in this run.

## Sample

- **N_gal:** 88
- **Coverage cut:** R_last ≥ 2.5 R_disc
- **Cut tallies:** `{"n_total": 175, "n_pass": 88, "fail_Q": 76, "fail_incl": 12, "fail_NR": 32, "fail_coverage": 14, "fail_meta": 0}`
- **Optional robustness cut (§3.2):** not applied as a rejection cut (no galaxies dropped by ΔV quality filter).

## Fitted global amplitude

- **A = 1.4125e-11 ± 7.0627e-13 m s^-2** (n_g = 0.4, Υ_* = 0.5)
- Outer-disk χ²_min = 31677.23
- Sensitivity (n_g = 0.45 packing index): A ≈ 1.3490e-11 m s^-2

## Median outer-disk residuals

| Model | med(ε) |
|-------|--------|
| UMM (single global A) | **0.0979** |
| Pure baryons (A=0) | **0.5068** |
| Simple NFW (1 mass/galaxy) | **0.0118** |

- Fraction with ε_umm < 0.25: **0.852**

## Optional robustness: Υ_* = 0.7

- A(Υ=0.7) = 1.2589e-11 ± 6.2946e-13 m s^-2
- med(ε_umm, Υ=0.7) = 0.0905
- **Δε = med(ε_umm,0.7) − med(ε_umm,0.5) = -0.0074**

## Deliverables

- `analysis/sparc_residuals/results/selection_manifest.json`
- `analysis/sparc_residuals/results/residuals_main.csv`
- `analysis/sparc_residuals/results/residuals_main.json`
- `analysis/sparc_residuals/results/residuals_summary.json`
- `analysis/sparc_residuals/run_pipeline.py` (reproducible entrypoint)

## Skeptic notes

- Selection frozen before residual ranking; no high-ε post-hoc drops.
- UMM column uses one global A for all galaxies.
- NFW baseline is one-parameter (M200) with Dutton–Macciò c(M), not full ΛCDM.
- Σ_gas uses Vgas thin-sheet proxy (SPARC_HI); not THINGS map re-reduction.

