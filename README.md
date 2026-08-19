# Universal Manifold Model (UMM)

Manuscript draft v23 and the frozen SPARC residual products that support it.

**Paper:** Carl Rothman, *Residual Geometric Charges and a Continuous Density-Dependent Contribution to Gravity from a Single Compact Manifold* (draft v23, 19 August 2026).

**PDF:** [`paper/UMM_Paper_Draft_v23.pdf`](paper/UMM_Paper_Draft_v23.pdf)

This public repository is the supplementary data and code location named in the manuscript.

**Zenodo preprint:** [https://doi.org/10.5281/zenodo.22013917](https://doi.org/10.5281/zenodo.22013917)

## Abstract

We present a geometric framework—the Universal Manifold Model (UMM)—in which residual operators of the algebraic type associated with spin and a continuous, density-dependent contribution to gravitational acceleration arise as projections of a single object: independent relational flux of amplitude Φ living on a continuous compact manifold that threads every point of ordinary three-space.

On the toy manifolds S¹ and S³ one obtains residual operators of the correct algebraic type and overall scale Φ/R; a packing argument yields a logistic amplitude whose intermediate-regime index is compatible with the observed galactic index; joint normalization produces a window consistent with horizon and low-density constraints. Recovery of the full Standard-Model spectrum and a pure first-principles conversion for pressure-supported systems remain open.

## Frozen SPARC contact

| Quantity | Value |
|----------|-------|
| Selection | pre-registered cuts on Lelli, McGaugh & Schombert (2016) mass models |
| Sample | N = 88 galaxies |
| Selection id | `ea6c92761ffd3c5a` |
| Global amplitude | A = (1.41 ± 0.07) × 10⁻¹¹ m s⁻² (Υ₊ = 0.5, n_g = 0.4) |
| med(ε_UMM) | 0.098 |
| med(ε_baryon) | 0.507 |
| med(ε_NFW) | 0.012 |
| Fraction with ε_UMM < 0.25 | 0.85 |

These numbers are measured catalog products of `analysis/sparc_residuals/`. They are not retuned in later notes.

## Layout

| Path | Contents |
|------|----------|
| `paper/UMM_Paper_Draft_v23.tex` | Current manuscript source |
| `paper/UMM_Paper_Draft_v23.pdf` | Compiled PDF |
| `supplement/UMM_Residual_Geometry_Layer.md` | Companion Residual-Geometry Layer note (pure mathematical interface) |
| `supplement/UMM_Residual_Geometry_Layer_Additive_Note_2026-08-11.md` | Additive domain-extension note (2026-08-11) |
| `figures/` | Figure TeX, PDF, PNG, and `umm_style.tex` |
| `Data/SPARC/` | Public SPARC parent catalog (third-party; see that README) |
| `analysis/sparc_residuals/` | Residual pipeline and frozen N=88 products |
| `analysis/udg_confrontation/` | Domain-extension confrontation for pressure-supported UDGs (Fig. 7) |
| `notes/UMM_SPARC_Residual_Table_Design.md` | Pipeline design contract |
| `tests/test_sparc_pipeline.py` | Tests for the SPARC pipeline |

## Compile the paper

Figures are loaded with `\graphicspath{{./}{../figures/}}`.

```bash
cd paper
pdflatex -interaction=nonstopmode UMM_Paper_Draft_v23.tex
pdflatex -interaction=nonstopmode UMM_Paper_Draft_v23.tex
```

To rebuild figures from TikZ sources:

```bash
cd figures
./build_figures.sh
```

Fig. 7 (`fig7_residual_ratio_conversion`) is shipped as PDF/PNG only.

## Residual-Geometry Layer (supplement)

v23 isolates the pure mathematical construction in a companion note so it can be read independently of the galactic phenomenology. The frozen interface and the 2026-08-11 additive paragraph (residual-geometry-constrained conversion for pressure-supported systems) are:

- [`supplement/UMM_Residual_Geometry_Layer.md`](supplement/UMM_Residual_Geometry_Layer.md)
- [`supplement/UMM_Residual_Geometry_Layer_Additive_Note_2026-08-11.md`](supplement/UMM_Residual_Geometry_Layer_Additive_Note_2026-08-11.md)

## Reproduce the SPARC residual table

Python 3.10+; standard library only (see `requirements.txt`).

```bash
python analysis/sparc_residuals/run_pipeline.py
```

Expected frozen products (already in the tree; the pipeline regenerates them):

- `analysis/sparc_residuals/results/selection_manifest.json`
- `analysis/sparc_residuals/results/residuals_main.csv`
- `analysis/sparc_residuals/results/residuals_summary.json`
- `analysis/sparc_residuals/results/production_report.md`

Sanity tests:

```bash
python -m unittest tests.test_sparc_pipeline
```

The shipped pipeline is a Python package (`run_pipeline.py` and modules). It is not a Jupyter notebook.

## Data provenance

`Data/SPARC/` redistributes the public SPARC catalog of Lelli, McGaugh & Schombert (2016) so the residual pipeline can be run without a separate download. Cite:

> Lelli, F., McGaugh, S. S., & Schombert, J. M. 2016, AJ, 152, 157.

Source: http://astroweb.case.edu/SPARC/

This repository does not claim ownership of the SPARC catalog. Details: [`Data/SPARC/README.md`](Data/SPARC/README.md).

## License

- Analysis code and tests: MIT (see `LICENSE`).
- Manuscript, figures, and `supplement/` notes: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- SPARC catalog files: third-party; see `Data/SPARC/README.md`.

## Citation

Rothman, C. (2026). *Residual Geometric Charges and a Continuous Density-Dependent Contribution to Gravity from a Single Compact Manifold* (v23). Zenodo. https://doi.org/10.5281/zenodo.22013917

See also `CITATION.cff`.

## AI assistance

The manuscript includes a *Note on Collaboration and AI Assistance*. Intensive formalization used the AI system Grok (xAI). Grok is not a co-author; scientific responsibility rests with the human author.
