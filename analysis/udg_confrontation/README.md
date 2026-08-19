# UMM confrontation with DM-deficient UDG kinematics

Illustrative first estimate of the frozen v14 intermediate-regime residual
acceleration against published kinematics of dark-matter-deficient ultra-diffuse
galaxies (NGC 1052-DF2, DF4, DF9, FCC 224).

**Honesty:** small sample; systematics remain significant; no re-fit of \(A\),
\(n_g\), \(\alpha_{\rm eff}\), or \(\xi_c\); no large-population claim.

## Layout

```
analysis/udg_confrontation/
  README.md                 # this file
  src/
    physics.py              # Wolf estimator, a_ind power-law & logistic
    manifest_data.py        # every adopted literature number + citation
    run_confrontation.py    # entry point: predictions + report
  results/
    confrontation_report.md # referee-readable main report
    literature_manifest.json
    selection_manifest.json
    udg_data_table.csv
    udg_predictions_three_models.csv
    residual_summary.json
    sigma_comparison.txt
  tests/
    test_udg_physics.py
```

## Frozen inputs (do not change)

| Quantity | Value | Source |
|----------|-------|--------|
| \(A\) | \((1.41\pm 0.07)\times 10^{-11}\,\mathrm{m\,s^{-2}}\) | SPARC residual fit, selection_id `ea6c92761ffd3c5a` |
| \(n_g\) | 0.4 | v14 / RAR comparison |
| \(\Sigma_{\rm ref}\) | \(1\,M_\odot\,\mathrm{pc^{-2}}\) | v14 galactic specialization |
| \(\Upsilon_*\) | 0.5 (SPARC context) | residual table design |
| \(\alpha_{\rm eff}\) | ≈1.8 | packing note / v14 |
| \(\Sigma_{\rm mid}\) window | {0.3, 1.0, 3.0} | joint-norm midpoint about \(\Sigma_{\rm ref}\) |

## Run

From the repository root:

```bash
python -m analysis.udg_confrontation.src.run_confrontation
python -m unittest analysis.udg_confrontation.tests.test_udg_physics -v
```

## What is not touched

- `paper/UMM_Paper_Draft_v14.tex`
- `analysis/sparc_residuals/results/*`

## Estimator (one-line)

Wolf (2010): \(M_{1/2}=4\sigma^2 R_e/G\), with \(a_{\rm tot}=a_{\rm bar}+a_{\rm ind}\) and
\(\sigma=\sqrt{a_{\rm tot} R_e/4}\). Mean \(\Sigma\) within \(R_e\) is used when no
published mass profile exists (limitation stated in the report).
