# SPARC catalog files (third-party)

This directory redistributes public SPARC (Spitzer Photometry and Accurate
Rotation Curves) catalog files so the N=88 residual pipeline in
`analysis/sparc_residuals/` can be run without a separate download.

These files are **not** authored or licensed by Carl Rothman. They remain
the SPARC team's data product.

## Citation (required)

If you use these files in a publication, cite the SPARC master paper:

> Lelli, F., McGaugh, S. S., & Schombert, J. M. 2016,
> *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry
> and Accurate Rotation Curves*, AJ, 152, 157.

ADS: https://ui.adsabs.harvard.edu/abs/2016AJ....152..157L

Original download site: http://astroweb.case.edu/SPARC/

Please consider citing the original HI / Hα rotation-curve sources listed
in `SPARC_Lelli2016c.mrt` when they are relevant.

## Contents

| Path | Description |
|------|-------------|
| `SPARC_Lelli2016c.mrt` | Table 1: galaxy sample metadata |
| `MassModels_Lelli2016c.mrt` | Combined mass-model rotation curves |
| `rotmod/` | Per-galaxy `*_rotmod.dat` mass models (175 galaxies) |
| `sfb/` | Per-galaxy `*.sfb` surface-brightness profiles |

The UMM residual analysis uses Table 1 plus the mass models. Selection
cuts are documented in `analysis/sparc_residuals/results/selection_manifest.json`
(selection_id `ea6c92761ffd3c5a`).
