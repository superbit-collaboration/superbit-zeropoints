# SuperBIT zeropoints

Computes SuperBIT's u, b, g magnitude zeropoints from Gaia DR3 XP spectra.
**Read `theory.md` first** for how and why this works.

## Contributors

Most of the utilities for estimating synthetic Vega magnitudes from Gaia
spectra were developed by **Thuy Vy Luu (Princeton)**, including the
self-calibration against the Gaia photometric catalog.

The AB magnitude zero points were estimated independently by
**Emaad Paracha (UToronto)** and **Ajay Gill (UToronto)**, the latter as
part of his thesis (see [Section 3.8](https://utoronto.scholaris.ca/items/9b6ae5fe-04cd-44bf-be63-67e6ea6dd68f)). The AB synthetic
magnitude utilities in this repo are imported from his photometry code
([link](https://github.com/superbit-collaboration/superbit_photometry/blob/master/photometry.py)). Both Emaad's and Ajay's zero points are
consistent with the final version implemented here.

This repository was created to assemble everyone's work into a single,
refactored pipeline for future reproducibility.

## Setup

```bash
git clone git@github.com:superbit-collaboration/superbit-zeropoints.git && cd superbit-zeropoints
pip install -r requirements.txt
bash scripts/download_ck04models.sh   # ~41MB stellar atmosphere grid
```

Then get the large data this repo doesn't ship in git (SuperBIT star
catalogs + Gaia XP spectra, ~200MB) into `external_data/`. This isn't
reachable from every filesystem, so grab it from either:

- hen: `/data/analysis/superbit_2023/sayan/zp_data/external_data`
- [Google Drive](https://drive.google.com/drive/folders/1t-GOyoJ9phlzvC2cKMEfg3B_TpKp-xte)

## Running it

```bash
python scripts/validate_against_gaia.py   # sanity check (optional but recommended)
python src/compute_sb_zeropoints.py       # main pipeline, ~25 min first run
```

Writes per-star tables and diagnostic plots to `outputs/`. Edit
`config.yaml` to change the flux column, aggregation method, or paths.

## Result

| Band | ZP (AB) | ZP (Vega) |
|------|---------|-----------|
| u | 28.687 ± 0.005 | 28.600 ± 0.005 |
| b | 30.228 ± 0.001 | 30.294 ± 0.001 |
| g | 29.614 ± 0.002 | 29.537 ± 0.002 |

(jackknifed across 30 target fields; see `theory.md` for why that error
bar, not a flat per-star one, is the one to trust)

![Final zeropoint diagnostic](docs/figures/sb_zeropoints_diagnostic.png)

