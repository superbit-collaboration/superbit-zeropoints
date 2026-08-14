"""
Build a single merged Gaia DR3 catalog from the per-target cached catalogs
used during development, keeping only the columns this pipeline actually
needs (the source per-target files carry ~200 Gaia columns; this keeps 9).

This script documents how data/gaia_dr3_merged.fits (referenced by
config.yaml) was produced. You only need to rerun it if you're adding a new
target -- point SOURCE_DIR at a directory of "<TARGET>_gaia_dr3.fits" files
(as returned by a Gaia DR3 cone-search query, e.g. via astroquery/Vizier).
"""
import sys
import numpy as np
from astropy.table import Table, vstack

SOURCE_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
OUT_FILE = sys.argv[2] if len(sys.argv) > 2 else "data/gaia_dr3_merged.fits"

TARGETS = [
    "1E0657_Bullet", "Abell1689", "Abell2163", "Abell2345", "Abell2384a",
    "Abell2384b", "Abell3365", "Abell3411", "Abell3526", "Abell3571",
    "Abell3667", "Abell3716S", "Abell3827", "AbellS0592", "AbellS780",
    "COSMOS113", "COSMOSk", "MACSJ0416d1m2403", "MACSJ0723d3_7327_JWST",
    "MACSJ1105d7m1014", "MACSJ1931d8m2635", "MS1008d1m1224", "MS2137m2353",
    "PLCKG287d0p32d9", "RXCJ1314d4m2515", "RXCJ1514d9m1523",
    "RXCJ2003d5m2323", "SMACSJ2031d8m4036", "SPTCLJ0411",
    "Z20_SPT_CLJ0135m5904",
]

KEEP_COLS = [
    "Source", "ALPHAWIN_J2000", "DELTAWIN_J2000", "XPcont",
    "Teff", "logg", "[Fe/H]", "Gmag", "BPmag", "RPmag",
]


def main():
    tables = []
    for target in TARGETS:
        path = f"{SOURCE_DIR}/{target}_gaia_dr3.fits"
        t = Table.read(path)[KEEP_COLS]
        t["TARGET"] = target
        tables.append(t)
        print(f"[{target}] {len(t)} rows")

    merged = vstack(tables, join_type="exact")
    merged.write(OUT_FILE, overwrite=True)
    print(f"\nWrote {OUT_FILE}: {len(merged)} rows, {len(merged.colnames)} columns")


if __name__ == "__main__":
    main()
