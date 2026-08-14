"""
Sanity check (run this before trusting the SuperBIT zeropoints): validate
synthetic Vega-system photometry -- vega_mag(), through Gaia's own G/BP/RP
passbands plus our rescaled reference Vega spectrum -- against real Gaia
catalog magnitudes (Gmag/BPmag/RPmag), for every star with a downloaded XP
continuous spectrum.

Deliberately does NOT use data/gaia_passbands/zeropt.dat -- this is the
"from first principles" check (spectrum + Gaia's own passband + our own
Vega spectrum only, all going through synthetic_photometry.py's vega_mag()),
not a shortcut through Gaia's own calibration constants. If this doesn't
come out close to a 1:1 line, don't trust the SuperBIT zeropoints either,
since they're computed with the exact same machinery just swapped to
SuperBIT's bandpasses.

Writes outputs/validate_against_gaia.fits and outputs/validate_against_gaia.png.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import time
import numpy as np
from astropy.table import Table
from gaiaxpy import calibrate

from config import load_config, REPO_ROOT
from xp_uv_extrapolation import build_stellar_param_lookup, extend_spectrum_uv
from synthetic_photometry import vega_mag
from gaia_passbands import load_gaia_passbands

BANDS = ["G", "BP", "RP"]
CAT_MAG_COL = {"G": "Gmag", "BP": "BPmag", "RP": "RPmag"}
OUTPUT_DIR = os.path.join(REPO_ROOT, "outputs")


def build_photometry_lookup(gaia_catalog):
    """source_id -> dict(Gmag, BPmag, RPmag) from the merged Gaia catalog."""
    t = Table.read(gaia_catalog)
    source = np.asarray(t["Source"], dtype=np.int64)
    mags = {band: np.ma.filled(np.ma.masked_invalid(t[col]), np.nan)
            for band, col in CAT_MAG_COL.items()}
    return {
        int(sid): {band: float(mags[band][i]) for band in BANDS}
        for i, sid in enumerate(source)
    }


def main():
    cfg = load_config()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Building stellar parameter lookup...")
    stellar_params = build_stellar_param_lookup(cfg["gaia_catalog"])
    print(f"  {len(stellar_params)} stars with valid Teff/logg/[Fe/H]")

    print("Building photometry lookup...")
    photometry = build_photometry_lookup(cfg["gaia_catalog"])

    passbands = load_gaia_passbands()

    vega_data = np.load(cfg["vega_rescaled_file"])
    vega_wave, vega_flux = vega_data["wave_nm"], vega_data["flux"]

    print("Calibrating XP spectra...")
    sampling = np.arange(336, 1050, 1)
    spectra, sampling = calibrate(cfg["xp_spectra_ecsv"], sampling=sampling, save_file=False)

    n = len(spectra)
    results = {"source_id": [], "extended": []}
    for band in BANDS:
        results[f"synth_{band}"] = []
        results[f"cat_{band}"] = []

    t0 = time.time()
    for i, row in spectra.iterrows():
        sid = int(row["source_id"])
        wave, flux = sampling, row["flux"]
        extended = False

        if sid in stellar_params:
            teff, feh, logg = stellar_params[sid]
            try:
                wave, flux = extend_spectrum_uv(
                    flux, sampling, teff, feh, logg,
                    uv_min=cfg["uv_extension_min_nm"],
                    overlap=tuple(cfg["uv_extension_overlap_nm"]),
                )
                extended = True
            except Exception as e:
                print(f"[WARNING] extension failed for {sid}: {e}")

        cat_mags = photometry.get(sid, {b: np.nan for b in BANDS})

        results["source_id"].append(sid)
        results["extended"].append(extended)
        for band in BANDS:
            bp_wave, bp_S = passbands[band]
            try:
                synth = vega_mag(wave, flux, bp_wave, bp_S, vega_wave, vega_flux)
            except Exception:
                synth = np.nan
            results[f"synth_{band}"].append(synth)
            results[f"cat_{band}"].append(cat_mags[band])

        if (i + 1) % 200 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (n - i - 1) / rate
            print(f"  {i+1}/{n}  ({rate:.1f} stars/s, ETA {eta/60:.1f} min)")

    out = Table(results)
    outfile = os.path.join(OUTPUT_DIR, "validate_against_gaia.fits")
    out.write(outfile, overwrite=True)
    print(f"Wrote {outfile} ({len(out)} rows)")

    make_plot(out)


def make_plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["text.usetex"] = False  # gaiaxpy forces this True on import

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, band in zip(axes, BANDS):
        synth = np.asarray(out[f"synth_{band}"])
        cat = np.asarray(out[f"cat_{band}"])
        mask = np.isfinite(synth) & np.isfinite(cat)
        synth, cat = synth[mask], cat[mask]
        resid = synth - cat
        med, std = np.median(resid), np.std(resid)

        ax.scatter(cat, resid, s=8, alpha=0.4, color="tab:green")
        ax.axhline(0, color="k", lw=1, ls="--")
        pad = max(4 * std, 0.01)
        ax.set_ylim(med - pad, med + pad)
        ax.set_xlabel(f"Catalog {band} mag")
        ax.set_ylabel(f"Synthetic Vega mag - Catalog {band} mag")
        ax.set_title(
            f"{band}: N={mask.sum()}\n"
            f"median resid={med:+.3f}, std={std:.3f}"
        )

    plt.tight_layout()
    outname = os.path.join(OUTPUT_DIR, "validate_against_gaia.png")
    plt.savefig(outname, dpi=150)
    print(f"Wrote {outname}")


if __name__ == "__main__":
    main()
