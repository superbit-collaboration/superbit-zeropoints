"""
Extend Gaia XP continuous spectra shortward of their reliable ~336 nm floor
using a Teff/logg/[Fe/H]-matched Castelli-Kurucz (ck04models) synthetic
spectrum, rescaled to match the real XP flux in an overlap window.

Needed because SuperBIT's u-band throughput extends into the ~300-336 nm
range where Gaia XP spectra are not considered reliable.
"""
import os
import numpy as np
from astropy.table import Table

from config import load_config

os.environ.setdefault("PYSYN_CDBS", load_config()["cdbs_dir"])
import astropy.units as u  # noqa: E402
import stsynphot  # noqa: E402  (import after PYSYN_CDBS is set)

FLUX_UNIT = u.Unit("W m-2 nm-1")


def build_stellar_param_lookup(gaia_catalog=None, targets=None):
    """
    source_id -> (teff, feh, logg), from the merged Gaia catalog
    (config.yaml's gaia_catalog). Pass `targets` to restrict to a subset of
    the catalog's TARGET column; default is every row in the catalog.
    """
    if gaia_catalog is None:
        gaia_catalog = load_config()["gaia_catalog"]

    t = Table.read(gaia_catalog)
    if targets is not None:
        t = t[np.isin(t["TARGET"], list(targets))]

    source = np.asarray(t["Source"], dtype=np.int64)
    teff = np.ma.filled(np.ma.masked_invalid(t["Teff"]), np.nan)
    logg = np.ma.filled(np.ma.masked_invalid(t["logg"]), np.nan)
    feh = np.ma.filled(np.ma.masked_invalid(t["[Fe/H]"]), np.nan)

    valid = ~np.isnan(teff) & ~np.isnan(logg) & ~np.isnan(feh)

    lookup = {}
    for sid, tf, fh, lg in zip(source[valid], teff[valid], feh[valid], logg[valid]):
        lookup[int(sid)] = (float(tf), float(fh), float(lg))
    return lookup


def get_model_spectrum(teff, feh, logg):
    """ck04models SourceSpectrum for the given stellar parameters."""
    return stsynphot.grid_to_spec("ck04models", teff, feh, logg)


def get_model_flux(teff, feh, logg, wave_nm):
    """ck04models flux (W/m^2/nm) at the given wavelengths (nm)."""
    sp = get_model_spectrum(teff, feh, logg)
    return sp(np.asarray(wave_nm) * u.nm, flux_unit=FLUX_UNIT).value


def extend_spectrum_uv(flux, sampling, teff, feh, logg,
                        uv_min=300, overlap=(336, 360)):
    """
    Splice a ck04models spectrum (scaled to match real flux in `overlap`)
    onto the blue end of a real Gaia XP spectrum.

    Returns (full_wave, full_flux) covering [uv_min, sampling.max()].
    """
    sampling = np.asarray(sampling)
    flux = np.asarray(flux)

    uv_wave = np.arange(uv_min, sampling.min(), 1)
    if len(uv_wave) == 0:
        return sampling, flux

    lo, hi = overlap
    overlap_mask = (sampling >= lo) & (sampling <= hi)
    if not overlap_mask.any():
        raise ValueError(f"No real data in overlap window {overlap}")

    # one grid_to_spec call per star instead of two -- grid interpolation
    # (reading + interpolating 8 corner FITS spectra) is the expensive part
    sp = get_model_spectrum(teff, feh, logg)
    overlap_wave = sampling[overlap_mask]
    model_overlap = sp(overlap_wave * u.nm, flux_unit=FLUX_UNIT).value
    scale = np.median(flux[overlap_mask]) / np.median(model_overlap)

    model_uv = sp(uv_wave * u.nm, flux_unit=FLUX_UNIT).value * scale

    full_wave = np.concatenate([uv_wave, sampling])
    full_flux = np.concatenate([model_uv, flux])
    return full_wave, full_flux
