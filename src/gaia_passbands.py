"""
Gaia EDR3 official G/BP/RP passbands and zero points (Riello et al. 2021),
for validating synthetic photometry against real Gaia catalog magnitudes.

Data in <repo>/data/gaia_passbands/ (path from config.yaml's
gaia_passband_dir), originally from:
https://cdn.gea.esac.esa.int/Gaia/gedr3/Photometry/GaiaEDR3_passbands_zeropoints_version2
"""
import os
import numpy as np
import pandas as pd

from config import load_config


def load_gaia_passbands(data_dir=None):
    """Returns a dict: band -> (wave_nm, transmission), NaNs dropped, for G/BP/RP."""
    if data_dir is None:
        data_dir = load_config()["gaia_passband_dir"]

    pb = pd.read_csv(
        os.path.join(data_dir, "passband.dat"), sep=r"\s+", header=None,
        names=["wave_nm", "G", "e_G", "BP", "e_BP", "RP", "e_RP"],
    ).replace(99.99, np.nan)

    passbands = {}
    for band in ["G", "BP", "RP"]:
        sub = pb.dropna(subset=[band])
        passbands[band] = (sub["wave_nm"].values, sub[band].values)
    return passbands


def load_gaia_zeropoints(data_dir=None):
    """Returns a dict: system ('AB' or 'VEGAMAG') -> {band: zp}."""
    if data_dir is None:
        data_dir = load_config()["gaia_passband_dir"]

    zp = pd.read_csv(
        os.path.join(data_dir, "zeropt.dat"), sep=r"\s+", header=None,
        names=["G", "e_G", "BP", "e_BP", "RP", "e_RP", "system"],
    )
    return {
        row["system"]: {"G": row["G"], "BP": row["BP"], "RP": row["RP"]}
        for _, row in zp.iterrows()
    }


def vega_ab_offset(band, zeropoints=None):
    """
    mag_AB - mag_VEGAMAG for the given Gaia band. Unit-independent (the
    absolute flux-unit convention Riello et al. assumed for their zp cancels
    out), so this is a reliable way to convert an AB magnitude computed from
    first principles (e.g. via synthetic_photometry.mean_fnu/ab_mag using the
    real Gaia passband) into the Vega-like system that Gaia catalog
    magnitudes (Gmag, BPmag, RPmag) are actually reported in.
    """
    if zeropoints is None:
        zeropoints = load_gaia_zeropoints()
    return zeropoints["AB"][band] - zeropoints["VEGAMAG"][band]
