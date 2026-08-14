"""
Documents how data/vega/vega_rescaled_Wm2nm.npz was produced from the raw
CALSPEC Vega model (data/vega/alpha_lyr_mod_002.fits, Bohlin 2014). Ships
already-generated; rerun this only if you want to change
vega_target_f550_Wm2nm in config.yaml.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from astropy.io import fits
from config import load_config


def main():
    cfg = load_config()

    with fits.open(cfg["vega_spectrum_file"]) as hdul:
        data = hdul["SCI"].data

    wave_A = data["WAVELENGTH"]
    flux_flam = data["FLUX"]  # erg/s/cm^2/A (CALSPEC FLAM convention)

    wave_nm = wave_A / 10.0
    flux_Wm2nm = flux_flam * 1e-2  # erg/s/cm^2/A -> W/m^2/nm

    current_550 = np.interp(550.0, wave_nm, flux_Wm2nm)
    target_550 = cfg["vega_target_f550_Wm2nm"]
    scale = target_550 / current_550
    flux_rescaled = flux_Wm2nm * scale

    print(f"current f(550nm)  = {current_550:.6e} W/m^2/nm")
    print(f"target f(550nm)   = {target_550:.6e} W/m^2/nm")
    print(f"scale factor       = {scale:.6f}")

    np.savez(cfg["vega_rescaled_file"], wave_nm=wave_nm, flux=flux_rescaled)
    print(f"Wrote {cfg['vega_rescaled_file']}")


if __name__ == "__main__":
    main()
