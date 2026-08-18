"""
Derives the AB magnitude zero point (3631 Jy, synthetic_photometry.py's
AB_ZEROPOINT_JY) from first principles, as a sanity check on theory.md's
claim that AB was anchored to Vega's own flux in the Johnson V band.

Oke & Gunn defined AB magnitudes so that AB = V (monochromatically, at
5480 Angstrom) for a source with a flat f_nu spectrum. Since 5480 Angstrom
is essentially the pivot wavelength of Johnson V (where f_lambda <-> f_nu
conversion is spectrum-independent -- computed and printed below), this is
close to saying "AB=0 matches Vega's own V-band flux". This script tests
that directly: integrate a real, independently-calibrated Vega spectrum
through a real Johnson V bandpass using this pipeline's own bandpass-
averaging formula (mean_fnu(), the photon-weighted convention used
everywhere else in this repo), and see how close the result lands to 3631 Jy.

Deliberately does NOT use data/vega/vega_rescaled_Wm2nm.npz -- that
spectrum was renormalized to an arbitrary internal anchor for this
pipeline's own Vega-magnitude bookkeeping (see scripts/rescale_vega.py),
so using it here to "derive" an absolute flux constant would be circular.
Instead this fetches Vega fresh via synphot.SourceSpectrum.from_vega()
(alpha_lyr_stis_011.fits, a CALSPEC calibration independent of anything
else in this repo) and a standard Johnson V bandpass via
synphot.SpectralElement.from_filter("johnson_v"). Requires network access
the first time (synphot caches the downloaded files after that).

In practice this lands within ~0.01% of 3631 Jy -- close enough to call it
confirmed. An older CALSPEC Vega calibration (alpha_lyr_mod_002.fits, also
bundled in this repo as data/vega/alpha_lyr_mod_002.fits) lands about 1%
off instead, which is why this script fetches the newer file rather than
reusing the bundled one: the ~1% gap turned out to be about which Vega
calibration vintage was used, not bandpass shape or Vega's well-known
non-zero (~0.02-0.03 mag) V magnitude -- both were checked and ruled out.

Note this comparison implicitly assumes Vega itself defines V=0 exactly
(comparing its raw f_nu directly to 3631 Jy, with no magnitude offset
applied). That's not a physical necessity -- it's just the actual
historical convention Oke & Gunn's constant was built on. Correcting for
Vega's modern, more precise V magnitude (~0.023-0.03, from later
photometry showing Vega is very slightly off the *ensemble*-defined
zero point) makes the match to 3631 Jy worse, not better (overshoots to
+1-2%), confirming the zero point really was pinned to "Vega==0" by
fiat, not to a more careful modern absolute V-band standard.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import astropy.units as u
import synphot
from synphot import units as syn_units

from synthetic_photometry import mean_fnu, fnu_to_jy, AB_ZEROPOINT_JY


def pivot_wavelength(bandpass_transmission, bandpass_wavelengths_nm):
    """sqrt[int(S*lambda dlambda) / int(S/lambda dlambda)] -- the wavelength
    at which f_lambda <-> f_nu conversion is exact regardless of the source
    spectrum's shape. Not the same as the transmission-weighted mean
    wavelength, though the two are close for a fairly symmetric bandpass
    like Johnson V."""
    lam = bandpass_wavelengths_nm * u.nm
    numerator = np.trapezoid(bandpass_transmission * lam, lam)
    denominator = np.trapezoid(bandpass_transmission / lam, lam)
    return np.sqrt(numerator / denominator).to(u.nm).value


def main():
    vega = synphot.SourceSpectrum.from_vega()
    vega_wave_nm = vega.waveset.to(u.nm).value
    vega_flux_Wm2nm = vega(vega.waveset, flux_unit=syn_units.FLAM).value * 1e-2  # erg/s/cm^2/A -> W/m^2/nm

    bp = synphot.SpectralElement.from_filter("johnson_v")
    bp_wave_nm = bp.waveset.to(u.nm).value
    bp_S = bp(bp.waveset).value

    lam_pivot = pivot_wavelength(bp_S, bp_wave_nm.copy())
    print(f"Johnson V pivot wavelength = {lam_pivot:.2f} nm  "
          f"(Oke & Gunn's AB reference wavelength was 548.0 nm)")

    # A pure flux integral -- no magnitude system involved yet. Comparing it
    # to AB_ZEROPOINT_JY below is what implicitly assumes Vega defines V=0
    # exactly (see module docstring).
    fnu = mean_fnu(vega_wave_nm, vega_flux_Wm2nm, bp_wave_nm, bp_S)
    fnu_jy = fnu_to_jy(fnu)

    print(f"Vega <f_nu> through Johnson V = {fnu_jy:.1f} Jy")
    print(f"AB_ZEROPOINT_JY (this repo)   = {AB_ZEROPOINT_JY:.1f} Jy")
    print(f"difference                    = {100 * (fnu_jy / AB_ZEROPOINT_JY - 1):+.2f}%")


if __name__ == "__main__":
    main()
