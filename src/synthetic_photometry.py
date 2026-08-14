"""
Bandpass-averaged synthetic photometry from an f_lambda spectrum.
"""
import numpy as np

C_NM_PER_S = 2.99792458e17  # speed of light, nm/s
AB_ZEROPOINT_JY = 3631.0    # AB mag=0 flux density, Jy


def mean_fnu(wave_nm, flux_lambda, bp_wave_nm, bp_transmission):
    """
    Bandpass-averaged mean flux density in frequency units:

        <f_nu> = int f_lambda(lambda) S(lambda) lambda dlambda
                 / int S(lambda) (c/lambda) dlambda

    Parameters
    ----------
    wave_nm : array
        Wavelength grid for `flux_lambda`, in nm.
    flux_lambda : array
        f_lambda(lambda) on `wave_nm`, in W/m^2/nm.
    bp_wave_nm : array
        Wavelength grid of the bandpass throughput, in nm.
    bp_transmission : array
        Dimensionless throughput S(lambda) on `bp_wave_nm`.

    Returns
    -------
    float
        <f_nu>, in W/m^2/Hz.
    """
    wave_nm = np.asarray(wave_nm)
    flux_lambda = np.asarray(flux_lambda)
    bp_wave_nm = np.asarray(bp_wave_nm)
    bp_transmission = np.asarray(bp_transmission)

    lo = max(wave_nm.min(), bp_wave_nm.min())
    hi = min(wave_nm.max(), bp_wave_nm.max())
    mask = (bp_wave_nm >= lo) & (bp_wave_nm <= hi)
    if not mask.any():
        raise ValueError("Spectrum and bandpass wavelength ranges do not overlap")

    lam = bp_wave_nm[mask]
    S = bp_transmission[mask]
    f_lam = np.interp(lam, wave_nm, flux_lambda)

    numerator = np.trapezoid(f_lam * S * lam, lam)
    denominator = np.trapezoid(S * (C_NM_PER_S / lam), lam)
    return numerator / denominator


def mean_flambda(wave_nm, flux_lambda, bp_wave_nm, bp_transmission):
    """
    Bandpass-averaged mean flux density in wavelength units:

        <f_lambda> = int f_lambda(lambda) S(lambda) lambda dlambda
                     / int S(lambda) lambda dlambda

    Note the denominator here is int S*lambda dlambda (no speed-of-light
    factor) -- this is NOT the same normalization as mean_fnu(), since the
    integrated energy is being kept in f_lambda units rather than converted
    to f_nu.

    Parameters
    ----------
    wave_nm : array
        Wavelength grid for `flux_lambda`, in nm.
    flux_lambda : array
        f_lambda(lambda) on `wave_nm`, in W/m^2/nm.
    bp_wave_nm : array
        Wavelength grid of the bandpass throughput, in nm.
    bp_transmission : array
        Dimensionless throughput S(lambda) on `bp_wave_nm`.

    Returns
    -------
    float
        <f_lambda>, in the same units as `flux_lambda` (e.g. W/m^2/nm).
    """
    wave_nm = np.asarray(wave_nm)
    flux_lambda = np.asarray(flux_lambda)
    bp_wave_nm = np.asarray(bp_wave_nm)
    bp_transmission = np.asarray(bp_transmission)

    lo = max(wave_nm.min(), bp_wave_nm.min())
    hi = min(wave_nm.max(), bp_wave_nm.max())
    mask = (bp_wave_nm >= lo) & (bp_wave_nm <= hi)
    if not mask.any():
        raise ValueError("Spectrum and bandpass wavelength ranges do not overlap")

    lam = bp_wave_nm[mask]
    S = bp_transmission[mask]
    f_lam = np.interp(lam, wave_nm, flux_lambda)

    numerator = np.trapezoid(f_lam * S * lam, lam)
    denominator = np.trapezoid(S * lam, lam)
    return numerator / denominator


def vega_mag(wave_nm, flux_lambda, bp_wave_nm, bp_transmission,
             vega_wave_nm, vega_flux_lambda):
    """
    Synthetic Vega-system magnitude:

        m_VEGA = -2.5 log10(<f_lambda>)
                 + 2.5 log10( int f_lambda^Vega(lambda) S(lambda) lambda dlambda
                              / int S(lambda) lambda dlambda )
               = -2.5 log10(<f_lambda> / <f_lambda_Vega>)

    i.e. mean_flambda() evaluated for both the target spectrum and a
    reference Vega spectrum, through the same bandpass.

    Parameters
    ----------
    wave_nm, flux_lambda : array
        Target spectrum, same convention as mean_flambda().
    bp_wave_nm, bp_transmission : array
        Bandpass throughput, same convention as mean_flambda().
    vega_wave_nm, vega_flux_lambda : array
        Reference Vega spectrum, in the same flux units as `flux_lambda`
        (e.g. W/m^2/nm), on its own wavelength grid (nm).

    Returns
    -------
    float
        Synthetic Vega-system magnitude.
    """
    f_star = mean_flambda(wave_nm, flux_lambda, bp_wave_nm, bp_transmission)
    f_vega = mean_flambda(vega_wave_nm, vega_flux_lambda, bp_wave_nm, bp_transmission)
    return -2.5 * np.log10(f_star / f_vega)


LN10 = np.log(10.0)


def _trapz_weights(x):
    """
    Quadrature weights w such that sum(w * y) == np.trapezoid(y, x).

    Needed for error propagation: mean_fnu()/mean_flambda() are just
    weighted sums of f_lambda, so once expressed as sum(w_i * f_i) instead
    of np.trapezoid(...), standard linear error propagation applies exactly.
    """
    x = np.asarray(x, dtype=float)
    w = np.zeros_like(x)
    w[0] = (x[1] - x[0]) / 2
    w[-1] = (x[-1] - x[-2]) / 2
    w[1:-1] = (x[2:] - x[:-2]) / 2
    return w


def _overlap_mask(wave_nm, bp_wave_nm):
    lo = max(wave_nm.min(), bp_wave_nm.min())
    hi = min(wave_nm.max(), bp_wave_nm.max())
    mask = (bp_wave_nm >= lo) & (bp_wave_nm <= hi)
    if not mask.any():
        raise ValueError("Spectrum and bandpass wavelength ranges do not overlap")
    return mask


def mean_fnu_with_error(wave_nm, flux_lambda, flux_lambda_err, bp_wave_nm, bp_transmission):
    """
    mean_fnu(), plus its propagated uncertainty.

    ONLY approximate: treats `flux_lambda_err` as the diagonal of the flux
    covariance matrix (i.e. assumes each wavelength bin's error is
    independent). Real Gaia XP flux errors are strongly correlated across
    wavelength (the spectrum is reconstructed from ~55 basis coefficients,
    not independent per-pixel samples), so this is a known *underestimate*
    of the true uncertainty -- getting it right requires the coefficient
    covariance matrix, not just the marginal per-sample errors used here.

    Returns
    -------
    fnu, sigma_fnu : float
        <f_nu> and its propagated 1-sigma uncertainty, both in W/m^2/Hz.
    """
    wave_nm = np.asarray(wave_nm)
    flux_lambda = np.asarray(flux_lambda)
    flux_lambda_err = np.asarray(flux_lambda_err)
    bp_wave_nm = np.asarray(bp_wave_nm)
    bp_transmission = np.asarray(bp_transmission)

    mask = _overlap_mask(wave_nm, bp_wave_nm)
    lam = bp_wave_nm[mask]
    S = bp_transmission[mask]

    f_lam = np.interp(lam, wave_nm, flux_lambda)
    sigma_lam = np.interp(lam, wave_nm, flux_lambda_err)  # linear-interp the errors too (approximate)

    w = _trapz_weights(lam) * S * lam
    numerator = np.sum(w * f_lam)
    sigma_numerator = np.sqrt(np.sum((w * sigma_lam) ** 2))

    denominator = np.trapezoid(S * (C_NM_PER_S / lam), lam)

    fnu = numerator / denominator
    sigma_fnu = sigma_numerator / denominator
    return fnu, sigma_fnu


def mean_flambda_with_error(wave_nm, flux_lambda, flux_lambda_err, bp_wave_nm, bp_transmission):
    """
    mean_flambda(), plus its propagated uncertainty. Same independent-bin
    caveat as mean_fnu_with_error() -- see its docstring.

    Returns
    -------
    flambda, sigma_flambda : float
        <f_lambda> and its propagated 1-sigma uncertainty, both in the same
        units as `flux_lambda` (e.g. W/m^2/nm).
    """
    wave_nm = np.asarray(wave_nm)
    flux_lambda = np.asarray(flux_lambda)
    flux_lambda_err = np.asarray(flux_lambda_err)
    bp_wave_nm = np.asarray(bp_wave_nm)
    bp_transmission = np.asarray(bp_transmission)

    mask = _overlap_mask(wave_nm, bp_wave_nm)
    lam = bp_wave_nm[mask]
    S = bp_transmission[mask]

    f_lam = np.interp(lam, wave_nm, flux_lambda)
    sigma_lam = np.interp(lam, wave_nm, flux_lambda_err)

    w = _trapz_weights(lam) * S * lam
    numerator = np.sum(w * f_lam)
    sigma_numerator = np.sqrt(np.sum((w * sigma_lam) ** 2))

    denominator = np.trapezoid(S * lam, lam)

    flambda = numerator / denominator
    sigma_flambda = sigma_numerator / denominator
    return flambda, sigma_flambda


def ab_mag_with_error(fnu_w_m2_hz, sigma_fnu_w_m2_hz):
    """ab_mag(), plus its propagated 1-sigma uncertainty (both in mag)."""
    mag = ab_mag(fnu_w_m2_hz)
    sigma_mag = (2.5 / LN10) * (sigma_fnu_w_m2_hz / fnu_w_m2_hz)
    return mag, sigma_mag


def vega_mag_with_error(f_star, sigma_f_star, f_vega):
    """
    vega_mag(), plus its propagated 1-sigma uncertainty (both in mag).
    Takes the already-computed mean_flambda() values (star and Vega) rather
    than raw spectra, since the Vega reference spectrum is treated as having
    no error budget of its own here.
    """
    mag = -2.5 * np.log10(f_star / f_vega)
    sigma_mag = (2.5 / LN10) * (sigma_f_star / f_star)
    return mag, sigma_mag


def fnu_to_jy(fnu_w_m2_hz):
    """W/m^2/Hz -> Jy (1 Jy = 1e-26 W/m^2/Hz)."""
    return fnu_w_m2_hz / 1e-26


def ab_mag(fnu_w_m2_hz):
    """W/m^2/Hz -> AB magnitude."""
    return -2.5 * np.log10(fnu_to_jy(fnu_w_m2_hz) / AB_ZEROPOINT_JY)
