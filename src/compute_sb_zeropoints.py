"""
Compute SuperBIT u, b, g magnitude zeropoints using stars with a downloaded
Gaia XP continuous spectrum.

For each band, independently (since each band's SB catalog has a different
number of detections):
  1. Load that band's SB star catalog (config.yaml's starcat_u/b/g)
  2. Match each target's detections to the merged Gaia catalog by sky
     position (config.yaml's gaia_catalog)
  3. Keep only matches with a downloaded XP continuous spectrum
     (config.yaml's xp_spectra_ecsv)
  4. Compute synthetic AB and Vega magnitudes (+ their propagated errors)
     through the SuperBIT bandpass (extending the spectrum below 336 nm
     where stellar params allow -- see xp_uv_extrapolation.py)
  5. Per-star zeropoint: ZP_i = m_synth,i + 2.5*log10(FLUX_COL_i)
     Per-star zeropoint error: sigma_ZP,i = sqrt(sigma_synth,i^2 + sigma_inst,i^2)
       - sigma_synth,i: propagated from the SED's flux_error (see
         synthetic_photometry.py's *_with_error functions -- a known
         underestimate, since it ignores cross-wavelength correlations in
         the Gaia XP reconstruction, and the UV-extended part of the
         spectrum has no error term at all)
       - sigma_inst,i = (2.5/ln10) * FLUXERR_COL_i / FLUX_COL_i, from
         SuperBIT's own SExtractor photometry -- an independent noise
         source, so combined in quadrature with sigma_synth,i
  6. Aggregate two ways (config.yaml's `method` picks which is primary,
     both always printed):
       - "simple": flat sigma-clipped median/std, no error propagation
       - "weighted": inverse-variance weighted mean using the propagated
         per-star errors, with chi2-rescaling so the sigma-clip and final
         uncertainty are self-consistent
     Then a leave-one-target-out jackknife over both -- the recommended
     number for reporting, since stars within the same target likely share
     correlated systematics, so treating ~4000 stars as independent (as
     both aggregations above do) understates the true uncertainty.

Writes outputs/sb_zeropoints_<band>_<FLUX_COL>.fits (per-star results) and
outputs/sb_zeropoints_<FLUX_COL>.png (diagnostic plot).
"""
import os
import pickle
import time
import numpy as np
import pandas as pd
from astropy.table import Table, vstack
from astropy.stats import sigma_clip
from gaiaxpy import calibrate

from config import load_config, REPO_ROOT
from skycoord_matcher import SkyCoordMatcher
from xp_uv_extrapolation import build_stellar_param_lookup, extend_spectrum_uv
from synthetic_photometry import (
    mean_fnu, ab_mag, mean_flambda, vega_mag,
    mean_fnu_with_error, ab_mag_with_error,
    mean_flambda_with_error, vega_mag_with_error,
)

LN10 = np.log(10.0)
OUTPUT_DIR = os.path.join(REPO_ROOT, "outputs")


def build_extended_spectra_cache(cfg, rebuild=False):
    """
    source_id -> (wave_nm, flux_Wm2nm) for every star with a downloaded XP
    spectrum, extended below 336 nm where Teff/logg/[Fe/H] allow. Cached to
    disk (under external_data_dir/cache/) since the extension step
    (stsynphot grid interpolation) is slow (~25 min for ~4500 stars).
    """
    cache_file = os.path.join(cfg["external_data_dir"], "cache", "extended_spectra_cache.pkl")
    if not rebuild and os.path.exists(cache_file):
        with open(cache_file, "rb") as f:
            cache = pickle.load(f)
        print(f"Loaded cached spectra for {len(cache)} stars from {cache_file}")
        return cache

    stellar_params = build_stellar_param_lookup(cfg["gaia_catalog"])

    sampling = np.arange(336, 1050, 1)
    spectra, sampling = calibrate(
        cfg["xp_spectra_ecsv"], sampling=sampling, save_file=False
    )

    cache = {}
    t0 = time.time()
    n = len(spectra)
    for i, row in spectra.iterrows():
        sid = int(row["source_id"])
        wave, flux = sampling, row["flux"]
        if sid in stellar_params:
            teff, feh, logg = stellar_params[sid]
            try:
                wave, flux = extend_spectrum_uv(
                    flux, sampling, teff, feh, logg,
                    uv_min=cfg["uv_extension_min_nm"],
                    overlap=tuple(cfg["uv_extension_overlap_nm"]),
                )
            except Exception as e:
                print(f"[WARNING] extension failed for {sid}: {e}")
        cache[sid] = (wave, flux)

        if (i + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (n - i - 1) / rate
            print(f"  {i+1}/{n}  ({rate:.1f} stars/s, ETA {eta/60:.1f} min)")

    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    with open(cache_file, "wb") as f:
        pickle.dump(cache, f)
    print(f"Cached {len(cache)} extended spectra to {cache_file}")
    return cache


def build_extended_spectra_err_cache(cfg, spectra_cache, rebuild=False):
    """
    source_id -> (wave_nm, flux_Wm2nm, flux_err_Wm2nm), reusing the already-
    extended (wave, flux) from `spectra_cache` (avoids redoing the slow
    stsynphot extension) and re-deriving flux_error from a fresh (fast)
    calibrate() call.

    The UV-extended portion of the spectrum (below 336 nm, from ck04models)
    has no formal error -- it's assigned flux_err=0 there, meaning
    sigma_synth will be an underestimate specifically for bands whose
    integral leans on that extrapolated region (most notably u).
    """
    cache_file = os.path.join(cfg["external_data_dir"], "cache", "extended_spectra_err_cache.pkl")
    if not rebuild and os.path.exists(cache_file):
        with open(cache_file, "rb") as f:
            cache = pickle.load(f)
        print(f"Loaded cached spectra+errors for {len(cache)} stars from {cache_file}")
        return cache

    sampling = np.arange(336, 1050, 1)
    spectra, sampling = calibrate(
        cfg["xp_spectra_ecsv"], sampling=sampling, save_file=False
    )
    err_lookup = {
        int(row["source_id"]): np.asarray(row["flux_error"])
        for _, row in spectra.iterrows()
    }

    cache = {}
    for sid, (wave, flux) in spectra_cache.items():
        real_err = err_lookup.get(sid)
        if real_err is None:
            continue
        n_extra = len(wave) - len(real_err)
        flux_err = np.concatenate([np.zeros(n_extra), real_err])
        cache[sid] = (wave, flux, flux_err)

    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    with open(cache_file, "wb") as f:
        pickle.dump(cache, f)
    print(f"Cached {len(cache)} spectra+errors to {cache_file}")
    return cache


def match_band_to_gaia(cfg, band, spectra_cache):
    """
    For one SuperBIT band, match its SB catalog to the merged Gaia catalog
    per target, and keep only matches with a downloaded XP spectrum.

    Returns a Table with columns: source_id, TARGET, <FLUX_COL>, plus the
    matched SB row's other columns.
    """
    sb_cat = Table.read(cfg[f"starcat_{band}"])
    gaia_cat = Table.read(cfg["gaia_catalog"])
    targets = sorted(set(sb_cat["TARGET"]))

    matched_rows = []
    for target in targets:
        sb_sub = sb_cat[sb_cat["TARGET"] == target]
        gaia = gaia_cat[gaia_cat["TARGET"] == target]
        if len(gaia) == 0:
            continue

        matcher = SkyCoordMatcher(
            sb_sub, gaia,
            cat1_ratag="ALPHAWIN_J2000", cat1_dectag="DELTAWIN_J2000",
            cat2_ratag="ALPHAWIN_J2000", cat2_dectag="DELTAWIN_J2000",
            return_idx=True, match_radius=cfg["match_tolerance_deg"], verbose=False,
        )
        matched_sb, matched_gaia, _, _ = matcher.get_matched_pairs()
        if len(matched_sb) == 0:
            continue

        source_ids = np.asarray(matched_gaia["Source"], dtype=np.int64)
        in_cache = np.array([sid in spectra_cache for sid in source_ids])
        if not in_cache.any():
            continue

        sub = matched_sb[in_cache]
        sub["source_id"] = source_ids[in_cache]
        matched_rows.append(sub)

    if not matched_rows:
        raise RuntimeError(f"No matches with downloaded spectra for band {band}")

    return vstack(matched_rows, join_type="outer", metadata_conflicts="silent")


def compute_band_zeropoints(cfg, band, spectra_err_cache, vega_wave, vega_flux):
    flux_col, fluxerr_col = cfg["flux_col"], cfg["fluxerr_col"]

    bp = pd.read_csv(os.path.join(cfg["bandpass_dir"], f"{band}_2023.csv"))
    bp_wave, bp_S = bp["wavelengths"].values, bp["transmission"].values

    f_vega_ref = mean_flambda(vega_wave, vega_flux, bp_wave, bp_S)

    matched = match_band_to_gaia(cfg, band, spectra_err_cache)
    print(f"[{band}] {len(matched)} SB detections matched to cached spectra")

    flux = np.asarray(matched[flux_col], dtype=float)
    fluxerr = np.asarray(matched[fluxerr_col], dtype=float)
    good = np.isfinite(flux) & (flux > 0) & np.isfinite(fluxerr) & (fluxerr > 0)

    m_ab = np.full(len(matched), np.nan)
    m_vega = np.full(len(matched), np.nan)
    sigma_synth_ab = np.full(len(matched), np.nan)
    sigma_synth_vega = np.full(len(matched), np.nan)

    for i in np.where(good)[0]:
        sid = int(matched["source_id"][i])
        wave, flx, flx_err = spectra_err_cache[sid]

        fnu, sigma_fnu = mean_fnu_with_error(wave, flx, flx_err, bp_wave, bp_S)
        m_ab[i], sigma_synth_ab[i] = ab_mag_with_error(fnu, sigma_fnu)

        f_star, sigma_f_star = mean_flambda_with_error(wave, flx, flx_err, bp_wave, bp_S)
        m_vega[i], sigma_synth_vega[i] = vega_mag_with_error(f_star, sigma_f_star, f_vega_ref)

    inst_mag = -2.5 * np.log10(np.where(good, flux, np.nan))
    sigma_inst = (2.5 / LN10) * np.where(good, fluxerr / flux, np.nan)

    zp_ab = m_ab + 2.5 * np.log10(np.where(good, flux, np.nan))
    zp_vega = m_vega + 2.5 * np.log10(np.where(good, flux, np.nan))

    sigma_zp_ab = np.sqrt(sigma_synth_ab**2 + sigma_inst**2)
    sigma_zp_vega = np.sqrt(sigma_synth_vega**2 + sigma_inst**2)

    out = Table({
        "source_id": matched["source_id"],
        "TARGET": matched["TARGET"],
        flux_col: flux,
        fluxerr_col: fluxerr,
        "inst_mag": inst_mag,
        "sigma_inst": sigma_inst,
        "synth_AB": m_ab,
        "sigma_synth_AB": sigma_synth_ab,
        "synth_VEGA": m_vega,
        "sigma_synth_VEGA": sigma_synth_vega,
        "zp_AB": zp_ab,
        "sigma_zp_AB": sigma_zp_ab,
        "zp_VEGA": zp_vega,
        "sigma_zp_VEGA": sigma_zp_vega,
    })
    return out


def robust_stats(x):
    """Flat sigma-clipped median/std -- kept for comparison with the
    inverse-variance-weighted result."""
    x = x[np.isfinite(x)]
    clipped = sigma_clip(x, sigma=4, maxiters=5)
    return np.ma.median(clipped), np.ma.std(clipped), clipped.count()


def weighted_zeropoint(zp, sigma_zp, sigma_clip_thresh=4.0, max_iter=15):
    """
    Inverse-variance-weighted mean zeropoint + formal uncertainty.

    sigma_zp is a known underestimate (ignores XP cross-wavelength
    correlations and the UV-extension's missing error term), so clipping
    directly in units of the raw sigma_zp is wrong -- normal star-to-star
    scatter (~0.03-0.1 mag) is a huge number of "sigma" relative to a
    per-star error of ~0.0003 mag, so nearly everything looks like an
    outlier. Instead, self-consistently inflate sigma_zp by a single global
    factor that drives the reduced chi^2 of the kept points to ~1 (the
    standard "rescale errors until chi2/dof=1" approach), and use that
    inflated sigma for both the sigma-clip threshold and the final weights.

    Returns
    -------
    zp_weighted, sigma_zp_weighted, n_used, reduced_chi2, inflation
        sigma_zp_weighted already includes the inflation factor. `inflation`
        is how many times larger the true per-star scatter is than the
        propagated sigma_zp suggested -- report it, since it's a direct
        measure of how much the analytic error propagation was missing.
    """
    mask = np.isfinite(zp) & np.isfinite(sigma_zp) & (sigma_zp > 0)
    zp, sigma_zp = zp[mask], sigma_zp[mask]

    keep = np.ones(len(zp), dtype=bool)
    inflation = 1.0
    for _ in range(max_iter):
        eff_sigma = sigma_zp * inflation
        w = 1.0 / eff_sigma[keep] ** 2
        zp_w = np.sum(w * zp[keep]) / np.sum(w)

        dof = keep.sum() - 1
        if dof < 1:
            break
        reduced_chi2 = np.sum(((zp[keep] - zp_w) / eff_sigma[keep]) ** 2) / dof
        # let inflation shrink as well as grow -- an outlier-driven overshoot
        # on an early iteration must be able to self-correct once that
        # outlier is clipped out, or inflation gets stuck too large forever
        inflation *= np.sqrt(reduced_chi2)

        new_keep = np.abs(zp - zp_w) / (sigma_zp * inflation) < sigma_clip_thresh
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep

    eff_sigma = sigma_zp * inflation
    w = 1.0 / eff_sigma[keep] ** 2
    zp_weighted = np.sum(w * zp[keep]) / np.sum(w)
    sigma_zp_weighted = 1.0 / np.sqrt(np.sum(w))

    reduced_chi2_final = np.sum(((zp[keep] - zp_weighted) / eff_sigma[keep]) ** 2) / (keep.sum() - 1)

    return zp_weighted, sigma_zp_weighted, int(keep.sum()), reduced_chi2_final, inflation


def aggregate_zeropoint(zp, sigma_zp, method):
    """
    Dispatch to either aggregation method, with a common return signature
    so callers (jackknife, plotting, printing) don't need to care which one
    was used.

    Returns
    -------
    value, error, n_used, info
        `info` is a dict of extra diagnostics: empty for "simple";
        {"reduced_chi2", "inflation"} for "weighted".
    """
    if method == "simple":
        value, error, n_used = robust_stats(zp)
        return value, error, n_used, {}
    elif method == "weighted":
        value, error, n_used, chi2, inflation = weighted_zeropoint(zp, sigma_zp)
        return value, error, n_used, {"reduced_chi2": chi2, "inflation": inflation}
    else:
        raise ValueError(f"Unknown method: {method!r} (expected 'simple' or 'weighted')")


def jackknife_zeropoint(zp, sigma_zp, target_labels, method):
    """
    Leave-one-target-out jackknife uncertainty on the aggregate zeropoint.

    Both robust_stats() and weighted_zeropoint() implicitly treat every
    star as an independent measurement, giving an uncertainty that shrinks
    as ~1/sqrt(N_stars). But stars within the same target field likely
    share correlated systematics (exposure conditions, extinction,
    residual calibration effects) -- exactly the kind of scatter driving
    the multi-x error-inflation factors seen per band. If that's true, the
    true number of independent "trials" is closer to N_targets (~30), not
    N_stars (~4000), and this is the standard way to get an uncertainty
    that reflects that: recompute the aggregate zeropoint N_targets times,
    each time leaving one target's stars out entirely, and use the spread
    of those leave-one-out estimates as the error bar.

    Returns
    -------
    theta_full, jackknife_se, n_targets, loo_estimates
        theta_full: the point estimate using every target (report this as
        the central value); jackknife_se: the recommended uncertainty.
    """
    targets = np.asarray(target_labels)
    mask = np.isfinite(zp) & np.isfinite(sigma_zp) & (sigma_zp > 0)
    zp, sigma_zp, targets = zp[mask], sigma_zp[mask], targets[mask]

    theta_full, _, _, _ = aggregate_zeropoint(zp, sigma_zp, method=method)

    unique_targets = sorted(set(targets))
    n = len(unique_targets)
    loo = np.zeros(n)
    for i, t in enumerate(unique_targets):
        keep = targets != t
        loo[i], _, _, _ = aggregate_zeropoint(zp[keep], sigma_zp[keep], method=method)

    theta_bar = np.mean(loo)
    jackknife_se = np.sqrt((n - 1) / n * np.sum((loo - theta_bar) ** 2))

    return theta_full, jackknife_se, n, loo


def main():
    cfg = load_config()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    spectra_cache = build_extended_spectra_cache(cfg)
    spectra_err_cache = build_extended_spectra_err_cache(cfg, spectra_cache)

    vega_data = np.load(cfg["vega_rescaled_file"])
    vega_wave, vega_flux = vega_data["wave_nm"], vega_data["flux"]

    results = {}
    for band in cfg["bands"]:
        out = compute_band_zeropoints(cfg, band, spectra_err_cache, vega_wave, vega_flux)
        out.write(os.path.join(OUTPUT_DIR, f"sb_zeropoints_{band}_{cfg['flux_col']}.fits"), overwrite=True)
        results[band] = out

        med_ab, std_ab, n_ab = robust_stats(np.asarray(out["zp_AB"]))
        med_vega, std_vega, n_vega = robust_stats(np.asarray(out["zp_VEGA"]))
        print(f"[{band}] flat median   ZP_AB   = {med_ab:.4f} +/- {std_ab:.4f}  (N={n_ab})")
        print(f"[{band}] flat median   ZP_VEGA = {med_vega:.4f} +/- {std_vega:.4f}  (N={n_vega})")

        zp_ab_w, sig_ab_w, n_ab_w, chi2_ab, infl_ab = weighted_zeropoint(
            np.asarray(out["zp_AB"]), np.asarray(out["sigma_zp_AB"])
        )
        zp_vega_w, sig_vega_w, n_vega_w, chi2_vega, infl_vega = weighted_zeropoint(
            np.asarray(out["zp_VEGA"]), np.asarray(out["sigma_zp_VEGA"])
        )
        print(f"[{band}] weighted mean ZP_AB   = {zp_ab_w:.4f} +/- {sig_ab_w:.4f}  "
              f"(N={n_ab_w}, reduced chi2={chi2_ab:.2f}, error inflation={infl_ab:.1f}x)")
        print(f"[{band}] weighted mean ZP_VEGA = {zp_vega_w:.4f} +/- {sig_vega_w:.4f}  "
              f"(N={n_vega_w}, reduced chi2={chi2_vega:.2f}, error inflation={infl_vega:.1f}x)")

    print("\n=== Final (jackknife) zeropoints: simple vs weighted ===")
    for band in cfg["bands"]:
        out = results[band]
        target_labels = np.asarray(out["TARGET"])
        print(f"[{band}]")
        for method in ("simple", "weighted"):
            theta_ab, jack_se_ab, n_targets, _ = jackknife_zeropoint(
                np.asarray(out["zp_AB"]), np.asarray(out["sigma_zp_AB"]), target_labels, method=method
            )
            theta_vega, jack_se_vega, _, _ = jackknife_zeropoint(
                np.asarray(out["zp_VEGA"]), np.asarray(out["sigma_zp_VEGA"]), target_labels, method=method
            )
            print(f"  {method:8s} ZP_AB   = {theta_ab:.4f} +/- {jack_se_ab:.4f}  (N_targets={n_targets})")
            print(f"  {method:8s} ZP_VEGA = {theta_vega:.4f} +/- {jack_se_vega:.4f}  (N_targets={n_targets})")

    make_diagnostic_plot(cfg, results)


def make_diagnostic_plot(cfg, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["text.usetex"] = False  # gaiaxpy forces this True on import

    method = cfg["method"]
    bands = cfg["bands"]

    fig, axes = plt.subplots(1, len(bands), figsize=(16, 5), sharey=False)
    for ax, band in zip(axes, bands):
        out = results[band]
        inst_mag = np.asarray(out["inst_mag"])
        zp_ab = np.asarray(out["zp_AB"])
        zp_vega = np.asarray(out["zp_VEGA"])
        sigma_zp_ab = np.asarray(out["sigma_zp_AB"])
        sigma_zp_vega = np.asarray(out["sigma_zp_VEGA"])
        target_labels = np.asarray(out["TARGET"])

        zp_ab_w, _, n_ab_w, _ = aggregate_zeropoint(zp_ab, sigma_zp_ab, method=method)
        zp_vega_w, _, n_vega_w, _ = aggregate_zeropoint(zp_vega, sigma_zp_vega, method=method)

        theta_ab, jack_se_ab, _, _ = jackknife_zeropoint(zp_ab, sigma_zp_ab, target_labels, method=method)
        theta_vega, jack_se_vega, _, _ = jackknife_zeropoint(zp_vega, sigma_zp_vega, target_labels, method=method)

        ax.scatter(inst_mag, zp_ab, s=6, alpha=0.35, color="tab:blue", label=f"AB (N={n_ab_w})")
        ax.scatter(inst_mag, zp_vega, s=6, alpha=0.35, color="tab:orange", label=f"Vega (N={n_vega_w})")
        ax.axhline(zp_ab_w, color="tab:blue", lw=1.2, ls="--")
        ax.axhline(zp_vega_w, color="tab:orange", lw=1.2, ls="--")

        ax.set_xlabel(f"Instrumental mag (-2.5log10({cfg['flux_col']}))")
        ax.set_ylabel("Per-star zeropoint")
        ax.set_title(
            f"{band}-band\n"
            f"AB: {theta_ab:.3f} +/- {jack_se_ab:.3f}   "
            f"Vega: {theta_vega:.3f} +/- {jack_se_vega:.3f}",
            fontsize=10,
        )
        ax.legend(fontsize=8, loc="best")

    plt.tight_layout()
    outname = os.path.join(OUTPUT_DIR, f"sb_zeropoints_{cfg['flux_col']}.png")
    plt.savefig(outname, dpi=150)
    print(f"Wrote {outname}")


if __name__ == "__main__":
    main()
