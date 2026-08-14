"""
Step 1 (optional -- skip if you're using the provided xp_spectra_ecsv):
find which Gaia source_ids need an XP continuous spectrum downloaded.

Matches the b-band SuperBIT star catalog to the merged Gaia DR3 catalog
(config.yaml's gaia_catalog) per target, and keeps the source_ids that both
matched a SuperBIT detection AND have XPcont==1 (i.e. Gaia actually has a
continuous spectrum for them). Writes xp_continuous_source_ids.txt, which
scripts/02_fetch_xp_spectra.py consumes.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from astropy.table import Table
from config import load_config
from skycoord_matcher import SkyCoordMatcher


def gather_xp_ids(cfg):
    mega_cat = Table.read(cfg["starcat_b"])
    gaia_cat = Table.read(cfg["gaia_catalog"])
    targets = sorted(set(mega_cat["TARGET"]))

    all_ids = set()
    per_target_counts = {}

    for target in targets:
        cat_sb = mega_cat[mega_cat["TARGET"] == target]
        gaia_stars = gaia_cat[gaia_cat["TARGET"] == target]

        if len(gaia_stars) == 0:
            print(f"[{target}] WARNING: no Gaia rows for this target, skipping")
            per_target_counts[target] = None
            continue

        matcher = SkyCoordMatcher(
            cat_sb, gaia_stars,
            cat1_ratag="ALPHAWIN_J2000", cat1_dectag="DELTAWIN_J2000",
            cat2_ratag="ALPHAWIN_J2000", cat2_dectag="DELTAWIN_J2000",
            return_idx=True, match_radius=cfg["match_tolerance_deg"], verbose=False,
        )
        _, matched_gaia, _, _ = matcher.get_matched_pairs()

        xp_ids = matched_gaia["Source"][matched_gaia["XPcont"] == 1].tolist()
        per_target_counts[target] = len(xp_ids)
        all_ids.update(xp_ids)

        print(f"[{target}] matched={len(matched_gaia)}  with XPcont: {len(xp_ids)}")

    return all_ids, per_target_counts


if __name__ == "__main__":
    cfg = load_config()
    all_ids, counts = gather_xp_ids(cfg)

    print("\n=== Summary ===")
    for target, n in counts.items():
        print(f"  {target:25s} {n}")
    print(f"\nTotal unique source_ids with continuous spectra: {len(all_ids)}")

    outfile = os.path.join(cfg["external_data_dir"], "xp_continuous_source_ids.txt")
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w") as f:
        for sid in sorted(all_ids):
            f.write(f"{sid}\n")
    print(f"Wrote id list to {outfile}")
