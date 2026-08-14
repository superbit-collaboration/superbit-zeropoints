"""
Step 2 (optional -- skip if you're using the provided xp_spectra_ecsv):
download Gaia XP continuous mean spectra for the source_ids found by
01_gather_xp_ids.py, directly from ESA's bulk CDN.

Needs network access to cdn.gea.esac.esa.int. Writes to config.yaml's
xp_spectra_ecsv path.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from config import load_config
from fetch_xp_continuous import download_and_extract

if __name__ == "__main__":
    cfg = load_config()

    id_file = os.path.join(cfg["external_data_dir"], "xp_continuous_source_ids.txt")
    with open(id_file) as f:
        ids = [int(line.strip()) for line in f if line.strip()]
    print(f"Downloading spectra for {len(ids)} source_ids...")

    os.makedirs(os.path.dirname(cfg["xp_spectra_ecsv"]), exist_ok=True)
    download_and_extract(ids, output_path=cfg["xp_spectra_ecsv"])
