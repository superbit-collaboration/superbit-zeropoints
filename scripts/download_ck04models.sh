#!/bin/bash
# Downloads the Castelli-Kurucz (ck04models) stellar atmosphere grid used
# to extend Gaia XP spectra below their ~336nm reliable floor (~41MB).
#
# Usage: ./scripts/download_ck04models.sh [destination_dir]
#   destination_dir defaults to external_data/cdbs, matching config.yaml's
#   default cdbs_dir. If you change cdbs_dir in config.yaml, pass the
#   matching directory here too.
#
# Source: ssb.stsci.edu (the archive.stsci.edu mirror blocks /hlsps/ in
# robots.txt; this host doesn't, so we use it instead).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${1:-$REPO_ROOT/external_data/cdbs}"

mkdir -p "$DEST/grid"
cd "$DEST/grid"

echo "Downloading ck04models grid to $DEST/grid/ck04models ..."
wget -q -r -np -nH --cut-dirs=2 -A "*.fits,AA_README" \
  "https://ssb.stsci.edu/trds/grid/ck04models/"

n_files=$(find ck04models -type f | wc -l)
echo "Done: $n_files files in $DEST/grid/ck04models"
echo "Set PYSYN_CDBS to $DEST (config.yaml's cdbs_dir already points here by default)."
