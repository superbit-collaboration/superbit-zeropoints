"""
Download Gaia DR3 XP continuous mean spectra for specific source_ids directly
from ESA's static bulk CDN, bypassing the archive's data-server/datalink API.

Useful when the Gaia archive's on-demand retrieval service (data-server,
used internally by astroquery's GaiaClass.load_data / gaiaxpy's calibrate()
with a list input) is unavailable, since this CDN is a separate system.

The output .ecsv file can be passed straight to gaiaxpy.calibrate().
"""
import gzip
import re
import sys
import urllib.request

BASE_URL = 'https://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum/'
MD5SUM_URL = BASE_URL + '_MD5SUM.txt'
CHUNK_RE = re.compile(r'XpContinuousMeanSpectrum_(\d+)-(\d+)\.csv\.gz')


def source_id_to_healpix8(source_id: int) -> int:
    """Gaia source_id encodes a level-12 HEALPix index in its top bits;
    the bulk CDN files are chunked by level-8 HEALPix ranges (4**4 coarser)."""
    return source_id >> 43


def get_chunk_ranges():
    with urllib.request.urlopen(MD5SUM_URL) as resp:
        text = resp.read().decode()
    ranges = []
    for line in text.splitlines():
        m = CHUNK_RE.search(line)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            ranges.append((lo, hi, m.group(0)))
    ranges.sort()
    return ranges


def find_chunk_for_source(source_id: int, ranges):
    pix = source_id_to_healpix8(source_id)
    for lo, hi, filename in ranges:
        if lo <= pix <= hi:
            return filename
    raise ValueError(f'No chunk file found for source_id={source_id} (healpix8={pix})')


def download_and_extract(source_ids, output_path='xp_continuous_selected.ecsv'):
    ranges = get_chunk_ranges()
    remaining = set(str(sid) for sid in source_ids)
    needed_chunks = {}
    for sid in source_ids:
        filename = find_chunk_for_source(sid, ranges)
        needed_chunks.setdefault(filename, []).append(str(sid))

    header_lines = None
    matched_rows = []

    for filename, ids_in_chunk in needed_chunks.items():
        url = BASE_URL + filename
        print(f'Downloading {filename} (needed for {len(ids_in_chunk)} source(s))...', file=sys.stderr)
        wanted = set(ids_in_chunk)
        with urllib.request.urlopen(url) as resp:
            with gzip.GzipFile(fileobj=resp) as gz:
                collecting_header = header_lines is None
                local_header = []
                for raw_line in gz:
                    line = raw_line.decode()
                    if collecting_header:
                        local_header.append(line)
                        if line.startswith('source_id'):
                            collecting_header = False
                            if header_lines is None:
                                header_lines = local_header
                        continue
                    sid_str = line.split(',', 1)[0]
                    if sid_str in wanted:
                        matched_rows.append(line)
                        remaining.discard(sid_str)
                        wanted.discard(sid_str)
                        if not wanted:
                            break  # rows are sorted by source_id; nothing more to find in this chunk

    if remaining:
        print(f'Warning: source_id(s) not found in their expected chunk: {sorted(remaining)}', file=sys.stderr)

    with open(output_path, 'w') as f:
        f.writelines(header_lines)
        f.writelines(matched_rows)

    print(f'Wrote {len(matched_rows)} spectra to {output_path}', file=sys.stderr)
    return output_path


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python fetch_xp_continuous.py SOURCE_ID [SOURCE_ID ...]', file=sys.stderr)
        sys.exit(1)
    ids = [int(x) for x in sys.argv[1:]]
    out_path = download_and_extract(ids)
    print(out_path)
