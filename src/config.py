"""
Loads config.yaml and resolves paths relative to the repo root, so every
other module just calls `load_config()` and gets back a dict of ready-to-use
absolute paths.
"""
import os
import yaml

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))  # .../repo/src
REPO_ROOT = os.path.dirname(_THIS_DIR)                   # .../repo
DEFAULT_CONFIG_PATH = os.path.join(REPO_ROOT, "config.yaml")

_PATH_KEYS = {
    "external_data_dir", "starcat_u", "starcat_b", "starcat_g",
    "xp_spectra_ecsv", "gaia_catalog", "bandpass_dir", "gaia_passband_dir",
    "vega_spectrum_file", "vega_rescaled_file", "cdbs_dir",
}


def load_config(path=None):
    path = path or DEFAULT_CONFIG_PATH
    with open(path) as f:
        cfg = yaml.safe_load(f)

    # resolve simple "${other_key}" substitutions (e.g. external_data_dir
    # used inside starcat_u/starcat_b/.../cdbs_dir)
    for key in list(cfg.keys()):
        value = cfg[key]
        if isinstance(value, str):
            for other_key, other_val in cfg.items():
                if isinstance(other_val, str):
                    value = value.replace(f"${{{other_key}}}", other_val)
            cfg[key] = value

    # resolve relative paths against the directory containing config.yaml
    base_dir = os.path.dirname(os.path.abspath(path))
    for key in _PATH_KEYS:
        if key in cfg and not os.path.isabs(cfg[key]):
            cfg[key] = os.path.normpath(os.path.join(base_dir, cfg[key]))

    return cfg
