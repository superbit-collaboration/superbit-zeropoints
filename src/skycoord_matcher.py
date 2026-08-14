"""
Standalone copy of superbit_lensing.match.SkyCoordMatcher, so this repo
doesn't require installing the full superbit-lensing package. Nearest-
neighbor sky-position matching between two catalogs via a KD-tree on unit
sphere Cartesian coordinates, with mutual (one-to-one) reassignment so each
catalog-2 object is only claimed by its single closest catalog-1 match.
"""
import numpy as np
from astropy.table import hstack
from scipy.spatial import cKDTree


class SkyCoordMatcher:
    def __init__(self, cat1, cat2,
                 cat1_ratag='ALPHAWIN_J2000', cat1_dectag='DELTAWIN_J2000',
                 cat2_ratag=None, cat2_dectag=None, return_idx=False,
                 match_radius=1.0/3600, verbose=True):
        """
        Initialize the SkyCoordMatcher with catalogs and matching parameters.

        Parameters:
        - cat1, cat2: Input catalogs (structured arrays or Astropy Tables).
        - cat1_ratag, cat1_dectag: Column names for RA and DEC in cat1.
        - cat2_ratag, cat2_dectag: Column names for RA and DEC in cat2 (defaults to same as cat1 if not provided).
        - match_radius: Matching radius in degrees.
        """
        self.cat1 = cat1
        self.cat2 = cat2

        self.cat1_ratag = cat1_ratag
        self.cat1_dectag = cat1_dectag
        self.cat2_ratag = cat2_ratag if cat2_ratag is not None else cat1_ratag
        self.cat2_dectag = cat2_dectag if cat2_dectag is not None else cat1_dectag

        self.match_radius = match_radius  # in degrees
        self.return_idx = return_idx
        self.verbose = verbose

        self.matched_cat1 = None
        self.matched_cat2 = None
        self.matched_cat = None
        self.Nobjs = 0

        self._match()

    def _sph_to_cart(self, ra_rad, dec_rad):
        """Convert spherical coordinates to 3D Cartesian unit vectors."""
        x = np.cos(dec_rad) * np.cos(ra_rad)
        y = np.cos(dec_rad) * np.sin(ra_rad)
        z = np.sin(dec_rad)
        return np.stack((x, y, z), axis=-1)

    def _match(self):
        """
        Perform the matching between the two catalogs using KD-tree.
        """
        # Convert RA/DEC to radians
        ra1 = np.radians(self.cat1[self.cat1_ratag])
        dec1 = np.radians(self.cat1[self.cat1_dectag])
        ra2 = np.radians(self.cat2[self.cat2_ratag])
        dec2 = np.radians(self.cat2[self.cat2_dectag])

        # Convert to unit vectors
        vec1 = self._sph_to_cart(ra1, dec1)
        vec2 = self._sph_to_cart(ra2, dec2)

        # Convert match radius (deg) to chord length on unit sphere
        tolerance_rad = np.radians(self.match_radius)
        chord_tol = 2 * np.sin(tolerance_rad / 2)

        tree = cKDTree(vec2)

        dists, idx = tree.query(vec1, k=1, distance_upper_bound=chord_tol)

        # Handle unmatched points
        valid = (idx < len(self.cat2)) & (dists != np.inf)
        matches = -np.ones(len(vec1), dtype=int)
        matches[valid] = idx[valid]

        # Convert chord length to angular separation
        dists = np.clip(dists, 0, 2)
        angular_dists = 2 * np.arcsin(dists / 2)
        distances = np.full(len(vec1), np.inf)
        distances[valid] = angular_dists[valid]
        if self.verbose:
            print(f"Number of matches after initial pass: {np.sum(valid)}")

        # Enforce one-to-one matching (closest wins)
        best_match_for_2 = {}
        for i, match in enumerate(matches):
            if match != -1:
                if match not in best_match_for_2 or distances[i] < distances[best_match_for_2[match]]:
                    best_match_for_2[match] = i

        final_matches1 = list(best_match_for_2.values())
        final_matches2 = list(best_match_for_2.keys())
        if self.verbose:
            print(f"Number of matches after reassignment: {len(final_matches1)}")

        self.matched_cat1 = self.cat1[final_matches1]
        self.matched_cat2 = self.cat2[final_matches2]
        self.matched_cat2['separation'] = distances[final_matches1]
        self.matched_cat = hstack([self.matched_cat1, self.matched_cat2])
        self.match_idx1 = final_matches1
        self.match_idx2 = final_matches2
        self.Nobjs = len(self.matched_cat)

    def get_matched_catalog(self):
        """Return the matched catalog."""
        return self.matched_cat

    def get_matched_pairs(self):
        """
        Return the matched pairs of objects.

        Returns:
        - matched_cat1, matched_cat2: Matched subsets of the input catalogs.
        """
        if not self.return_idx:
            return self.matched_cat1, self.matched_cat2
        return self.matched_cat1, self.matched_cat2, self.match_idx1, self.match_idx2
