"""First-order entropy proxies for the waters in a hydration site.

Both estimators are k-nearest-neighbour (Kozachenko-Leonenko) estimates of the
entropy of a site's water population *relative to bulk water*, which is what matters
for displacement: a water that is already as disordered as bulk gains nothing when
pushed out. Values are per water molecule, in units of the gas constant, and are
negative for ordered sites.

This is the first-order inhomogeneous-solvation-theory approximation used by
WaterMap and SSTMap: translational and orientational contributions are treated
separately and water-water correlations are ignored. It is much cheaper than a full
GIST grid and, for ranking sites against each other, usually good enough.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

#: Number density of bulk water at 300 K, in molecules per cubic angstrom.
BULK_WATER_DENSITY = 0.0329

#: Gas constant in kcal / (mol K).
GAS_CONSTANT = 1.987204259e-3

#: Below this many observations the nearest-neighbour estimator is meaningless.
MIN_SAMPLES = 10


def water_orientations(oxygen: np.ndarray, hydrogen: np.ndarray) -> Rotation:
    """Orientation of each water as a rotation of its molecular frame.

    The frame is the dipole bisector, the in-plane perpendicular, and the plane
    normal. Degenerate geometries yield an identity rotation and are filtered out by
    the caller through the returned validity mask of :func:`orientational_entropy`.
    """
    oh1 = hydrogen[:, 0] - oxygen
    oh2 = hydrogen[:, 1] - oxygen
    bisector = _normalise(oh1 + oh2)
    normal = _normalise(np.cross(oh1, oh2))
    tangent = np.cross(normal, bisector)
    return Rotation.from_matrix(np.stack([bisector, tangent, normal], axis=-1))


def translational_entropy(
    positions: np.ndarray,
    density: float = BULK_WATER_DENSITY,
    min_samples: int = MIN_SAMPLES,
) -> float:
    """Excess translational entropy per water, in units of the gas constant.

    Zero means the water is spread out exactly as in bulk; negative means confined.
    """
    positions = np.asarray(positions, dtype=np.float64)
    if positions.shape[0] < min_samples:
        return float("nan")

    distances, _ = cKDTree(positions).query(positions, k=2)
    nearest = distances[:, 1]
    valid = nearest > 0
    if np.count_nonzero(valid) < min_samples:
        return float("nan")

    n_pairs = np.count_nonzero(valid) - 1
    ball_volume = 4.0 / 3.0 * np.pi * nearest[valid] ** 3
    return float(np.mean(np.log(n_pairs * density * ball_volume)) + np.euler_gamma)


def orientational_entropy(
    orientations: Rotation,
    min_samples: int = MIN_SAMPLES,
) -> float:
    """Excess orientational entropy per water, in units of the gas constant.

    Measured against uniformly random orientations, i.e. bulk water, so the value is
    at most zero. The two hydrogens are treated as indistinguishable: without that
    symmetry every water would look half as ordered as it is.
    """
    quaternions = np.atleast_2d(orientations.as_quat())
    n_samples = quaternions.shape[0]
    if n_samples < min_samples:
        return float("nan")

    # a half turn about the body x axis (the dipole bisector) swaps the two hydrogens
    swapped = (orientations * Rotation.from_euler("x", np.pi)).as_quat()

    overlap = np.maximum(np.abs(quaternions @ quaternions.T), np.abs(quaternions @ swapped.T))
    np.fill_diagonal(overlap, -np.inf)
    nearest = np.arccos(np.clip(overlap.max(axis=1), -1.0, 1.0)) * 2.0

    valid = nearest > 0
    if np.count_nonzero(valid) < min_samples:
        return float("nan")

    theta = nearest[valid]
    # the reference space is SO(3) folded by the hydrogen-swap symmetry, hence the factor 2
    ball_measure = np.minimum(2.0 * (theta - np.sin(theta)) / np.pi, 1.0)
    n_pairs = np.count_nonzero(valid) - 1
    return float(np.mean(np.log(n_pairs * ball_measure)) + np.euler_gamma)


def minus_t_delta_s(entropy: float, temperature: float = 300.0) -> float:
    """Convert an entropy in units of the gas constant into -T dS, in kcal/mol."""
    return -temperature * GAS_CONSTANT * entropy


def _normalise(vectors: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(vectors, axis=-1, keepdims=True)
    return vectors / np.where(norms > 0, norms, 1.0)
