import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from water_entropy.entropy import (
    BULK_WATER_DENSITY,
    GAS_CONSTANT,
    minus_t_delta_s,
    orientational_entropy,
    translational_entropy,
    water_orientations,
)


def _water_geometries(orientations, rng):
    """Rigid water molecules with the given orientations, at random positions."""
    n = len(orientations)
    oxygen = rng.normal(scale=5.0, size=(n, 3))
    local = np.array([[0.586, 0.757, 0.0], [0.586, -0.757, 0.0]])
    hydrogen = np.stack([orientations.apply(local[i]) for i in range(2)], axis=1)
    return oxygen, oxygen[:, None, :] + hydrogen


def test_translational_entropy_is_zero_at_bulk_density():
    rng = np.random.default_rng(0)
    side = (1.0 / BULK_WATER_DENSITY) ** (1 / 3)
    positions = rng.uniform(0, side, size=(4000, 3))
    assert translational_entropy(positions) == pytest.approx(0.0, abs=0.1)


def test_translational_entropy_follows_the_volume_scaling_law():
    rng = np.random.default_rng(1)
    positions = rng.normal(scale=0.4, size=(2000, 3))
    wide = translational_entropy(positions)
    narrow = translational_entropy(positions * 0.5)
    # halving every length must cost exactly 3 ln 2 per water
    assert wide - narrow == pytest.approx(3 * np.log(2), abs=0.02)


def test_tighter_sites_are_more_ordered():
    rng = np.random.default_rng(2)
    loose = translational_entropy(rng.normal(scale=0.6, size=(500, 3)))
    tight = translational_entropy(rng.normal(scale=0.2, size=(500, 3)))
    assert tight < loose < 0


def test_translational_entropy_needs_enough_samples():
    rng = np.random.default_rng(3)
    assert np.isnan(translational_entropy(rng.normal(size=(5, 3))))


def test_duplicate_positions_do_not_blow_up():
    positions = np.zeros((40, 3))
    assert np.isnan(translational_entropy(positions))


def test_orientational_entropy_is_zero_for_random_orientations():
    orientations = Rotation.random(4000, rng=1)
    assert orientational_entropy(orientations) == pytest.approx(0.0, abs=0.1)


def test_orientational_entropy_is_negative_for_aligned_waters():
    rng = np.random.default_rng(4)
    loose = Rotation.from_rotvec(rng.normal(scale=0.4, size=(2000, 3)))
    tight = Rotation.from_rotvec(rng.normal(scale=0.2, size=(2000, 3)))
    assert orientational_entropy(tight) < orientational_entropy(loose) < 0


def test_orientational_entropy_ignores_hydrogen_labelling():
    rng = np.random.default_rng(5)
    orientations = Rotation.from_rotvec(rng.normal(scale=0.3, size=(400, 3)))
    oxygen, hydrogen = _water_geometries(orientations, rng)

    reference = orientational_entropy(water_orientations(oxygen, hydrogen))
    swapped = orientational_entropy(water_orientations(oxygen, hydrogen[:, ::-1]))

    half = hydrogen.copy()
    half[::2] = half[::2, ::-1]
    mixed = orientational_entropy(water_orientations(oxygen, half))

    assert swapped == pytest.approx(reference, abs=1e-9)
    assert mixed == pytest.approx(reference, abs=1e-9)


def test_orientational_entropy_needs_enough_samples():
    assert np.isnan(orientational_entropy(Rotation.random(4, rng=0)))


def test_water_orientations_build_an_orthonormal_frame():
    rng = np.random.default_rng(6)
    orientations = Rotation.random(50, rng=7)
    oxygen, hydrogen = _water_geometries(orientations, rng)
    matrices = water_orientations(oxygen, hydrogen).as_matrix()

    identity = matrices @ np.transpose(matrices, (0, 2, 1))
    assert np.allclose(identity, np.eye(3), atol=1e-6)
    assert np.allclose(np.linalg.det(matrices), 1.0, atol=1e-6)


def test_minus_t_delta_s_converts_to_kcal_per_mole():
    assert minus_t_delta_s(-1.0, 300.0) == pytest.approx(300.0 * GAS_CONSTANT)
    assert minus_t_delta_s(0.0) == 0.0
