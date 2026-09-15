import numpy as np
import pytest

from hydrarank.clustering import (
    HydrationSites,
    bulk_equivalent_count,
    cluster_hydration_sites,
    cluster_positions,
)
from hydrarank.data import WaterObservations
from hydrarank.exceptions import HydraRankError


def _blobs(centers, per_blob=40, sigma=0.25, seed=0):
    rng = np.random.default_rng(seed)
    points = [rng.normal(center, sigma, size=(per_blob, 3)) for center in centers]
    return np.concatenate(points)


def _observations(positions, n_frames=40):
    n = positions.shape[0]
    return WaterObservations(
        oxygen=positions,
        hydrogen=np.repeat(positions[:, None, :], 2, axis=1) + 0.9,
        frame=np.arange(n) % n_frames,
        resid=np.arange(n),
        frames=np.arange(n_frames),
        times=np.arange(n_frames) * 10.0,
        fit_rmsd=np.zeros(n_frames),
        ligand_reference=np.zeros((1, 3)),
    )


def test_three_blobs_give_three_sites():
    centers = [(0.0, 0.0, 0.0), (6.0, 0.0, 0.0), (0.0, 6.0, 0.0)]
    positions = _blobs(centers)

    found, labels = cluster_positions(positions, radius=1.0, min_count=10)

    assert found.shape[0] == 3
    order = np.argsort(found[:, 0] * 10 + found[:, 1])
    assert np.allclose(found[order], np.array(centers)[[0, 2, 1]], atol=0.2)
    assert np.count_nonzero(labels >= 0) > 0.9 * positions.shape[0]


def test_sites_are_ordered_by_population():
    positions = np.concatenate([_blobs([(0.0, 0.0, 0.0)], per_blob=10), _blobs([(6.0, 0.0, 0.0)])])
    found, labels = cluster_positions(positions, radius=1.0, min_count=5)
    counts = np.bincount(labels[labels >= 0])
    assert counts[0] >= counts[1]
    assert np.allclose(found[0], [6.0, 0.0, 0.0], atol=0.2)


def test_no_point_belongs_to_two_sites():
    positions = _blobs([(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)], per_blob=60, sigma=0.3)
    _, labels = cluster_positions(positions, radius=1.0, min_count=5)
    assert labels.shape == (positions.shape[0],)
    # labels are scalar per point by construction; sites must not overlap in membership
    assert np.count_nonzero(labels == 0) + np.count_nonzero(labels == 1) <= positions.shape[0]


def test_sparse_cloud_yields_no_sites():
    rng = np.random.default_rng(3)
    positions = rng.uniform(-30, 30, size=(200, 3))
    found, labels = cluster_positions(positions, radius=1.0, min_count=10)
    assert found.shape == (0, 3)
    assert np.all(labels == -1)


def test_max_sites_is_respected():
    positions = _blobs([(0.0, 0.0, 0.0), (6.0, 0.0, 0.0), (0.0, 6.0, 0.0)])
    found, _ = cluster_positions(positions, radius=1.0, min_count=5, max_sites=2)
    assert found.shape[0] == 2


def test_empty_input():
    found, labels = cluster_positions(np.empty((0, 3)))
    assert found.shape == (0, 3)
    assert labels.size == 0


@pytest.mark.parametrize(
    "positions,radius",
    [(np.zeros((4, 2)), 1.0), (np.zeros((4, 3)), 0.0)],
)
def test_invalid_input_rejected(positions, radius):
    with pytest.raises(HydraRankError):
        cluster_positions(positions, radius=radius)


def test_bulk_threshold_scales_with_frames_and_volume():
    assert bulk_equivalent_count(100, 1.0) < bulk_equivalent_count(1000, 1.0)
    assert bulk_equivalent_count(1000, 1.0) < bulk_equivalent_count(1000, 2.0)
    assert bulk_equivalent_count(1, 1.0) >= 1


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n_frames": 0, "radius": 1.0},
        {"n_frames": 1, "radius": 0.0},
        {"n_frames": 1, "radius": 1.0, "density_factor": 0.0},
    ],
)
def test_bulk_threshold_rejects_invalid_inputs(kwargs):
    with pytest.raises(HydraRankError):
        bulk_equivalent_count(**kwargs)


@pytest.mark.parametrize("kwargs", [{"min_count": 0}, {"max_sites": 0}])
def test_clustering_rejects_invalid_limits(kwargs):
    with pytest.raises(HydraRankError):
        cluster_positions(np.zeros((2, 3)), **kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"centers": np.zeros((1, 2))},
        {"labels": np.array([1])},
        {"n_frames": 0},
        {"min_count": 0},
    ],
)
def test_hydration_sites_validate_their_data(kwargs):
    values = {
        "centers": np.zeros((1, 3)),
        "labels": np.array([0]),
        "radius": 1.0,
        "n_frames": 1,
        "min_count": 1,
    }
    values.update(kwargs)
    with pytest.raises(HydraRankError):
        HydrationSites(**values)


def test_hydration_sites_report_occupancy_and_spread():
    positions = _blobs([(0.0, 0.0, 0.0), (6.0, 0.0, 0.0)], per_blob=40)
    observations = _observations(positions, n_frames=40)

    sites = cluster_hydration_sites(observations, radius=1.0, min_count=10)

    assert sites.n_sites == 2
    occupancy = sites.frame_occupancy(observations.frame)
    assert np.all(occupancy > 0.9)
    assert np.all(sites.spread(observations.oxygen) < 0.6)
    assert np.allclose(sites.mean_waters(), 1.0, atol=0.1)

    rows = sites.summary(observations)
    assert len(rows) == 2
    assert rows[0]["n_waters"] >= rows[1]["n_waters"]


def test_default_threshold_removes_bulk_like_density():
    rng = np.random.default_rng(7)
    positions = rng.uniform(-10, 10, size=(300, 3))
    observations = _observations(positions, n_frames=300)
    sites = cluster_hydration_sites(observations, radius=1.0, density_factor=2.0)
    assert sites.n_sites == 0
