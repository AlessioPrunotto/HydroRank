import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from hydrarank.analysis import analyse_sites, format_analysis, residence_stats
from hydrarank.clustering import HydrationSites
from hydrarank.data import WaterObservations
from hydrarank.exceptions import HydraRankError

LOCAL_HYDROGENS = np.array([[0.586, 0.757, 0.0], [0.586, -0.757, 0.0]])


def _site_cloud(center, n, position_sigma, orientation_sigma, rng):
    oxygen = rng.normal(center, position_sigma, size=(n, 3))
    if orientation_sigma is None:
        orientations = Rotation.random(n, rng=rng.integers(1 << 30))
    else:
        orientations = Rotation.from_rotvec(rng.normal(scale=orientation_sigma, size=(n, 3)))
    hydrogen = np.stack([orientations.apply(LOCAL_HYDROGENS[i]) for i in range(2)], axis=1)
    return oxygen, oxygen[:, None, :] + hydrogen


def _two_site_system(n_frames=120, ordered_sigma=0.15, disordered_sigma=0.45):
    rng = np.random.default_rng(0)
    ordered = _site_cloud((0.0, 0.0, 0.0), n_frames, ordered_sigma, 0.15, rng)
    disordered = _site_cloud((6.0, 0.0, 0.0), n_frames, disordered_sigma, None, rng)

    oxygen = np.concatenate([ordered[0], disordered[0]])
    hydrogen = np.concatenate([ordered[1], disordered[1]])
    frame = np.concatenate([np.arange(n_frames), np.arange(n_frames)])
    resid = np.concatenate([np.full(n_frames, 1), np.full(n_frames, 2)])

    observations = WaterObservations(
        oxygen=oxygen,
        hydrogen=hydrogen,
        frame=frame,
        resid=resid,
        frames=np.arange(n_frames),
        times=np.arange(n_frames) * 2.0,
        fit_rmsd=np.zeros(n_frames),
        ligand_reference=np.zeros((1, 3)),
    )
    sites = HydrationSites(
        centers=np.array([[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]]),
        labels=np.concatenate([np.zeros(n_frames, int), np.ones(n_frames, int)]),
        radius=1.0,
        n_frames=n_frames,
        min_count=10,
    )
    return observations, sites


def test_residence_stats_splits_visits_at_gaps():
    frames = np.array([0, 1, 2, 5, 6])
    stats = residence_stats(frames, np.ones(5, dtype=int), n_frames=7)
    assert stats.n_episodes == 2
    assert stats.mean_frames == pytest.approx(2.5)
    assert stats.max_frames == 3
    assert stats.n_censored == 2  # one touches the first frame, one the last


def test_residence_stats_tolerates_a_gap():
    frames = np.array([0, 1, 3, 4])
    assert residence_stats(frames, np.ones(4, dtype=int), 10, max_gap=1).n_episodes == 1
    assert residence_stats(frames, np.ones(4, dtype=int), 10, max_gap=0).n_episodes == 2


def test_residence_stats_separates_water_molecules():
    frames = np.array([0, 1, 2, 0, 1, 2])
    resids = np.array([1, 1, 1, 2, 2, 2])
    stats = residence_stats(frames, resids, n_frames=10)
    assert stats.n_episodes == 2
    assert stats.mean_frames == pytest.approx(3.0)


def test_residence_stats_on_an_empty_site():
    stats = residence_stats(np.empty(0, int), np.empty(0, int), 10)
    assert stats.n_episodes == 0
    assert stats.mean_frames == 0.0


def test_ordered_site_has_lower_entropy_than_the_disordered_one():
    observations, sites = _two_site_system()
    analysis = analyse_sites(observations, sites)

    assert analysis.s_trans[0] < analysis.s_trans[1] < 0
    assert analysis.s_orient[0] < analysis.s_orient[1]
    assert analysis.s_orient[1] == pytest.approx(0.0, abs=0.3)  # random orientations = bulk
    assert analysis.minus_t_delta_s[0] > analysis.minus_t_delta_s[1] > 0


def test_occupancy_and_residence_are_reported_per_site():
    observations, sites = _two_site_system(n_frames=120)
    analysis = analyse_sites(observations, sites)

    assert np.allclose(analysis.occupancy, 1.0)
    assert analysis.residence[0].n_episodes == 1
    assert analysis.residence[0].mean_ps(observations.dt_ps) == pytest.approx(240.0)


def test_rows_and_table_are_consistent():
    observations, sites = _two_site_system()
    analysis = analyse_sites(observations, sites)
    rows = analysis.rows()

    assert [row["site"] for row in rows] == [1, 2]
    assert rows[0]["s_total"] == pytest.approx(analysis.s_trans[0] + analysis.s_orient[0])

    table = format_analysis(analysis)
    assert "Hydration sites" in table
    assert "S_tr" in table and "-TdS" in table
    data_rows = [line.split() for line in table.splitlines()]
    assert "encl" in table and "hb" in table
    assert sum(1 for cells in data_rows if len(cells) == 12 and cells[0].isdigit()) == 2


def test_small_sites_yield_nan_rather_than_a_fake_number():
    observations, sites = _two_site_system(n_frames=120)
    sites.labels[:115] = -1  # leave only 5 waters in site 0
    analysis = analyse_sites(observations, sites)

    assert np.isnan(analysis.s_trans[0])
    assert np.isnan(analysis.s_orient[0])
    assert not np.isnan(analysis.s_trans[1])


def test_coarse_sampling_is_flagged_in_the_table():
    observations, sites = _two_site_system()
    observations.times = np.arange(observations.n_frames) * 100.0
    table = format_analysis(analyse_sites(observations, sites))
    assert "residence times are unresolved" in table


def test_analysis_requires_matching_observations_and_sites():
    observations, sites = _two_site_system()
    sites.labels = sites.labels[:-1]
    with pytest.raises(HydraRankError, match="one entry per"):
        analyse_sites(observations, sites)

    observations, sites = _two_site_system()
    sites.n_frames += 1
    with pytest.raises(HydraRankError, match="frame count"):
        analyse_sites(observations, sites)
