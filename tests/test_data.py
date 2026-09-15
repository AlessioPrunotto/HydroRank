import json
from dataclasses import replace

import numpy as np
import pytest

from water_entropy.config import PreprocessConfig
from water_entropy.data import WaterObservations
from water_entropy.exceptions import WaterEntropyError


def _observations(n_obs=6, n_frames=3):
    rng = np.random.default_rng(0)
    return WaterObservations(
        oxygen=rng.normal(size=(n_obs, 3)),
        hydrogen=rng.normal(size=(n_obs, 2, 3)),
        frame=np.arange(n_obs) % n_frames,
        resid=np.arange(n_obs),
        frames=np.arange(n_frames),
        times=np.arange(n_frames) * 5.0,
        fit_rmsd=np.zeros(n_frames),
        ligand_reference=rng.normal(size=(4, 3)),
        metadata={"water_cutoff": 5.0},
    )


def test_roundtrip(tmp_path):
    original = _observations()
    loaded = WaterObservations.load(original.save(tmp_path / "obs.npz"))

    assert np.allclose(loaded.oxygen, original.oxygen)
    assert np.allclose(loaded.hydrogen, original.hydrogen)
    assert loaded.metadata == original.metadata
    assert loaded.n_observations == original.n_observations
    assert loaded.n_frames == original.n_frames
    assert np.array_equal(loaded.water_id, original.water_id)


def test_format_two_cache_remains_readable(tmp_path):
    original = _observations()
    path = tmp_path / "v2.npz"
    np.savez_compressed(
        path,
        oxygen=original.oxygen,
        hydrogen=original.hydrogen,
        frame=original.frame,
        resid=original.resid,
        frames=original.frames,
        times=original.times,
        fit_rmsd=original.fit_rmsd,
        ligand_reference=original.ligand_reference,
        hb_solute=original.hb_solute,
        enclosure=original.enclosure,
        metadata=json.dumps({"format_version": 2}),
    )

    loaded = WaterObservations.load(path)
    assert np.array_equal(loaded.water_id, loaded.resid)


def test_dipole_is_normalised():
    dipole = _observations().dipole
    assert np.allclose(np.linalg.norm(dipole, axis=1), 1.0)


def test_dipole_points_from_oxygen_to_hydrogens():
    observations = WaterObservations(
        oxygen=np.zeros((1, 3)),
        hydrogen=np.array([[[1.0, 1.0, 0.0], [1.0, -1.0, 0.0]]]),
        frame=np.zeros(1, dtype=int),
        resid=np.zeros(1, dtype=int),
        frames=np.zeros(1, dtype=int),
        times=np.zeros(1),
        fit_rmsd=np.zeros(1),
        ligand_reference=np.zeros((1, 3)),
    )
    assert np.allclose(observations.dipole, [[1.0, 0.0, 0.0]])


def test_dt_is_derived_from_times():
    assert _observations().dt_ps == 5.0


def test_shape_mismatch_is_rejected():
    with pytest.raises(WaterEntropyError, match="hydrogen has shape"):
        WaterObservations(
            oxygen=np.zeros((3, 3)),
            hydrogen=np.zeros((2, 2, 3)),
            frame=np.zeros(3, dtype=int),
            resid=np.zeros(3, dtype=int),
            frames=np.zeros(1, dtype=int),
            times=np.zeros(1),
            fit_rmsd=np.zeros(1),
            ligand_reference=np.zeros((1, 3)),
        )


def test_frame_arrays_must_agree():
    with pytest.raises(WaterEntropyError, match="same length"):
        WaterObservations(
            oxygen=np.zeros((1, 3)),
            hydrogen=np.zeros((1, 2, 3)),
            frame=np.zeros(1, dtype=int),
            resid=np.zeros(1, dtype=int),
            frames=np.zeros(2, dtype=int),
            times=np.zeros(1),
            fit_rmsd=np.zeros(1),
            ligand_reference=np.zeros((1, 3)),
        )


def test_oxygen_and_ligand_reference_must_be_three_dimensional():
    observations = _observations()
    with pytest.raises(WaterEntropyError, match="oxygen has shape"):
        replace(observations, oxygen=np.zeros((observations.n_observations, 2)))
    with pytest.raises(WaterEntropyError, match="ligand_reference"):
        replace(observations, ligand_reference=np.zeros((4, 2)))


def test_observation_frames_must_be_in_range():
    observations = _observations()
    frames = observations.frame.copy()
    frames[-1] = observations.n_frames
    with pytest.raises(WaterEntropyError, match="outside"):
        replace(observations, frame=frames)


@pytest.mark.parametrize("times", [np.array([0.0, 5.0, 4.0]), np.array([0.0, 5.0, 11.0])])
def test_trajectory_times_must_be_increasing_and_regular(times):
    with pytest.raises(WaterEntropyError, match="trajectory times"):
        replace(_observations(), times=times)


def test_non_finite_coordinates_are_rejected():
    observations = _observations()
    oxygen = observations.oxygen.copy()
    oxygen[0, 0] = np.nan
    with pytest.raises(WaterEntropyError, match="non-finite"):
        replace(observations, oxygen=oxygen)


def test_frames_and_water_ids_must_be_integers():
    with pytest.raises(WaterEntropyError, match="frame must contain integers"):
        replace(_observations(), frame=np.zeros(6, dtype=float))
    with pytest.raises(WaterEntropyError, match="water_id must contain integers"):
        replace(_observations(), water_id=np.zeros(6, dtype=float))


def test_counts_cannot_be_negative():
    with pytest.raises(WaterEntropyError, match="cannot be negative"):
        replace(_observations(), hb_solute=-np.ones(6))


def test_config_yaml_roundtrip_includes_clustering_options(tmp_path):
    config = PreprocessConfig(
        topology=tmp_path / "top.psf",
        output_dir=tmp_path / "out",
        site_radius=1.2,
        density_factor=3.0,
    )
    config.to_yaml(tmp_path / "c.yaml")
    assert PreprocessConfig.from_yaml(tmp_path / "c.yaml") == config
