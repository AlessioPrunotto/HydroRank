import numpy as np
import pytest

from hydrarank.config import PreprocessConfig
from hydrarank.data import WaterObservations
from hydrarank.exceptions import HydraRankError
from hydrarank.preprocess import prepare_system, run_preprocess


@pytest.fixture
def config(tmp_path):
    return PreprocessConfig(
        topology=tmp_path / "top.psf",
        ligand_selection="resname 547",
        water_cutoff=5.0,
        pocket_cutoff=10.0,
    )


def test_prepare_system_sets_everything_up(pocket_system, config):
    system = prepare_system(pocket_system, config)

    assert system.ligand.n_atoms == 4
    assert system.ligand_heavy.n_atoms == 4
    assert {"O1", "N1", "S1"} <= set(system.ligand_heavy.names)
    assert system.water.n_waters == 2
    assert system.oxygens.n_atoms == 2
    assert system.fit_group.n_atoms == 4
    assert system.reference_frame == 0
    assert system.ligand_reference.shape == (4, 3)
    assert pocket_system.trajectory.transformations


def test_observations_have_one_row_per_water_and_frame(pocket_system, config):
    observations = run_preprocess(pocket_system, config)

    assert observations.n_frames == 4
    assert observations.n_observations == 8  # 2 waters x 4 frames
    assert observations.hydrogen.shape == (8, 2, 3)
    assert set(observations.frame) == {0, 1, 2, 3}
    assert set(observations.resid) == {6, 7}
    assert len(np.unique(observations.water_id)) == 2
    assert np.allclose(observations.fit_rmsd, 0.0, atol=1e-6)


def test_water_identity_remains_unique_when_resids_repeat(pocket_system, config):
    pocket_system.residues[-2:].resids = 99
    observations = run_preprocess(pocket_system, config)

    assert set(observations.resid) == {99}
    assert len(np.unique(observations.water_id)) == 2


def test_virtual_site_never_enters_the_observations(pocket_system, config):
    observations = run_preprocess(pocket_system, config)
    # the OPC virtual site sits on top of the oxygen; hydrogens must be ~0.9 A away
    distances = np.linalg.norm(observations.hydrogen - observations.oxygen[:, None, :], axis=2)
    assert np.all(distances > 0.5)


def test_shell_is_defined_against_the_reference_ligand_pose(pocket_system, config):
    observations = run_preprocess(pocket_system, config)
    distances = np.linalg.norm(
        observations.oxygen[:, None, :] - observations.ligand_reference[None], axis=2
    )
    assert np.all(distances.min(axis=1) <= config.water_cutoff + 1e-6)


def test_tight_cutoff_yields_no_observations(pocket_system, tmp_path):
    config = PreprocessConfig(
        topology=tmp_path / "top.psf",
        ligand_selection="resname 547",
        water_cutoff=0.5,
        pocket_cutoff=10.0,
    )
    observations = run_preprocess(pocket_system, config)
    assert observations.n_observations == 0
    assert observations.oxygen.shape == (0, 3)
    assert observations.n_frames == 4


def test_frame_range_is_honoured(pocket_system, tmp_path):
    config = PreprocessConfig(
        topology=tmp_path / "top.psf",
        ligand_selection="resname 547",
        pocket_cutoff=10.0,
        start=1,
        step=2,
    )
    observations = run_preprocess(pocket_system, config)
    assert list(observations.frames) == [1, 3]
    assert set(observations.frame) <= {0, 1}


def test_progress_callback_receives_each_processed_frame(pocket_system, config):
    updates = []
    run_preprocess(
        pocket_system, config, progress=lambda current, total: updates.append((current, total))
    )
    assert updates == [(1, 4), (2, 4), (3, 4), (4, 4)]


def test_observations_survive_a_roundtrip(pocket_system, config, tmp_path):
    observations = run_preprocess(pocket_system, config)
    loaded = WaterObservations.load(observations.save(tmp_path / "obs.npz"))
    assert np.allclose(loaded.oxygen, observations.oxygen)
    assert loaded.metadata["water_cutoff"] == 5.0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"water_cutoff": 25.0, "pocket_cutoff": 25.0},
        {"pocket_cutoff": 25.0},
        {"hbond_distance": 25.0},
        {"enclosure_radius": 25.0},
    ],
)
def test_every_spatial_cutoff_is_validated_against_the_box(pocket_system, tmp_path, kwargs):
    values = {
        "topology": tmp_path / "top.psf",
        "ligand_selection": "resname 547",
        "pocket_cutoff": 10.0,
    }
    values.update(kwargs)
    config = PreprocessConfig(**values)
    with pytest.raises(HydraRankError, match="minimum-image"):
        prepare_system(pocket_system, config)
