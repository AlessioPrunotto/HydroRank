import numpy as np
import pytest

from water_entropy.config import PreprocessConfig
from water_entropy.qc import run_preprocess_qc

N_RESIDUES = 4


@pytest.fixture
def config(tmp_path):
    return PreprocessConfig(
        topology=tmp_path / "top.psf",
        ligand_selection="resname 547",
        water_cutoff=5.0,
        pocket_cutoff=10.0,
    )


def test_qc_on_a_rigid_system(pocket_system, config):
    qc = run_preprocess_qc(pocket_system, config)

    assert qc.frames.size == 4
    assert qc.n_fit_atoms == N_RESIDUES
    assert qc.reference_frame == 0
    assert qc.solute_always_whole
    assert qc.fit_rmsd.max() < 1e-5
    assert np.all(qc.n_shell_waters == 2)
    assert "solute whole" in str(qc)


def test_qc_honours_the_frame_range(pocket_system, tmp_path):
    config = PreprocessConfig(
        topology=tmp_path / "top.psf",
        ligand_selection="resname 547",
        pocket_cutoff=10.0,
        start=1,
        step=2,
    )
    qc = run_preprocess_qc(pocket_system, config)
    assert list(qc.frames) == [1, 3]
    assert qc.reference_frame == 1


def test_qc_serialises_to_json_friendly_types(pocket_system, config):
    payload = run_preprocess_qc(pocket_system, config).to_dict()
    assert payload["n_frames"] == 4
    assert payload["solute_always_whole"] is True
    assert payload["n_shell_waters"]["mean"] == 2.0
