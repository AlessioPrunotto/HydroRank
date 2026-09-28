import pytest
import yaml

from hydrarank.config import PreprocessConfig
from hydrarank.exceptions import HydraRankError


def test_roundtrip_yaml(tmp_path):
    config = PreprocessConfig(
        topology=tmp_path / "top.psf",
        trajectory=[tmp_path / "traj.xtc"],
        ligand_selection="resname LIG",
        step=2,
        output_dir=tmp_path / "output",
    )
    path = tmp_path / "config.yaml"
    config.to_yaml(path)
    assert PreprocessConfig.from_yaml(path) == config


def test_relative_paths_resolved_against_config_location(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({"topology": "top.psf", "trajectory": "traj.xtc"}))
    config = PreprocessConfig.from_yaml(path)
    assert config.topology == tmp_path / "top.psf"
    assert config.trajectory == [tmp_path / "traj.xtc"]


def test_apo_reference_paths_roundtrip_and_resolve(tmp_path):
    path = tmp_path / "apo.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "topology": "apo.tpr",
                "trajectory": ["apo.xtc"],
                "reference_structure": "bound.tpr",
                "reference_coordinates": "bound.xtc",
                "reference_ligand_selection": "resname LIG",
                "output_dir": "output",
            }
        )
    )
    config = PreprocessConfig.from_yaml(path)
    assert config.reference_structure == tmp_path / "bound.tpr"
    assert config.reference_coordinates == tmp_path / "bound.xtc"
    written = tmp_path / "roundtrip.yaml"
    config.to_yaml(written)
    assert PreprocessConfig.from_yaml(written) == config


@pytest.mark.parametrize(
    "kwargs",
    [
        {"reference_structure": "bound.pdb"},
        {"reference_ligand_selection": "resname LIG"},
        {"reference_coordinates": "bound.xtc"},
        {
            "reference_structure": "bound.tpr",
            "reference_ligand_selection": "resname LIG",
        },
        {
            "ligand_selection": "resname LIG",
            "reference_structure": "bound.pdb",
            "reference_ligand_selection": "resname LIG",
        },
    ],
)
def test_apo_reference_configuration_must_be_complete_and_exclusive(tmp_path, kwargs):
    with pytest.raises(HydraRankError):
        PreprocessConfig(topology=tmp_path / "apo.tpr", **kwargs)


def test_unknown_key_rejected(tmp_path):
    with pytest.raises(HydraRankError, match="unknown configuration keys"):
        PreprocessConfig.from_dict({"topology": "a.psf", "typo": 1}, base_dir=tmp_path)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"step": 0},
        {"start": -1},
        {"start": 2, "stop": 2},
        {"min_ligand_heavy_atoms": 0},
        {"hbond_penalty": -0.1},
        {"water_cutoff": 0.0},
        {"water_cutoff": 10.0, "pocket_cutoff": 5.0},
    ],
)
def test_invalid_values_rejected(tmp_path, kwargs):
    with pytest.raises(HydraRankError):
        PreprocessConfig(topology=tmp_path / "top.psf", **kwargs)
