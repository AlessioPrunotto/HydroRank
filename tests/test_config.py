import pytest
import yaml

from water_entropy.config import PreprocessConfig
from water_entropy.exceptions import WaterEntropyError


def test_roundtrip_yaml(tmp_path):
    config = PreprocessConfig(
        topology=tmp_path / "top.psf",
        trajectory=[tmp_path / "traj.xtc"],
        ligand_selection="resname 547",
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


def test_unknown_key_rejected(tmp_path):
    with pytest.raises(WaterEntropyError, match="unknown configuration keys"):
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
    with pytest.raises(WaterEntropyError):
        PreprocessConfig(topology=tmp_path / "top.psf", **kwargs)
