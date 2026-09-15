import pytest

from water_entropy.config import PreprocessConfig
from water_entropy.exceptions import TrajectoryError
from water_entropy.io import describe_system, frame_slice


@pytest.fixture
def config(tmp_path):
    return PreprocessConfig(topology=tmp_path / "top.psf", ligand_selection="resname 547")


def test_describe_system(opc_system, config):
    report, water_topology = describe_system(opc_system, config)
    assert report.n_waters == 5
    assert report.water_resname == "OPC"
    assert report.n_ligand_heavy_atoms == 10
    assert report.is_orthorhombic
    assert water_topology.n_sites == 4
    assert "OPC" in str(report)


def test_report_warns_about_short_trajectory(opc_system, config):
    report, _ = describe_system(opc_system, config)
    assert any("frames selected" in message for message in report.warnings)


def test_frame_slice_respects_bounds(tmp_path):
    config = PreprocessConfig(topology=tmp_path / "top.psf", start=1, stop=100, step=2)
    assert frame_slice(config, 10) == slice(1, 10, 2)


def test_frame_slice_rejects_empty_range(tmp_path):
    config = PreprocessConfig(topology=tmp_path / "top.psf", start=20)
    with pytest.raises(TrajectoryError, match="empty frame range"):
        frame_slice(config, 10)
