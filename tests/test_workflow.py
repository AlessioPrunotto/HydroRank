import json

import numpy as np

import hydrarank.workflow as workflow
from hydrarank.config import PreprocessConfig
from hydrarank.data import WaterObservations
from hydrarank.io import SystemReport


def _observations(n=20):
    rng = np.random.default_rng(7)
    oxygen = rng.normal(scale=0.08, size=(n, 3))
    local_hydrogen = np.array([[0.586, 0.757, 0.0], [0.586, -0.757, 0.0]])
    return WaterObservations(
        oxygen=oxygen,
        hydrogen=oxygen[:, None, :] + local_hydrogen,
        frame=np.arange(n),
        resid=np.ones(n, dtype=int),
        water_id=np.ones(n, dtype=int),
        frames=np.arange(n),
        times=np.arange(n) * 2.0,
        fit_rmsd=np.full(n, 0.5),
        ligand_reference=np.zeros((1, 3)),
        hb_solute=np.zeros(n),
        enclosure=np.ones(n),
        metadata={
            "water_cutoff": 5.0,
            "reference_frame": 0,
            "n_fit_atoms": 12,
            "qc_longest_bond": [1.5] * n,
            "qc_n_shell_waters": [1] * n,
        },
    )


def _system_report(n_frames=20):
    return SystemReport(
        n_atoms=100,
        n_frames=n_frames,
        dt_ps=2.0,
        n_selected_frames=n_frames,
        box=(40.0, 40.0, 40.0, 90.0, 90.0, 90.0),
        is_orthorhombic=True,
        has_bonds=True,
        residue_counts={"HOH": 20},
        n_waters=20,
        water_model="3-site",
        water_resname="HOH",
        ligand_label="resname LIG",
        n_ligand_atoms=10,
        n_ligand_heavy_atoms=8,
    )


def test_analysis_writes_bundle_and_reuses_compatible_cache(monkeypatch, tmp_path):
    topology = tmp_path / "system.psf"
    trajectory = tmp_path / "trajectory.xtc"
    topology.write_text("topology")
    trajectory.write_text("trajectory")
    config = PreprocessConfig(
        topology=topology,
        trajectory=[trajectory],
        ligand_selection="resname LIG",
        density_factor=0.1,
        output_dir=tmp_path / "results",
    )
    calls = []

    monkeypatch.setattr(workflow, "load_universe", lambda config: object())
    monkeypatch.setattr(
        workflow, "describe_system", lambda universe, config: (_system_report(), None)
    )

    def preprocess(universe, config, progress=None, collect_qc=False):
        calls.append(collect_qc)
        return _observations()

    monkeypatch.setattr(workflow, "run_preprocess", preprocess)
    first = workflow.run_analysis(config)

    assert calls == [True]
    assert first.cache_reused is False
    assert first.ranking.analysis.n_sites == 1
    expected = (
        "observations.npz",
        "config.yaml",
        "ranking.csv",
        "sites.pdb",
        "results.json",
        "report.txt",
    )
    for name in expected:
        assert (config.output_dir / name).is_file()
    payload = json.loads((config.output_dir / "results.json").read_text())
    assert payload["hydrarank_version"]
    assert payload["preprocessing_qc"]["solute_always_whole"] is True

    monkeypatch.setattr(
        workflow,
        "run_preprocess",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("cache not reused")),
    )
    second = workflow.run_analysis(config)
    assert second.cache_reused is True


def test_force_rebuilds_compatible_cache(monkeypatch, tmp_path):
    topology = tmp_path / "system.psf"
    topology.write_text("topology")
    config = PreprocessConfig(
        topology=topology,
        ligand_selection="resname LIG",
        density_factor=0.1,
        output_dir=tmp_path / "results",
    )
    calls = []
    monkeypatch.setattr(workflow, "load_universe", lambda config: object())
    monkeypatch.setattr(
        workflow, "describe_system", lambda universe, config: (_system_report(), None)
    )
    monkeypatch.setattr(
        workflow,
        "run_preprocess",
        lambda universe, config, progress=None, collect_qc=False: (
            calls.append(collect_qc) or _observations()
        ),
    )

    workflow.run_analysis(config)
    rebuilt = workflow.run_analysis(config, force=True)

    assert calls == [True, True]
    assert rebuilt.cache_reused is False
