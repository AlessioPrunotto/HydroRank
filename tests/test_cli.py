import json

import numpy as np

import water_entropy.cli as cli
from water_entropy.data import WaterObservations


class _Report:
    def to_dict(self):
        return {"ok": True}

    def __str__(self):
        return "report ok"


def _observations(n=20):
    rng = np.random.default_rng(4)
    oxygen = rng.normal(scale=0.1, size=(n, 3))
    local_hydrogen = np.array([[0.586, 0.757, 0.0], [0.586, -0.757, 0.0]])
    return WaterObservations(
        oxygen=oxygen,
        hydrogen=oxygen[:, None, :] + local_hydrogen,
        frame=np.arange(n),
        resid=np.ones(n, dtype=int),
        water_id=np.ones(n, dtype=int),
        frames=np.arange(n),
        times=np.arange(n) * 2.0,
        fit_rmsd=np.zeros(n),
        ligand_reference=np.zeros((1, 3)),
        hb_solute=np.zeros(n),
        enclosure=np.ones(n),
    )


def test_info_command_emits_json(monkeypatch, capsys):
    monkeypatch.setattr(cli, "load_universe", lambda config: object())
    monkeypatch.setattr(cli, "describe_system", lambda universe, config: (_Report(), None))

    assert cli.main(["info", "-s", "top.psf", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == {"ok": True}


def test_check_command_emits_report(monkeypatch, capsys):
    monkeypatch.setattr(cli, "load_universe", lambda config: object())
    monkeypatch.setattr(cli, "run_preprocess_qc", lambda universe, config: _Report())

    assert cli.main(["check", "-s", "top.psf"]) == 0
    assert capsys.readouterr().out.strip() == "report ok"


def test_preprocess_command_writes_cache(monkeypatch, tmp_path):
    observations = _observations()
    destination = tmp_path / "observations.npz"
    monkeypatch.setattr(cli, "load_universe", lambda config: object())
    monkeypatch.setattr(cli, "run_preprocess", lambda universe, config, progress=None: observations)

    assert cli.main(["preprocess", "-s", "top.psf", "-o", str(destination)]) == 0
    assert WaterObservations.load(destination).n_observations == observations.n_observations


def test_sites_and_rank_commands_reuse_a_cache(tmp_path, capsys):
    cache = _observations().save(tmp_path / "observations.npz")
    csv_path = tmp_path / "sites.csv"
    pdb_path = tmp_path / "sites.pdb"
    shared = [
        "-s",
        "unused.psf",
        "--observations",
        str(cache),
        "--density-factor",
        "0.1",
        "--json",
    ]

    assert (
        cli.main(["sites", *shared, "--csv", str(csv_path), "--site-coordinates", str(pdb_path)])
        == 0
    )
    sites = json.loads(capsys.readouterr().out)
    assert len(sites["sites"]) == 1
    assert csv_path.is_file() and pdb_path.is_file()

    assert cli.main(["rank", *shared, "--top", "1"]) == 0
    ranking = json.loads(capsys.readouterr().out)
    assert len(ranking["sites"]) == 1
    assert "category" in ranking["sites"][0]


def test_command_line_options_override_config_file(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "topology: top.psf\ntrajectory: [traj.xtc]\nstep: 2\nwater_cutoff: 4.0\n"
    )
    args = cli.build_parser().parse_args(
        ["sites", "-c", str(config_path), "--step", "3", "--water-cutoff", "6.0"]
    )

    config = cli._config_from_args(args)
    assert config.step == 3
    assert config.water_cutoff == 6.0
    assert config.trajectory == [tmp_path / "traj.xtc"]


def test_top_must_be_positive():
    parser = cli.build_parser()
    try:
        parser.parse_args(["rank", "-s", "top.psf", "--top", "0"])
    except SystemExit as error:
        assert error.code == 2
    else:
        raise AssertionError("--top 0 should be rejected")


def test_progress_reporter_is_throttled_and_can_be_disabled(capsys):
    assert cli._progress_reporter(True) is None
    report = cli._progress_reporter(False)
    for current in range(1, 21):
        report(current, 20)
    output = capsys.readouterr().err
    assert "1/20" in output
    assert "20/20 frames (100%)" in output
