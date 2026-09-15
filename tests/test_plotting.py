import numpy as np

from hydrarank.analysis import analyse_sites
from hydrarank.clustering import HydrationSites
from hydrarank.data import WaterObservations
from hydrarank.plotting import write_analysis_plots, write_ranking_plot
from hydrarank.ranking import rank_sites


def _analysis():
    rng = np.random.default_rng(8)
    n = 20
    oxygen = rng.normal(scale=0.2, size=(n, 3))
    local_hydrogen = np.array([[0.586, 0.757, 0.0], [0.586, -0.757, 0.0]])
    observations = WaterObservations(
        oxygen=oxygen,
        hydrogen=oxygen[:, None, :] + local_hydrogen,
        frame=np.arange(n),
        resid=np.ones(n, dtype=int),
        water_id=np.ones(n, dtype=int),
        frames=np.arange(n),
        times=np.arange(n) * 2.0,
        fit_rmsd=np.zeros(n),
        ligand_reference=np.zeros((1, 3)),
    )
    sites = HydrationSites(
        centers=np.array([oxygen.mean(axis=0)]),
        labels=np.zeros(n, dtype=int),
        radius=1.0,
        n_frames=n,
        min_count=1,
    )
    return analyse_sites(observations, sites)


def test_all_plots_are_written(tmp_path, monkeypatch):
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "matplotlib"))
    analysis = _analysis()

    paths = write_analysis_plots(analysis, tmp_path / "plots")
    paths.append(write_ranking_plot(rank_sites(analysis), tmp_path / "plots"))

    assert {path.name for path in paths} == {
        "site_occupancy.png",
        "residence_distributions.png",
        "entropy_convergence.png",
        "ranked_site_map.png",
    }
    assert all(path.stat().st_size > 0 for path in paths)
