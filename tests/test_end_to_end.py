from hydrarank.analysis import analyse_sites
from hydrarank.clustering import cluster_hydration_sites
from hydrarank.config import PreprocessConfig
from hydrarank.preprocess import run_preprocess
from hydrarank.ranking import rank_sites


def test_compact_in_memory_pipeline_runs_end_to_end(pocket_system, tmp_path):
    config = PreprocessConfig(
        topology=tmp_path / "top.psf",
        ligand_selection="resname 547",
        pocket_cutoff=10.0,
    )

    observations = run_preprocess(pocket_system, config)
    sites = cluster_hydration_sites(observations, min_count=1, max_sites=2)
    analysis = analyse_sites(observations, sites)
    ranking = rank_sites(analysis)

    assert observations.n_observations == 8
    assert sites.n_sites > 0
    assert len(analysis.rows()) == sites.n_sites
    assert len(ranking.rows()) == sites.n_sites
