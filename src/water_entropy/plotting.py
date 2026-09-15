"""Optional diagnostic plots for hydration-site analyses."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from water_entropy.analysis import SiteAnalysis
from water_entropy.entropy import minus_t_delta_s, orientational_entropy, translational_entropy
from water_entropy.exceptions import WaterEntropyError
from water_entropy.ranking import SiteRanking


def write_analysis_plots(analysis: SiteAnalysis, output_dir: str | Path) -> list[Path]:
    """Write occupancy, residence-distribution, and entropy-convergence plots."""
    plt = _pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [
        _plot_occupancy(analysis, output_dir / "site_occupancy.png", plt),
        _plot_residence_distributions(analysis, output_dir / "residence_distributions.png", plt),
        _plot_convergence(analysis, output_dir / "entropy_convergence.png", plt),
    ]
    return paths


def write_ranking_plot(ranking: SiteRanking, output_dir: str | Path) -> Path:
    """Write a 3D map of site centres coloured by displacement score."""
    plt = _pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "ranked_site_map.png"
    rows = ranking.rows()
    figure = plt.figure(figsize=(7, 6))
    axis = figure.add_subplot(111, projection="3d")
    if rows:
        xyz = np.array([[row[axis_name] for axis_name in ("x", "y", "z")] for row in rows])
        scores = np.array([row["score"] for row in rows])
        points = axis.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=scores, cmap="viridis", s=55)
        for row, point in zip(rows, xyz, strict=True):
            axis.text(*point, str(row["site"]), fontsize=8)
        figure.colorbar(points, ax=axis, label="displacement score (kcal/mol)", shrink=0.7)
    axis.set(xlabel="x (Å)", ylabel="y (Å)", zlabel="z (Å)", title="Ranked hydration sites")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _plot_occupancy(analysis, path, plt):
    site = np.arange(1, analysis.n_sites + 1)
    figure, axis = plt.subplots(figsize=(max(7, analysis.n_sites * 0.35), 4))
    axis.bar(site, analysis.occupancy)
    axis.set(
        xlabel="site", ylabel="fraction of frames", title="Hydration-site occupancy", ylim=(0, 1)
    )
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _plot_residence_distributions(analysis, path, plt):
    distributions = []
    labels = []
    for index in range(analysis.n_sites):
        members = analysis.sites.labels == index
        lengths = _episode_lengths(
            analysis.observations.frame[members], analysis.observations.water_id[members]
        )
        if lengths:
            distributions.append(np.asarray(lengths) * analysis.observations.dt_ps)
            labels.append(str(index + 1))
    figure, axis = plt.subplots(figsize=(max(7, len(labels) * 0.35), 4))
    if distributions:
        axis.boxplot(distributions, tick_labels=labels, showfliers=False)
    axis.set(xlabel="site", ylabel="episode duration (ps)", title="Residence-time distributions")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _plot_convergence(analysis, path, plt):
    observations = analysis.observations
    checkpoints = np.unique(
        np.linspace(1, observations.n_frames, min(8, observations.n_frames), dtype=int)
    )
    values = np.full((checkpoints.size, analysis.n_sites), np.nan)
    orientations = None
    if observations.n_observations:
        from water_entropy.entropy import water_orientations

        orientations = water_orientations(observations.oxygen, observations.hydrogen)
    for row, stop in enumerate(checkpoints):
        before = observations.frame < stop
        for site in range(analysis.n_sites):
            members = before & (analysis.sites.labels == site)
            s_trans = translational_entropy(observations.oxygen[members])
            s_orient = (
                orientational_entropy(orientations[members]) if orientations is not None else np.nan
            )
            values[row, site] = minus_t_delta_s(s_trans + s_orient, analysis.temperature)
    time_ns = checkpoints * observations.dt_ps / 1000.0
    figure, axis = plt.subplots(figsize=(7, 4))
    for site in range(min(analysis.n_sites, 10)):
        axis.plot(time_ns, values[:, site], marker=".", label=f"site {site + 1}")
    axis.set(xlabel="analysed time (ns)", ylabel="-TΔS (kcal/mol)", title="Entropy convergence")
    if analysis.n_sites:
        axis.legend(fontsize=7, ncol=2)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _episode_lengths(frames: np.ndarray, water_id: np.ndarray) -> list[int]:
    lengths = []
    for identity in np.unique(water_id):
        visits = np.unique(frames[water_id == identity])
        breaks = np.flatnonzero(np.diff(visits) > 1) + 1
        lengths.extend(int(part[-1] - part[0] + 1) for part in np.split(visits, breaks))
    return lengths


def _pyplot():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise WaterEntropyError(
            "plotting requires the optional dependencies; install with 'uv sync --extra plots'"
        ) from None
    return plt
