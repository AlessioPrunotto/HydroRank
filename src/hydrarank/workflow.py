"""The unified, reproducible HydraRank analysis workflow."""

from __future__ import annotations

import hashlib
import platform
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from hydrarank import __version__
from hydrarank.analysis import SiteAnalysis, analyse_sites
from hydrarank.clustering import cluster_hydration_sites
from hydrarank.config import PreprocessConfig
from hydrarank.data import WaterObservations
from hydrarank.exceptions import HydraRankError
from hydrarank.export import write_csv, write_site_coordinates
from hydrarank.io import SystemReport, describe_system, load_universe
from hydrarank.jsonio import dumps as json_dumps
from hydrarank.plotting import write_analysis_plots, write_ranking_plot
from hydrarank.preprocess import run_preprocess
from hydrarank.qc import PreprocessQC, qc_from_observations
from hydrarank.ranking import SiteRanking, format_ranking, rank_sites


@dataclass
class AnalyseResult:
    """In-memory results and artifacts produced by :func:`run_analysis`."""

    config: PreprocessConfig
    system: SystemReport
    qc: PreprocessQC
    analysis: SiteAnalysis
    ranking: SiteRanking
    output_dir: Path
    cache_reused: bool
    artifacts: dict[str, Path] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    created_at: str = field(default_factory=lambda: datetime.now(UTC).isoformat())

    def to_dict(self) -> dict[str, Any]:
        return {
            "hydrarank_version": __version__,
            "created_at": self.created_at,
            "python_version": platform.python_version(),
            "configuration": self.config.to_dict(),
            "system": self.system.to_dict(),
            "preprocessing_qc": self.qc.to_dict(),
            "cache_reused": self.cache_reused,
            "warnings": self.warnings,
            "ranking": self.ranking.to_dict(),
            "artifacts": {name: str(path) for name, path in self.artifacts.items()},
        }

    def format(self, top: int | None = None) -> str:
        lines = [
            f"HydraRank {__version__} analysis",
            f"Results directory: {self.output_dir}",
            f"Observations cache: {'reused' if self.cache_reused else 'created'}",
            "",
            str(self.system),
            "",
            str(self.qc),
            "",
            format_ranking(self.ranking, top=top),
        ]
        workflow_warnings = [item for item in self.warnings if item not in self.system.warnings]
        if workflow_warnings:
            lines.extend(["", "Workflow warnings", *[f"  - {item}" for item in workflow_warnings]])
        if self.artifacts:
            artifact_lines = [f"  {name:<18}: {path}" for name, path in self.artifacts.items()]
            lines.extend(["", "Saved results", *artifact_lines])
        return "\n".join(lines)


def run_analysis(
    config: PreprocessConfig,
    *,
    output_dir: str | Path | None = None,
    max_sites: int | None = None,
    hbond_penalty: float | None = None,
    entropy_threshold: float = 2.0,
    hbond_threshold: float = 1.0,
    force: bool = False,
    allow_qc_failures: bool = False,
    plots: bool = False,
    progress=None,
) -> AnalyseResult:
    """Run inspection, preprocessing, QC, site analysis, ranking, and export."""
    if output_dir is not None:
        config = replace(config, output_dir=Path(output_dir))
    destination = config.output_dir
    destination.mkdir(parents=True, exist_ok=True)
    cache_path = destination / "observations.npz"
    cache_key = _cache_key(config)

    universe = load_universe(config)
    system, _ = describe_system(universe, config)

    cache_reused = False
    if cache_path.is_file() and not force:
        candidate = WaterObservations.load(cache_path)
        if candidate.metadata.get("hydrarank_cache_key") == cache_key:
            observations = candidate
            cache_reused = True
        else:
            observations = _preprocess(universe, config, cache_key, progress)
    else:
        observations = _preprocess(universe, config, cache_key, progress)
    if not cache_reused:
        observations.save(cache_path)

    qc = qc_from_observations(observations)
    if not qc.solute_always_whole and not allow_qc_failures:
        raise HydraRankError(
            "preprocessing QC found a split solute, so ranking was stopped; check topology "
            "bonds/PBC handling, or pass --allow-qc-failures after reviewing the risk"
        )

    sites = cluster_hydration_sites(
        observations,
        radius=config.site_radius,
        density_factor=config.density_factor,
        max_sites=max_sites,
    )
    analysis = analyse_sites(
        observations,
        sites,
        temperature=config.temperature,
        max_gap=config.max_gap,
    )
    ranking = rank_sites(
        analysis,
        hbond_penalty=config.hbond_penalty if hbond_penalty is None else hbond_penalty,
        entropy_threshold=entropy_threshold,
        hbond_threshold=hbond_threshold,
    )

    warnings = list(system.warnings)
    if qc.fit_rmsd.max() > 2.5:
        warnings.append("binding-site fit RMSD exceeds 2.5 A; site disorder may be inflated")
    if not ranking.rows():
        warnings.append("no hydration sites passed the configured density threshold")

    result = AnalyseResult(
        config=config,
        system=system,
        qc=qc,
        analysis=analysis,
        ranking=ranking,
        output_dir=destination,
        cache_reused=cache_reused,
        warnings=warnings,
    )
    _write_artifacts(result, cache_path, plots=plots)
    return result


def _preprocess(universe, config, cache_key, progress):
    observations = run_preprocess(universe, config, progress=progress, collect_qc=True)
    observations.metadata["hydrarank_cache_key"] = cache_key
    return observations


def _write_artifacts(result: AnalyseResult, cache_path: Path, *, plots: bool) -> None:
    destination = result.output_dir
    result.artifacts["observations"] = cache_path

    config_path = destination / "config.yaml"
    result.config.to_yaml(config_path)
    result.artifacts["configuration"] = config_path

    rows = result.ranking.rows()
    if rows:
        result.artifacts["ranking_csv"] = write_csv(rows, destination / "ranking.csv")
        result.artifacts["site_coordinates"] = write_site_coordinates(
            rows, destination / "sites.pdb"
        )
    if plots:
        plot_paths = write_analysis_plots(result.analysis, destination / "plots")
        plot_paths.append(write_ranking_plot(result.ranking, destination / "plots"))
        for path in plot_paths:
            result.artifacts[f"plot_{path.stem}"] = path

    json_path = destination / "results.json"
    result.artifacts["results_json"] = json_path
    report_path = destination / "report.txt"
    result.artifacts["report"] = report_path
    json_path.write_text(json_dumps(result.to_dict(), indent=2) + "\n")
    report_path.write_text(result.format() + "\n")


def _cache_key(config: PreprocessConfig) -> str:
    preprocessing_parameters = (
        "ligand_selection",
        "water_resnames",
        "ion_resnames",
        "start",
        "stop",
        "step",
        "water_cutoff",
        "pocket_cutoff",
        "min_ligand_heavy_atoms",
        "hbond_distance",
        "hbond_angle",
        "enclosure_radius",
    )
    payload = {name: getattr(config, name) for name in preprocessing_parameters}
    inputs = [config.topology, *config.trajectory]
    payload["inputs"] = [_file_signature(Path(path)) for path in inputs]
    encoded = json_dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _file_signature(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    try:
        stat = resolved.stat()
    except FileNotFoundError:
        # load_universe will produce the user-facing error; this keeps key generation simple.
        return {"path": str(resolved), "missing": True}
    return {"path": str(resolved), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
