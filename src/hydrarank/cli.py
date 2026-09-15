"""HydraRank: hydration-site analysis and ligand-displacement ranking."""

from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

from hydrarank import __version__
from hydrarank.analysis import analyse_sites, format_analysis
from hydrarank.clustering import cluster_hydration_sites
from hydrarank.config import PreprocessConfig
from hydrarank.data import WaterObservations
from hydrarank.exceptions import HydraRankError
from hydrarank.export import write_csv, write_site_coordinates
from hydrarank.io import describe_system, load_universe
from hydrarank.jsonio import dumps as json_dumps
from hydrarank.plotting import write_analysis_plots, write_ranking_plot
from hydrarank.preprocess import run_preprocess
from hydrarank.qc import run_preprocess_qc
from hydrarank.ranking import (
    DEFAULT_ENTROPY_THRESHOLD,
    DEFAULT_HBOND_THRESHOLD,
    format_ranking,
    rank_sites,
)
from hydrarank.workflow import run_analysis


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="hydrarank", description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyse = subparsers.add_parser(
        "analyse",
        aliases=["analyze"],
        help="run the complete workflow and write a reproducible results directory",
    )
    _add_system_arguments(analyse)
    _add_analysis_parameters(analyse)
    _add_ranking_parameters(analyse)
    analyse.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="results directory (default: output_dir from config, or ./output)",
    )
    analyse.add_argument(
        "--force", action="store_true", help="ignore a compatible observations cache"
    )
    analyse.add_argument(
        "--allow-qc-failures",
        action="store_true",
        help="continue despite a split-solute QC failure (unsafe unless reviewed)",
    )
    analyse.add_argument(
        "--plots", action="store_true", help="write plots (requires the plots extra)"
    )
    analyse.set_defaults(func=_cmd_analyse)

    info = subparsers.add_parser(
        "info", help="load a system, run the selections and print what was found"
    )
    _add_system_arguments(info)
    info.set_defaults(func=_cmd_info)

    check = subparsers.add_parser(
        "check",
        help="apply the PBC and alignment stack and report the preprocessing diagnostics",
    )
    _add_system_arguments(check)
    check.add_argument("--water-cutoff", type=float, default=None)
    check.add_argument("--pocket-cutoff", type=float, default=None)
    check.set_defaults(func=_cmd_check)

    preprocess = subparsers.add_parser(
        "preprocess",
        help="extract the aligned first-shell waters and write them to disk",
    )
    _add_system_arguments(preprocess)
    preprocess.add_argument("--water-cutoff", type=float, default=None)
    preprocess.add_argument("--pocket-cutoff", type=float, default=None)
    preprocess.add_argument("-o", "--output", type=Path, default=None, help="destination .npz file")
    preprocess.set_defaults(func=_cmd_preprocess)

    sites = subparsers.add_parser(
        "sites", help="cluster the first-shell waters into hydration sites"
    )
    _add_system_arguments(sites)
    _add_site_arguments(sites)
    sites.set_defaults(func=_cmd_sites)

    rank = subparsers.add_parser("rank", help="rank hydration sites as ligand-displacement targets")
    _add_system_arguments(rank)
    _add_site_arguments(rank)
    _add_ranking_parameters(rank)
    rank.set_defaults(func=_cmd_rank)
    return parser


def _add_site_arguments(parser: argparse.ArgumentParser) -> None:
    _add_analysis_parameters(parser)
    parser.add_argument(
        "--observations",
        type=Path,
        default=None,
        help="reuse a .npz written by 'preprocess' instead of reading the trajectory",
    )
    parser.add_argument("--csv", type=Path, default=None, help="write the site table as CSV")
    parser.add_argument(
        "--site-coordinates",
        type=Path,
        default=None,
        help="write site centres as .pdb, .cif or .mmcif",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=None,
        help="write occupancy, residence and convergence plots to this directory",
    )


def _add_analysis_parameters(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--water-cutoff", type=float, default=None)
    parser.add_argument("--pocket-cutoff", type=float, default=None)
    parser.add_argument("--site-radius", type=float, default=None)
    parser.add_argument("--density-factor", type=float, default=None)
    parser.add_argument("--temperature", type=float, default=None)
    parser.add_argument("--max-gap", type=int, default=None)
    parser.add_argument("--max-sites", type=int, default=None)


def _add_ranking_parameters(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--hbond-penalty",
        type=float,
        default=None,
        help="kcal/mol subtracted per mean solute hydrogen bond (default: config or 1.0)",
    )
    parser.add_argument(
        "--entropy-threshold",
        type=float,
        default=DEFAULT_ENTROPY_THRESHOLD,
        help="minimum -TdS for an ordered site (default: %(default)s kcal/mol)",
    )
    parser.add_argument(
        "--hbond-threshold",
        type=float,
        default=DEFAULT_HBOND_THRESHOLD,
        help="mean solute hydrogen bonds requiring replacement (default: %(default)s)",
    )
    parser.add_argument("--top", type=_positive_int, default=None, help="show only the top N sites")


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _add_system_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-c", "--config", type=Path, help="YAML configuration file")
    group.add_argument("-s", "--topology", type=Path, help="topology file (PSF, TPR, PDB, ...)")
    parser.add_argument("-f", "--trajectory", type=Path, nargs="*", default=None)
    parser.add_argument("-l", "--ligand-selection", default=None)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--stop", type=int, default=None)
    parser.add_argument("--step", type=int, default=None)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    parser.add_argument(
        "--no-progress", action="store_true", help="disable preprocessing progress on stderr"
    )


def _config_from_args(args: argparse.Namespace) -> PreprocessConfig:
    if args.config is not None:
        config = PreprocessConfig.from_yaml(args.config)
    else:
        config = PreprocessConfig(topology=args.topology)

    overrides = {}
    for name in (
        "trajectory",
        "ligand_selection",
        "start",
        "stop",
        "step",
        "water_cutoff",
        "pocket_cutoff",
        "site_radius",
        "density_factor",
        "temperature",
        "max_gap",
        "hbond_penalty",
    ):
        if getattr(args, name, None) is not None:
            overrides[name] = getattr(args, name)
    return dataclasses.replace(config, **overrides) if overrides else config


def _emit(report, as_json: bool) -> int:
    print(json_dumps(report.to_dict(), indent=2) if as_json else report)
    return 0


def _cmd_analyse(args: argparse.Namespace) -> int:
    config = _config_from_args(args)
    if args.output_dir is not None:
        config = dataclasses.replace(config, output_dir=args.output_dir)
    result = run_analysis(
        config,
        max_sites=args.max_sites,
        hbond_penalty=args.hbond_penalty,
        entropy_threshold=args.entropy_threshold,
        hbond_threshold=args.hbond_threshold,
        force=args.force,
        allow_qc_failures=args.allow_qc_failures,
        plots=args.plots,
        progress=_progress_reporter(args.no_progress),
    )
    print(json_dumps(result.to_dict(), indent=2) if args.json else result.format(top=args.top))
    return 0


def _cmd_info(args: argparse.Namespace) -> int:
    config = _config_from_args(args)
    report, _ = describe_system(load_universe(config), config)
    return _emit(report, args.json)


def _cmd_check(args: argparse.Namespace) -> int:
    config = _config_from_args(args)
    return _emit(run_preprocess_qc(load_universe(config), config), args.json)


def _cmd_preprocess(args: argparse.Namespace) -> int:
    config = _config_from_args(args)
    observations = run_preprocess(
        load_universe(config), config, progress=_progress_reporter(args.no_progress)
    )
    destination = args.output or config.output_dir / "observations.npz"
    observations.save(destination)
    print(
        f"wrote {observations.n_observations} water observations from "
        f"{observations.n_frames} frames to {destination}"
    )
    return 0


def _cmd_sites(args: argparse.Namespace) -> int:
    config = _config_from_args(args)
    analysis = _analysis_from_args(args, config)
    if args.json:
        print(json_dumps(analysis.to_dict(), indent=2))
    else:
        print(format_analysis(analysis))
    _write_site_artifacts(analysis.rows(), analysis, args)
    return 0


def _analysis_from_args(args: argparse.Namespace, config: PreprocessConfig):
    if args.observations is not None:
        observations = WaterObservations.load(args.observations)
    else:
        observations = run_preprocess(
            load_universe(config), config, progress=_progress_reporter(args.no_progress)
        )

    sites = cluster_hydration_sites(
        observations,
        radius=config.site_radius,
        density_factor=config.density_factor,
        max_sites=args.max_sites,
    )
    return analyse_sites(
        observations, sites, temperature=config.temperature, max_gap=config.max_gap
    )


def _cmd_rank(args: argparse.Namespace) -> int:
    config = _config_from_args(args)
    analysis = _analysis_from_args(args, config)
    penalty = config.hbond_penalty if args.hbond_penalty is None else args.hbond_penalty
    ranking = rank_sites(
        analysis,
        hbond_penalty=penalty,
        entropy_threshold=args.entropy_threshold,
        hbond_threshold=args.hbond_threshold,
    )
    if args.json:
        print(json_dumps(ranking.to_dict(top=args.top), indent=2))
    else:
        print(format_ranking(ranking, top=args.top))
    rows = ranking.rows()
    if args.top is not None:
        rows = rows[: args.top]
    _write_site_artifacts(rows, analysis, args)
    if args.plot_dir is not None:
        path = write_ranking_plot(ranking, args.plot_dir)
        print(f"wrote {path}", file=sys.stderr)
    return 0


def _write_site_artifacts(rows, analysis, args: argparse.Namespace) -> None:
    if args.csv is not None:
        print(f"wrote {write_csv(rows, args.csv)}", file=sys.stderr)
    if args.site_coordinates is not None:
        path = write_site_coordinates(rows, args.site_coordinates)
        print(f"wrote {path}", file=sys.stderr)
    if args.plot_dir is not None:
        for path in write_analysis_plots(analysis, args.plot_dir):
            print(f"wrote {path}", file=sys.stderr)


def _progress_reporter(disabled: bool):
    if disabled:
        return None
    last_percent = -1

    def report(current: int, total: int) -> None:
        nonlocal last_percent
        percent = int(current * 100 / total)
        if percent >= last_percent + 5 or current == total:
            ending = "\n" if current == total else "\r"
            print(
                f"Preprocessing: {current}/{total} frames ({percent}%)",
                end=ending,
                file=sys.stderr,
                flush=True,
            )
            last_percent = percent

    return report


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (HydraRankError, OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
