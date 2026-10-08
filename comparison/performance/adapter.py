"""Run inside either existing environment; no production code is modified."""

import argparse
import importlib.metadata
import os
import platform
import sys
from pathlib import Path

from common import ROOT, Timings, digest, dump


def manifest():
    packages = {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()}
    files = list((ROOT / "src/hydrarank").glob("*.py"))
    if "sstmap" in packages or "SSTMap" in packages:
        import sstmap

        base = Path(sstmap.__file__).parent
        files += list(base.glob("*.py")) + list(base.glob("*.so"))
        files += list(base.parent.glob("*sstmap*.so"))
    return dict(
        python=sys.version,
        executable=sys.executable,
        machine=platform.machine(),
        packages=packages,
        files={str(p): digest(p) for p in files},
    )


def hydra(args, t):
    from dataclasses import replace

    import numpy as np

    import hydrarank.analysis as analysis
    import hydrarank.clustering as clustering
    import hydrarank.entropy as entropy
    import hydrarank.workflow as workflow
    from hydrarank import preprocess
    from hydrarank.config import PreprocessConfig

    if args.workers == "single":
        from scipy.spatial import cKDTree

        class SingleTree:
            def __init__(self, *a, **kw):
                self.tree = cKDTree(*a, **kw)

            def __getattr__(self, name):
                return getattr(self.tree, name)

            def query(self, *a, **kw):
                kw["workers"] = 1
                return self.tree.query(*a, **kw)

            def query_ball_point(self, *a, **kw):
                kw["workers"] = 1
                return self.tree.query_ball_point(*a, **kw)

        clustering.cKDTree = entropy.cKDTree = SingleTree
    for name in (
        "load_universe",
        "describe_system",
        "run_preprocess",
        "cluster_hydration_sites",
        "analyse_sites",
        "rank_sites",
        "_write_artifacts",
    ):
        t.wrap(workflow, name)
    t.wrap(preprocess, "prepare_system")
    for name in ("translational_entropy", "orientational_entropy"):
        t.wrap(analysis, name)
    config = replace(
        PreprocessConfig.from_yaml(ROOT / "validation/hsp90/hsp90_3rlp.yaml"),
        topology=ROOT / "sample_traj/3rlp/step5.tpr",
        trajectory=[Path(args.inputs) / "raw.xtc"],
        stop=args.frames,
        output_dir=Path(args.output) / "results",
    )
    result = workflow.run_analysis(config)
    assert result.qc.frames.size == args.frames
    assert result.qc.n_fit_atoms == 61
    assert np.all(np.isfinite(result.analysis.minus_t_delta_s))
    return dict(
        actual_frames=int(result.qc.frames.size),
        n_sites=result.analysis.n_sites,
        cache_reused=result.cache_reused,
        n_fit_atoms=result.qc.n_fit_atoms,
        mean_shell_waters=result.qc.n_shell_waters.mean().item(),
        warnings=result.warnings,
    )


def sst(args, t):
    import mdtraj as md
    import numpy as np
    from sstmap import GridWaterAnalysis, SiteWaterAnalysis

    inp = Path(args.inputs)
    params = dict(
        topology_file=str(ROOT / "sstmap_analysis/3rlp/input/3rlp_sstmap_topology.gro"),
        trajectory=str(inp / "aligned.xtc"),
        start_frame=0,
        num_frames=args.frames,
        ligand_file=str(inp / "ligand.pdb"),
        supporting_file=str(ROOT / "sample_traj/3rlp/charmm-gui-8947372302/gromacs/topol.top"),
        rho_bulk=0.0333,
        prefix="benchmark",
    )
    flags = dict(
        energy=args.mode in ("hsa_full", "gist_full"),
        entropy=True,
        hbonds=args.mode != "gist_entropy",
    )
    grid = None
    with t.stage("initialization"):
        if args.mode.startswith("hsa"):
            obj = SiteWaterAnalysis(hsa_region_radius=5.0, **params)
        else:
            assert args.spacing == 0.5, (
                "Installed GIST C voxel assignment hardcodes 0.5 A; "
                "other resolutions need a validated patch."
            )
            lig = md.load_pdb(params["ligand_file"])
            xyz = lig.xyz[0] * 10
            low = np.floor((xyz.min(axis=0) - args.padding) / args.spacing) * args.spacing
            high = np.ceil((xyz.max(axis=0) + args.padding) / args.spacing) * args.spacing
            dims = np.rint((high - low) / args.spacing).astype(int)
            obj = GridWaterAnalysis(
                grid_center=((low + high) / 2).tolist(),
                grid_dimensions=dims.tolist(),
                grid_resolution=[args.spacing] * 3,
                **params,
            )
            assert np.allclose(obj.origin, low)
            assert np.allclose(obj.voxeldata[:, 1:4].min(axis=0), low + args.spacing / 2)
            assert np.allclose(obj.voxeldata[:, 1:4].max(axis=0), high - args.spacing / 2)
            grid = dict(
                origin=low.tolist(),
                upper=high.tolist(),
                dimensions=dims.tolist(),
                spacing=args.spacing,
                volume_A3=float(np.prod(high - low)),
            )
            dump(Path(args.output) / "grid.json", dict(grid=grid, flags=flags))
    assert len(obj.all_atom_ids) == 62822
    assert len(obj.wat_oxygen_atom_ids) == 14799
    assert np.array_equal(obj.neighbor_ids, obj.wat_oxygen_atom_ids)
    assert all(obj.topology.atom(int(i)).residue.name == "OPC" for i in obj.wat_oxygen_atom_ids)
    if args.mode.startswith("hsa"):
        with t.stage("clustering"):
            obj.initialize_hydration_sites(clustering_density_cutoff=2.0)
        for name in (
            "generate_data_for_entropycalcs",
            "run_entropy_scripts",
            "normalize_site_quantities",
        ):
            t.wrap(obj, name)
        with t.stage("site_quantities"):
            obj.calculate_site_quantities(**flags)
        assert obj.num_frames == args.frames
        assert np.isfinite(obj.hsa_data[:, 14:17]).all()
        with t.stage("export"):
            obj.write_calculation_summary()
            obj.write_data()
        summary = dict(
            n_sites=int(obj.hsa_data.shape[0]),
            mean_site_waters=float(obj.hsa_data[:, 4].sum() / args.frames),
        )
    else:
        t.wrap(obj, "calculate_entropy")
        with t.stage("grid_quantities"):
            obj.calculate_grid_quantities(**flags)
        assert obj.num_frames == args.frames
        assert obj.voxeldata[:, 4].sum() > 0
        assert np.isfinite(obj.voxeldata[:, 5:13]).all()
        with t.stage("export"):
            obj.write_data()
            obj.generate_dx_files()
        summary = dict(grid=grid, mean_grid_waters=float(obj.voxeldata[:, 4].sum() / args.frames))
    return dict(actual_frames=int(obj.num_frames), flags=flags, **summary)


def prepare(args, t):
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.analysis.align import rotation_matrix

    from hydrarank.config import PreprocessConfig
    from hydrarank.io import load_universe
    from hydrarank.preprocess import prepare_system

    inp = Path(args.output)
    inp.mkdir(exist_ok=True)
    config = PreprocessConfig.from_yaml(ROOT / "validation/hsp90/hsp90_3rlp.yaml")
    raw = load_universe(config)
    comparison_topology = mda.Universe(
        str(ROOT / "sstmap_analysis/3rlp/input/3rlp_sstmap_topology.gro")
    )
    for attribute in ("names", "resnames", "resindices"):
        assert np.array_equal(
            getattr(raw.atoms, attribute), getattr(comparison_topology.atoms, attribute)
        ), f"SSTMap atom-order mismatch: {attribute}"
    # Reading original frames into a distinct file keeps source offset caches unchanged.
    with (
        t.stage("raw_prefix_export"),
        mda.Writer(str(inp / "raw.xtc"), n_atoms=len(raw.atoms)) as w,
    ):
        for _ts in raw.trajectory[: args.frames]:
            w.write(raw.atoms)
    with t.stage("alignment_initialization"):
        system = prepare_system(raw, config)
    assert system.fit_group.n_atoms == 61
    rmsds = []
    with (
        t.stage("alignment_and_export"),
        mda.Writer(str(inp / "aligned.xtc"), n_atoms=len(raw.atoms)) as w,
    ):
        for ts in raw.trajectory[: args.frames]:
            motion = system.aligner.fit(system.fit_group.positions)
            rmsds.append(motion.rmsd)
            raw.atoms.positions = motion.apply(raw.atoms.positions)
            w.write(raw.atoms)
            if ts.frame == 0:
                system.ligand_heavy.write(str(inp / "ligand.pdb"))
    # Check representative frames against the historical harmonized trajectory.
    reference = mda.Universe(
        str(ROOT / "sample_traj/3rlp/step5.tpr"),
        str(ROOT / "sstmap_analysis/3rlp/input/3rlp_hydrarank61_aligned.xtc"),
    )
    generated = mda.Universe(str(ROOT / "sample_traj/3rlp/step5.tpr"), str(inp / "aligned.xtc"))
    reference.trajectory[0]
    generated.trajectory[0]
    ix = system.fit_group.indices
    source = generated.atoms[ix].positions.astype(float)
    target = reference.atoms[ix].positions.astype(float)
    rotation, _ = rotation_matrix(source - source.mean(0), target - target.mean(0))
    errors = []
    for frame in sorted(set((0, args.frames // 2, args.frames - 1))):
        generated.trajectory[frame]
        reference.trajectory[frame]
        moved = (generated.atoms[ix].positions - source.mean(0)) @ rotation.T + target.mean(0)
        errors.append(
            float(np.sqrt(np.mean(np.sum((moved - reference.atoms[ix].positions) ** 2, axis=1))))
        )
    assert max(errors) < 0.02, errors
    dump(
        inp / "preparation.json",
        dict(
            n_frames=args.frames,
            n_fit_atoms=61,
            historical_fit_group_rmsd_A=errors,
            max_fit_rmsd_A=max(rmsds),
            stages=t.events,
        ),
    )
    return dict(actual_frames=args.frames, historical_fit_group_rmsd_A=errors)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", required=True)
    p.add_argument("--inputs")
    p.add_argument("--output", required=True)
    p.add_argument("--frames", type=int, default=101)
    p.add_argument("--workers", choices=("single", "native"), default="single")
    p.add_argument("--spacing", type=float, default=0.5)
    p.add_argument("--padding", type=float, default=5.0)
    args = p.parse_args()
    if args.inputs:
        args.inputs = str(Path(args.inputs).resolve())
    output = Path(args.output).resolve()
    output.mkdir(exist_ok=True)
    os.chdir(output)
    t = Timings()
    try:
        if args.mode == "manifest":
            dump(output / "environment.json", manifest())
            return
        if args.mode == "prepare":
            result = prepare(args, t)
        elif args.mode == "hydrarank":
            result = hydra(args, t)
        else:
            result = sst(args, t)
        dump(output / "adapter.json", dict(status="complete", result=result, stages=t.events))
    except BaseException as error:
        dump(output / "adapter.json", dict(status="failed", error=repr(error), stages=t.events))
        raise


if __name__ == "__main__":
    main()
