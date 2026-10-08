"""Sequential benchmark pilot. A full run requires an explicit --full-study flag."""

import argparse
import json
import os
import platform
import random
import shutil
import subprocess
import time
from pathlib import Path

from common import MODES, ROOT, digest, dump
from monitor import physical_ram, run

ADAPTER = Path(__file__).with_name("adapter.py")
HYDRA = ROOT / ".venv/bin/python"
SST = ROOT / "sstmap_analysis/3rlp/.env/bin/python"


def environment(workers, folder):
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["MPLCONFIGDIR"] = str(folder / "matplotlib-cache")
    env["XDG_CACHE_HOME"] = str(folder / "cache")
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONFAULTHANDLER"] = "1"
    if workers == "single":
        for key in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            env[key] = "1"
    else:
        for key in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            env.pop(key, None)
    return env


def adapter_command(mode, folder, inputs, frames, workers, spacing=0.5, padding=5):
    python = HYDRA if mode in ("hydrarank", "prepare") else SST
    cmd = [
        python,
        "-B",
        ADAPTER,
        "--mode",
        mode,
        "--output",
        folder,
        "--frames",
        str(frames),
        "--workers",
        workers,
        "--spacing",
        str(spacing),
        "--padding",
        str(padding),
    ]
    if inputs:
        cmd += ["--inputs", inputs]
    return cmd


def warm_prefix(inputs, frames):
    # Inputs are finite exported prefixes, not multi-gigabyte production files.
    for name in ("raw.xtc", "aligned.xtc", "ligand.pdb"):
        with (inputs / name).open("rb") as f:
            while f.read(1024 * 1024):
                pass


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--budget-seconds", type=float, default=7200)
    p.add_argument("--trial-timeout-seconds", type=float, default=1200)
    p.add_argument("--full-study", action="store_true")
    args = p.parse_args()
    if args.budget_seconds <= 0 or args.trial_timeout_seconds <= 0:
        p.error("budgets and timeouts must be positive")
    if not args.full_study and args.trial_timeout_seconds > 1200:
        p.error("pilot trials are limited to twenty minutes")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    ram = physical_ram()
    hardware = dict(
        platform=platform.platform(),
        cpu_count=os.cpu_count(),
        physical_ram_bytes=ram,
        model=subprocess.check_output(["/usr/sbin/sysctl", "-n", "hw.model"], text=True).strip(),
        git_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        git_diff=subprocess.check_output(
            ["git", "diff", "HEAD", "--", "src/hydrarank"], cwd=ROOT, text=True
        ),
        source_hashes={
            str(f.relative_to(ROOT)): digest(f) for f in Path(__file__).parent.glob("*.py")
        },
        filesystem_cache="untimed prefix read; no cache purge",
        sample_interval_seconds=0.5,
    )
    dump(root / "hardware.json", hardware)
    snapshot = root / "harness_snapshot"
    snapshot.mkdir()
    for source in hardware["source_hashes"]:
        shutil.copy2(ROOT / source, snapshot / Path(source).name)

    def trial(mode, name, frames, workers="single", inputs=None, spacing=0.5, padding=5):
        folder = root / name
        warm = inputs if inputs else None
        if warm:
            warm_prefix(warm, frames)
        remaining = args.budget_seconds - (time.monotonic() - started)
        if remaining <= 0:
            return None
        metrics = run(
            adapter_command(mode, folder, inputs, frames, workers, spacing, padding),
            folder,
            env=environment(workers, folder),
            timeout=min(args.trial_timeout_seconds, remaining),
            memory_limit=ram // 2,
        )
        metrics.update(
            mode=mode, frames=frames, workers=workers, name=name, spacing=spacing, padding=padding
        )
        if (folder / "adapter.json").exists():
            metrics["adapter"] = json.loads((folder / "adapter.json").read_text())
            if metrics["adapter"]["status"] != "complete":
                metrics["status"] = "failed"
        dump(folder / "resources.json", metrics)
        return metrics

    records = []
    # Environment inventories are outside timed scientific comparisons.
    for name, python in [("hydrarank", HYDRA), ("sstmap", SST)]:
        folder = root / ("environment_" + name)
        folder.mkdir()
        subprocess.run(
            [str(python), "-B", str(ADAPTER), "--mode", "manifest", "--output", str(folder)],
            env=environment("single", folder),
            check=True,
            stdout=subprocess.DEVNULL,
        )
    # Monitor overhead: paired fresh processes doing the same CPU workload.
    for repeat in range(3):
        for monitored in (False, True):
            folder = root / f"overhead_{repeat}_{monitored}"
            workload = [
                HYDRA,
                "-B",
                "-c",
                'import hashlib; b=b"x"*1048576; [hashlib.sha256(b).digest() for _ in range(1500)]',
            ]
            records.append(
                dict(
                    kind="monitor_overhead",
                    repeat=repeat,
                    **run(workload, folder, monitored=monitored, timeout=30),
                )
            )
    max_frames = 25001 if args.full_study else 1001
    prepared = trial("prepare", "preparation", max_frames)
    records.append(dict(kind="preparation", **prepared))
    inputs = root / "preparation"
    if prepared["status"] != "complete":
        dump(root / "trials.json", records)
        raise RuntimeError("Preparation validation failed; see logs")
    known_inputs = [
        ROOT / "sample_traj/3rlp/step5.tpr",
        ROOT / "sample_traj/3rlp/step5.xtc",
        ROOT / "sstmap_analysis/3rlp/input/3rlp_sstmap_topology.gro",
        ROOT / "sstmap_analysis/3rlp/input/3rlp_hydrarank61_aligned.xtc",
        ROOT / "sample_traj/3rlp/charmm-gui-8947372302/gromacs/topol.top",
    ]
    parameter_root = ROOT / "sample_traj/3rlp/charmm-gui-8947372302/gromacs"
    known_inputs += list((parameter_root / "toppar").glob("*")) + list(
        (parameter_root / "../toppar").glob("*")
    )
    known_inputs += list(inputs.glob("*.xtc")) + list(inputs.glob("*.pdb"))
    dump(
        root / "inputs.json",
        [
            dict(path=str(f), bytes=f.stat().st_size, sha256=digest(f))
            for f in known_inputs
            if f.is_file()
        ],
    )
    # Preflight runs are excluded from speed summaries.
    valid = []
    for mode in MODES:
        result = trial(mode, "preflight_" + mode, 101, inputs=inputs)
        if result:
            result["kind"] = "preflight"
            records.append(result)
            if result["status"] == "complete":
                valid.append(mode)
        dump(root / "trials.json", records)
    if args.full_study:
        schedule = [
            (n, rep, mode, "single", 0.5, 5)
            for n in (101, 501, 1001, 2501, 5001, 10001, 25001)
            for rep in range(3)
            for mode in valid
        ]
        schedule += [(n, 0, mode, "native", 0.5, 5) for n in (1001, 25001) for mode in valid]
        # 0.25 and 1 A are blocked in this installed GIST implementation, not silently benchmarked.
        schedule += [
            (1001, rep, mode, "single", 0.5, pad)
            for pad in (7.5, 10)
            for rep in range(3)
            for mode in valid
            if mode.startswith("gist")
        ]
    else:
        schedule = [
            (n, rep, mode, "single", 0.5, 5)
            for n, reps in ((101, 3), (501, 3), (1001, 1))
            for rep in range(reps)
            for mode in valid
        ]
        schedule += [(501, 0, mode, "native", 0.5, 5) for mode in valid]
    # Balanced deterministic random order within each frame size.
    random.Random(20261008).shuffle(schedule)
    schedule.sort(key=lambda row: row[0])
    dump(
        root / "schedule.json",
        [
            dict(frames=n, repeat=r, mode=m, workers=w, spacing=s, padding=pad)
            for n, r, m, w, s, pad in schedule
        ],
    )
    for index, (n, rep, mode, workers, spacing, padding) in enumerate(schedule):
        print(
            f"Trial {index + 1}/{len(schedule)}: {mode}, {n} frames, {workers}, repeat {rep + 1}",
            flush=True,
        )
        result = trial(
            mode, f"trial_{index:03d}_{mode}_{n}_{workers}", n, workers, inputs, spacing, padding
        )
        if result is None:
            break
        result.update(kind="benchmark", repeat=rep)
        records.append(result)
        dump(root / "trials.json", records)
    # Cached trials are always distinct and explicitly labelled.
    cache_counts = (1001, 25001) if args.full_study else (501,)
    for count in cache_counts:
        sources = [
            r
            for r in records
            if r.get("kind") == "benchmark"
            and r.get("mode") == "hydrarank"
            and r["frames"] == count
            and r["workers"] == "single"
            and r["status"] == "complete"
        ]
        if sources and time.monotonic() - started < args.budget_seconds:
            folder = root / f"cached_hydrarank_{count}"
            # Adapter output must be a new directory; preload only the observations cache.
            source = root / sources[0]["name"] / "results/observations.npz"
            # run() owns directory creation; the worker copies the cache before adapter start.
            command = [
                HYDRA,
                "-B",
                Path(__file__).with_name("cache_worker.py"),
                str(source),
                *adapter_command("hydrarank", folder, inputs, count, "single"),
            ]
            result = run(
                command,
                folder,
                env=environment("single", folder),
                timeout=min(
                    args.trial_timeout_seconds, args.budget_seconds - (time.monotonic() - started)
                ),
                memory_limit=ram // 2,
            )
            result.update(
                kind="cached", mode="hydrarank", frames=count, workers="single", name=folder.name
            )
            if (folder / "adapter.json").exists():
                result["adapter"] = json.loads((folder / "adapter.json").read_text())
            records.append(result)
            dump(root / "trials.json", records)
    dump(
        root / "completion.json",
        dict(
            elapsed_seconds=time.monotonic() - started,
            budget_seconds=args.budget_seconds,
            full_study=args.full_study,
            valid_modes=valid,
            scheduled=len(schedule),
            completed=sum(r.get("kind") == "benchmark" for r in records),
        ),
    )
    print(f"Campaign complete: {root}", flush=True)


if __name__ == "__main__":
    main()
