"""Summarize a completed pilot without claiming full-trajectory speedups."""
# ruff: noqa: E501 -- Report paragraphs and Markdown table rows are generated verbatim.

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

from common import dump


def summarize(records):
    groups = defaultdict(list)
    for row in records:
        if row.get("kind") == "benchmark" and row["status"] == "complete":
            groups[
                (row["mode"], row["frames"], row["workers"], row["spacing"], row["padding"])
            ].append(row)
    rows = []
    for key, values in sorted(groups.items()):
        row = dict(zip(("mode", "frames", "workers", "spacing", "padding"), key, strict=True))
        wall = [v.get("os_wall_seconds", v["wall_seconds"]) for v in values]
        row.update(
            repetitions=len(values),
            wall_median_s=statistics.median(wall),
            wall_min_s=min(wall),
            wall_max_s=max(wall),
            frames_per_second=key[1] / statistics.median(wall),
        )
        for field in (
            "user_seconds",
            "system_seconds",
            "os_peak_rss_bytes",
            "sampled_tree_peak_rss_bytes",
            "output_bytes",
            "sampled_peak_disk_bytes",
            "peak_threads",
            "peak_processes",
        ):
            found = [v[field] for v in values if v.get(field) is not None]
            row[field] = statistics.median(found) if found else None
        cpu = [
            v["user_seconds"] + v["system_seconds"]
            for v in values
            if "user_seconds" in v and "system_seconds" in v
        ]
        row["cpu_seconds_per_frame"] = statistics.median(cpu) / key[1] if cpu else None
        rows.append(row)
    return rows


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("campaign", type=Path)
    p.add_argument("--publish", type=Path, default=Path(__file__).parent / "results")
    args = p.parse_args()
    out = args.publish
    out.mkdir(parents=True, exist_ok=True)
    records = json.loads((args.campaign / "trials.json").read_text())
    hardware = json.loads((args.campaign / "hardware.json").read_text())
    rows = summarize(records)
    dump(out / "trials.json", records)
    dump(out / "summary.json", rows)
    dump(out / "hardware.json", hardware)
    for name in ("inputs.json", "completion.json", "schedule.json", "validation.json"):
        if (args.campaign / name).exists():
            (out / name).write_bytes((args.campaign / name).read_bytes())
    for tool in ("hydrarank", "sstmap"):
        (out / f"environment_{tool}.json").write_bytes(
            (args.campaign / f"environment_{tool}/environment.json").read_bytes()
        )
    flat = []
    for r in records:
        flat.append(
            {
                k: r.get(k)
                for k in (
                    "kind",
                    "mode",
                    "frames",
                    "workers",
                    "repeat",
                    "status",
                    "wall_seconds",
                    "os_wall_seconds",
                    "user_seconds",
                    "system_seconds",
                    "os_peak_rss_bytes",
                    "sampled_tree_peak_rss_bytes",
                    "output_bytes",
                    "sampled_peak_disk_bytes",
                    "peak_threads",
                    "peak_processes",
                    "name",
                )
            }
        )
        result = r.get("adapter", {}).get("result", {})
        flat[-1].update(
            {
                key: result.get(key)
                for key in (
                    "actual_frames",
                    "n_sites",
                    "n_fit_atoms",
                    "mean_shell_waters",
                    "mean_site_waters",
                    "mean_grid_waters",
                    "cache_reused",
                )
            }
        )
    write_csv(out / "trials.csv", flat)
    write_csv(out / "summary.csv", rows)
    stage_rows = []
    for r in records:
        events = r.get("adapter", {}).get("stages", [])
        aggregated = defaultdict(list)
        for e in events:
            aggregated[(e["name"], e["depth"])].append(e["seconds"])
        for (name, depth), seconds in aggregated.items():
            stage_rows.append(
                dict(
                    trial=r.get("name"),
                    kind=r.get("kind"),
                    mode=r.get("mode"),
                    frames=r.get("frames"),
                    stage=name,
                    depth=depth,
                    calls=len(seconds),
                    inclusive_seconds=sum(seconds),
                )
            )
    write_csv(out / "stages.csv", stage_rows)
    # Copy all logs and raw monitoring traces into the publishable evidence bundle, not simulation output.
    import shutil

    if (args.campaign / "harness_snapshot").exists():
        shutil.make_archive(
            str(out / "harness_snapshot"), "zip", args.campaign / "harness_snapshot"
        )
    for folder in args.campaign.iterdir():
        if folder.is_dir():
            target = out / "logs" / folder.name
            target.mkdir(parents=True, exist_ok=True)
            for name in (
                "stdout.log",
                "stderr.log",
                "resources.json",
                "samples.json",
                "adapter.json",
                "preparation.json",
                "grid.json",
            ):
                if (folder / name).exists():
                    shutil.copy2(folder / name, target / name)
    overhead = {
        m: [
            r.get("os_wall_seconds", r["wall_seconds"])
            for r in records
            if r.get("kind") == "monitor_overhead" and r["monitored"] == m
        ]
        for m in (True, False)
    }
    overhead_ratio = (
        statistics.median(overhead[True]) / statistics.median(overhead[False]) - 1
    ) * 100
    lines = [
        "# 3RLP performance pilot",
        "",
        "This report concerns short trajectory prefixes and computational performance. It does not establish scientific convergence or a full-trajectory speedup.",
        "",
        f"Host: {hardware['model']}, {hardware['cpu_count']} logical CPUs, {hardware['physical_ram_bytes'] / 2**30:.0f} GiB RAM. Separate installed Python/dependency environments are part of the comparison.",
        "",
        "## Measured analysis costs",
        "",
        "| Mode | Frames | Workers | Repeats | Median wall s | Range s | Tree peak MiB | OS peak MiB | Output MiB |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|",
    ]

    def memory(v):
        return f"{v / 2**20:.1f}" if v is not None else "unavailable"

    for r in rows:
        lines.append(
            f"| {r['mode']} | {r['frames']} | {r['workers']} | {r['repetitions']} | {r['wall_median_s']:.2f} | {r['wall_min_s']:.2f}–{r['wall_max_s']:.2f} | {memory(r['sampled_tree_peak_rss_bytes'])} | {memory(r['os_peak_rss_bytes'])} | {memory(r['output_bytes'])} |"
        )
    lines += [
        "",
        "## CPU use and stages",
        "",
        "| Mode at 1,001 frames, single worker | Frames/s | CPU s/frame | Peak processes | Peak threads |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in rows:
        if r["frames"] == 1001 and r["workers"] == "single":
            lines.append(
                f"| {r['mode']} | {r['frames_per_second']:.2f} | {r['cpu_seconds_per_frame']:.4f} | {r['peak_processes']} | {r['peak_threads']} |"
            )
    lines += [
        "",
        "Stage timings are inclusive; nested stages must not be added to their parents. External elapsed time is authoritative. Stage totals can differ from external elapsed time because timing boundaries and clocks differ.",
        "",
        "| Mode at 1,001 frames | Stage | Depth | Inclusive s |",
        "|---|---|---:|---:|",
    ]
    for r in stage_rows:
        if r["kind"] == "benchmark" and r["frames"] == 1001:
            lines.append(
                f"| {r['mode']} | {r['stage']} | {r['depth']} | {r['inclusive_seconds']:.3f} |"
            )
    lines += [
        "",
        "SSTMap's entropy-script timing includes its compiled subprocesses, but collecting entropy inputs also occurs within the site-calculation loop. The script timing alone is not the complete cost of obtaining entropy. These stage timings do not isolate equivalent entropy workloads across tools.",
        "",
        "## Cache reuse",
        "",
    ]
    for r in records:
        if r.get("kind") == "cached":
            elapsed = r.get("os_wall_seconds", r["wall_seconds"])
            reused = r.get("adapter", {}).get("result", {}).get("cache_reused")
            lines.append(
                f"HydraRank at {r['frames']} frames: {r['status']}, {elapsed:.2f} s, cache reuse verified: {reused}. This trial is excluded from fresh-analysis ratios and projections."
            )
    lines += ["", "## Validation and failures", ""]
    for r in records:
        if r.get("kind") == "preflight":
            error = r.get("adapter", {}).get("error", "")
            if r["exit_code"] == -11:
                error = "SIGSEGV in native voxel assignment; see ../GIST_BLOCKER.md and preserved diagnostic logs."
            lines.append(f"- {r['mode']}: {r['status']}. {error}")
    validation_path = args.campaign / "validation.json"
    if validation_path.exists():
        box = json.loads(validation_path.read_text()).get("gist_reference_box")
        if box:
            lines += [
                "",
                f"The intended 0.5 Å GIST box spans {box['lower_A']} to {box['upper_A']} Å, with dimensions {box['dimensions']} and volume {box['volume_A3']:.1f} Å³. An independent Cartesian count over {box['frames']} prepared frames found {box['reference_mean_waters']:.1f} oxygen atoms per frame on average (range {box['reference_min_waters']}–{box['reference_max_waters']}). These are reference box counts, not validated GIST voxel assignments; GIST sampled counts remain unavailable.",
            ]
    prepare = [r for r in records if r.get("kind") == "preparation"]
    if prepare:
        r = prepare[0]
        lines += [
            "",
            f"Preparation of the 1,001-frame raw and aligned prefixes plus validation took {r['wall_seconds']:.2f} s. Stage timings distinguish raw-prefix export, alignment initialization, and alignment/export. Raw-prefix export is a benchmarking convenience, not part of SSTMap's scientific preparation requirement.",
            "Preparation resource use is separate from analysis. The existing historical aligned trajectory was checked at representative frames in a common rigid coordinate frame.",
        ]
        alignment = sum(
            e["seconds"]
            for e in r.get("adapter", {}).get("stages", [])
            if e["name"] in ("alignment_initialization", "alignment_and_export")
        )
        # Match the complete prefix preparation to the largest pilot measurements.
        lines += [
            "",
            "## Preparation and combined cost",
            "",
            f"The measured alignment preparation stages for 1,001 frames total {alignment:.2f} s; this includes trajectory export and excludes raw-prefix export and validation. This is a single measurement, not a repeated preparation benchmark.",
            "",
            "| Mode, single worker | Analysis at 1,001 frames s | Alignment preparation s | Combined s |",
            "|---|---:|---:|---:|",
        ]
        for item in rows:
            if item["frames"] == 1001 and item["workers"] == "single":
                extra = 0 if item["mode"] == "hydrarank" else alignment
                lines.append(
                    f"| {item['mode']} | {item['wall_median_s']:.2f} | {extra:.2f} | {item['wall_median_s'] + extra:.2f} |"
                )
        lines += [
            "",
            "HydraRank already includes alignment. SSTMap initialization and Python startup are included in its analysis time; preparation-stage totals omit their own interpreter startup. The summed cost is therefore an explicitly defined estimate of the combined stages, not an independently measured end-to-end SSTMap command.",
        ]
    lines += [
        "",
        "## Projected full-study cost",
        "",
        "The following planning envelopes extrapolate from the largest completed single-worker prefix. They are not measured full-trajectory runtimes. Lower bounds scale wall time linearly with frame count; upper bounds use twice that projection. Nonlinear entropy/clustering costs, memory growth, and grid volume can invalidate these estimates.",
        "",
        "| Mode | Largest measured frames | Projected 25,001-frame run | Three-repeat frame-scaling campaign |",
        "|---|---:|---:|---:|",
    ]
    projections = []
    for mode in sorted({r["mode"] for r in rows}):
        selected = [
            r for r in rows if r["mode"] == mode and r["workers"] == "single" and r["padding"] == 5
        ]
        if not selected:
            continue
        r = max(selected, key=lambda x: x["frames"])
        full = r["wall_median_s"] * 25001 / r["frames"]
        campaign = (
            r["wall_median_s"] * sum((101, 501, 1001, 2501, 5001, 10001, 25001)) * 3 / r["frames"]
        )
        projections.append(
            dict(
                mode=mode,
                full_run_linear_seconds=full,
                full_run_planning_upper_seconds=2 * full,
                scaling_campaign_linear_seconds=campaign,
            )
        )
        lines.append(
            f"| {mode} | {r['frames']} | {full / 60:.1f}–{2 * full / 60:.1f} min | {campaign / 3600:.2f}–{2 * campaign / 3600:.2f} h |"
        )
    dump(out / "projections.json", projections)
    total = sum(r["scaling_campaign_linear_seconds"] for r in projections)
    lines += [
        "",
        f"Frame-scaling subtotal for validated modes: {total / 3600:.2f}–{2 * total / 3600:.2f} hours, excluding native-worker runs, cache trials, grid-region studies, preparation, and any blocked modes.",
        "Memory projections are not assumed linear. Full runs must retain the half-RAM safety limit and report any infeasible configurations.",
        "",
        "## Interpretation and limits",
        "",
        "- Historical full SSTMap HSA: 83 min 55 s for 25,001 frames; clustering 88.18 s; site calculations 4,935.28 s. These are separate historical measurements without resource monitoring or repetitions.",
        "- HSA and HydraRank use related site representations; the reduced HSA workflow is overlapping functionality, not identical functionality. Native clustering may give different site counts.",
        "- GIST targets the ligand bounding box plus padding, including more solvent than the 5 Å ligand-distance region. Both modes crashed before output metadata was finalized; sampled grid-water counts are unavailable and must be measured after the native blocker is resolved.",
        "- Analysis-only SSTMap measurements exclude alignment; preparation must be added for a raw-input comparison. HydraRank includes its local alignment during preprocessing.",
        "- Single-worker trials override HydraRank KD-tree workers and limit numerical libraries. Native trials remove those limits; actual observed thread peaks are reported.",
        "- OS peak RSS and sampled aggregate process-tree RSS measure different quantities; the latter can miss short-lived peaks between 0.5 s samples and double-count shared memory.",
        "- Output storage excludes the original simulation inputs and includes logs and tool-local caches, but excludes resource-monitor JSON written after the timed process exits. Sampled peak disk usage can miss temporary files removed between scans.",
        "- Shared trajectory-index caches were generated during input preparation and preflight; they are retained across trials. Fresh-cache comparisons refer to observations and per-trial runtime caches, not rebuilding shared trajectory indices.",
        "- CPU times include waited descendants. Page faults and block-operation counters are OS-specific; they do not measure physical storage throughput.",
        "- Each trial uses a new directory and a fresh process. Input prefixes are read before timing; no disk-cold or fully idle-host claim is made. Diagnostic and report-development activity occurred during parts of the pilot; repeat full-study measurements on an otherwise idle host.",
        f"- Monitor microbenchmark median elapsed difference: {overhead_ratio:+.1f}% across three paired runs. Scheduling noise and startup cost limit this estimate; inspect raw trials rather than treating it as a correction factor.",
        "- The installed GIST voxel-assignment extension hardcodes 0.5 Å. The 0.25/1.0 Å resolution study is blocked pending a documented, validated correction.",
        "- Cached HydraRank trials are labelled separately and excluded from primary summaries and projections.",
        "",
        "The full study remains pending user review. Do not replace these planning estimates with manuscript claims of full-trajectory speed.",
    ]
    lines += [
        "",
        "## Figures",
        "",
        "![Wall time versus frame count](figures/time.png)",
        "",
        "![Peak sampled memory versus frame count](figures/memory.png)",
        "",
        "![Output storage versus frame count](figures/storage.png)",
    ]
    (out / "report.md").write_text("\n".join(lines) + "\n")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = out / "figures"
    figdir.mkdir(exist_ok=True)
    for metric, label, name in [
        ("wall_median_s", "Wall time (s)", "time"),
        ("sampled_tree_peak_rss_bytes", "Sampled tree peak RSS (MiB)", "memory"),
        ("output_bytes", "Output storage (MiB)", "storage"),
    ]:
        fig, ax = plt.subplots(figsize=(8, 5))
        for mode in sorted({r["mode"] for r in rows}):
            selected = sorted(
                [r for r in rows if r["mode"] == mode and r["workers"] == "single"],
                key=lambda r: r["frames"],
            )
            x = [r["frames"] for r in selected]
            y = [r[metric] / (2**20 if name != "time" else 1) for r in selected]
            ax.plot(x, y, "o-", label=mode)
            if name == "time":
                ax.fill_between(
                    x,
                    [r["wall_min_s"] for r in selected],
                    [r["wall_max_s"] for r in selected],
                    alpha=0.12,
                )
        ax.set(
            xlabel="Trajectory frames (2 ps spacing)",
            ylabel=label,
            title="3RLP pilot, single-worker analysis",
        )
        ax.legend()
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(figdir / f"{name}.png", dpi=180)
        plt.close(fig)
    print(out / "report.md")


if __name__ == "__main__":
    main()
