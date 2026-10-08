"""External macOS process-tree monitor and OS timing, including descendants."""

import ctypes
import os
import re
import signal
import subprocess
import time
from contextlib import suppress
from pathlib import Path

from common import dump


class TaskInfo(ctypes.Structure):
    _fields_ = [
        (n, ctypes.c_uint64)
        for n in ("virtual", "resident", "user", "system", "threads_user", "threads_system")
    ]
    _fields_ += [
        (n, ctypes.c_int32)
        for n in (
            "policy",
            "faults",
            "pageins",
            "cow_faults",
            "sent",
            "received",
            "mach_calls",
            "unix_calls",
            "switches",
            "threads",
            "running",
            "priority",
        )
    ]


def physical_ram():
    return int(subprocess.check_output(["/usr/sbin/sysctl", "-n", "hw.memsize"], text=True))


def snapshot(root_pid):
    rows = subprocess.check_output(["/bin/ps", "-axo", "pid=,ppid=,rss="], text=True)
    processes = {
        int(p): (int(parent), int(rss) * 1024)
        for p, parent, rss in (row.split() for row in rows.splitlines())
    }
    selected = {root_pid}
    while True:
        expanded = selected | {p for p, (parent, _) in processes.items() if parent in selected}
        if expanded == selected:
            break
        selected = expanded
    lib = ctypes.CDLL("/usr/lib/libproc.dylib")
    threads = 0
    known = 0
    for pid in selected:
        info = TaskInfo()
        if lib.proc_pidinfo(pid, 4, 0, ctypes.byref(info), ctypes.sizeof(info)) == ctypes.sizeof(
            info
        ):
            threads += info.threads
            known += 1
    return dict(
        rss_bytes=sum(processes[p][1] for p in selected if p in processes),
        processes=sum(p in processes for p in selected),
        threads=threads if known else None,
    )


def normalize_rss(value, system):
    # Darwin's ru_maxrss/time -l uses bytes; Linux uses KiB.
    return int(value) if system == "Darwin" else int(value) * 1024


def parse_time(text):
    values = {}
    match = re.search(r"([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys", text)
    if match:
        values.update(
            os_wall_seconds=float(match[1]),
            user_seconds=float(match[2]),
            system_seconds=float(match[3]),
        )
    for label, key in [
        ("maximum resident set size", "os_peak_rss_bytes"),
        ("page faults", "page_faults"),
        ("page reclaims", "page_reclaims"),
        ("block input operations", "block_inputs"),
        ("block output operations", "block_outputs"),
    ]:
        match = re.search(r"(\d+)\s+" + label, text)
        if match:
            values[key] = int(match[1])
    return values


def terminate_group(pid):
    try:
        os.killpg(pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    time.sleep(0.1)
    with suppress(ProcessLookupError, PermissionError):
        os.killpg(pid, signal.SIGKILL)


def run(command, destination, env=None, timeout=1200, memory_limit=None, monitored=True):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    status = "complete"
    samples = []
    high_disk = 0
    with (
        (destination / "stdout.log").open("w") as out,
        (destination / "stderr.log").open("w") as err,
    ):
        process = subprocess.Popen(
            ["/usr/bin/time", "-l", *map(str, command)],
            cwd=destination,
            env=env,
            stdout=out,
            stderr=err,
            start_new_session=True,
        )
        while process.poll() is None:
            elapsed = time.monotonic() - start
            if monitored:
                sample = snapshot(process.pid)
                sample["elapsed_seconds"] = elapsed
                samples.append(sample)
                if memory_limit is not None and sample["rss_bytes"] > memory_limit:
                    status = "memory_limit"
                    terminate_group(process.pid)
                    break
            if elapsed > timeout:
                status = "timeout"
                terminate_group(process.pid)
                break
            # File scans are infrequent so monitoring remains lightweight.
            if len(samples) % 10 == 0:
                high_disk = max(high_disk, tree_size(destination))
            time.sleep(0.5 if monitored else 0.05)
        returncode = process.wait()
    if returncode and status == "complete":
        status = "failed"
    elapsed = time.monotonic() - start
    disk = tree_size(destination)
    metrics = dict(
        status=status,
        exit_code=returncode,
        wall_seconds=elapsed,
        command=list(map(str, command)),
        monitored=monitored,
        sampled_tree_peak_rss_bytes=max((s["rss_bytes"] for s in samples), default=None),
        peak_processes=max((s["processes"] for s in samples), default=None),
        peak_threads=max((s["threads"] for s in samples if s["threads"] is not None), default=None),
        output_bytes=disk,
        sampled_peak_disk_bytes=max(disk, high_disk),
        output_files=sum(p.is_file() for p in destination.rglob("*")),
    )
    metrics.update(parse_time((destination / "stderr.log").read_text(errors="replace")))
    dump(destination / "resources.json", metrics)
    dump(destination / "samples.json", samples)
    return metrics


def tree_size(path):
    return sum(p.stat().st_size for p in Path(path).rglob("*") if p.is_file())
