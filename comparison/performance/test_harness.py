"""Run these tests outside the sandbox: process inspection is required on macOS."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from monitor import normalize_rss, parse_time, run


def test_units_and_os_time():
    assert normalize_rss(100, "Darwin") == 100
    assert normalize_rss(100, "Linux") == 102400
    parsed = parse_time(
        "1.25 real 0.20 user 0.10 sys\n 1234 maximum resident set size\n 3 page faults"
    )
    assert parsed["os_peak_rss_bytes"] == 1234
    assert parsed["user_seconds"] == 0.20
    assert parsed["page_faults"] == 3


def test_completion_failure_and_isolation(tmp_path):
    result = run([sys.executable, "-c", 'print("ok")'], tmp_path / "success", monitored=False)
    assert result["status"] == "complete"
    assert result["os_peak_rss_bytes"] > 0
    assert (tmp_path / "success/stdout.log").read_text().strip() == "ok"
    with pytest.raises(FileExistsError):
        run([sys.executable, "-c", "pass"], tmp_path / "success", monitored=False)
    result = run(
        [sys.executable, "-c", "raise SystemExit(3)"], tmp_path / "failure", monitored=False
    )
    assert result["status"] == "failed"
    assert result["exit_code"] == 3


def test_children_and_timeout(tmp_path):
    command = [
        sys.executable,
        "-c",
        'import subprocess,sys; p=subprocess.Popen([sys.executable,"-c",'
        '"import time; b=bytearray(30_000_000); time.sleep(5)"]); p.wait()',
    ]
    result = run(command, tmp_path / "child", timeout=1.2)
    assert result["status"] == "timeout"
    assert result["peak_processes"] >= 3  # time, interpreter, child
    assert result["sampled_tree_peak_rss_bytes"] >= 30_000_000
    assert result["peak_threads"] >= 2


def test_memory_limit(tmp_path):
    result = run(
        [sys.executable, "-c", "import time; b=bytearray(30_000_000); time.sleep(5)"],
        tmp_path / "limit",
        timeout=4,
        memory_limit=1,
    )
    assert result["status"] == "memory_limit"
    assert json.loads((tmp_path / "limit/resources.json").read_text())["status"] == "memory_limit"
