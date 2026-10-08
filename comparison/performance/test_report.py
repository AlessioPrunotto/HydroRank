"""Aggregation must never mix cached or failed trials into speed comparisons."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from report import summarize


def test_summary_excludes_failed_cached_and_native_groups():
    base = dict(
        kind="benchmark",
        mode="hydrarank",
        frames=101,
        workers="single",
        spacing=0.5,
        padding=5,
        status="complete",
        wall_seconds=10,
        os_wall_seconds=9,
    )
    records = [
        base,
        {**base, "os_wall_seconds": 11},
        {**base, "kind": "cached", "os_wall_seconds": 1},
        {**base, "status": "failed", "os_wall_seconds": 1},
        {**base, "workers": "native", "os_wall_seconds": 4},
    ]
    result = summarize(records)
    single = next(r for r in result if r["workers"] == "single")
    assert single["wall_median_s"] == 10
    assert single["wall_min_s"] == 9
    assert single["wall_max_s"] == 11
    assert single["repetitions"] == 2
    assert len(result) == 2
