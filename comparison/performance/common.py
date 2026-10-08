"""Shared benchmark helpers; uses only the standard library."""

import hashlib
import json
import time
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MODES = ("hydrarank", "hsa_full", "hsa_reduced", "gist_full", "gist_entropy")


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


class Timings:
    def __init__(self):
        self.events = []
        self.depth = 0

    @contextmanager
    def stage(self, name):
        start = time.perf_counter()
        depth = self.depth
        self.depth += 1
        try:
            yield
        finally:
            self.depth -= 1
            self.events.append(dict(name=name, seconds=time.perf_counter() - start, depth=depth))

    def wrap(self, obj, name, label=None):
        original = getattr(obj, name)

        def wrapped(*args, **kwargs):
            with self.stage(label or name):
                return original(*args, **kwargs)

        setattr(obj, name, wrapped)
