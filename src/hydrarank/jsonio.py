"""Strict JSON serialization helpers for scientific results."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np


def dumps(value: Any, **kwargs) -> str:
    """Serialize scientific data as strict JSON, mapping non-finite values to null."""
    return json.dumps(to_json_compatible(value), allow_nan=False, **kwargs)


def to_json_compatible(value: Any) -> Any:
    """Recursively convert paths, NumPy scalars, and non-finite floats for JSON."""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): to_json_compatible(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
        return [to_json_compatible(item) for item in value]
    return value
