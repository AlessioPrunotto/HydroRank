"""Light-weight hydration-site analysis around a ligand."""

from importlib.metadata import PackageNotFoundError, version

from hydrarank.config import PreprocessConfig
from hydrarank.exceptions import (
    AmbiguousSelectionError,
    EmptySelectionError,
    HydraRankError,
    WaterModelError,
)

try:
    __version__ = version("hydrarank")
except PackageNotFoundError:  # source tree imported without installation
    __version__ = "0+unknown"

from hydrarank.workflow import AnalyseResult, run_analysis  # noqa: E402

__all__ = [
    "AmbiguousSelectionError",
    "AnalyseResult",
    "EmptySelectionError",
    "PreprocessConfig",
    "HydraRankError",
    "WaterModelError",
    "__version__",
    "run_analysis",
]
