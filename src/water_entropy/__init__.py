"""Light-weight hydration-site analysis around a ligand."""

from importlib.metadata import PackageNotFoundError, version

from water_entropy.config import PreprocessConfig
from water_entropy.exceptions import (
    AmbiguousSelectionError,
    EmptySelectionError,
    WaterEntropyError,
    WaterModelError,
)

try:
    __version__ = version("water-entropy")
except PackageNotFoundError:  # source tree imported without installation
    __version__ = "0+unknown"

__all__ = [
    "AmbiguousSelectionError",
    "EmptySelectionError",
    "PreprocessConfig",
    "WaterEntropyError",
    "WaterModelError",
    "__version__",
]
