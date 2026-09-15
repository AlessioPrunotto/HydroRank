"""Exception hierarchy for water-entropy."""


class WaterEntropyError(Exception):
    """Base class for all errors raised by water-entropy."""


class EmptySelectionError(WaterEntropyError):
    """A selection that must not be empty matched zero atoms."""


class AmbiguousSelectionError(WaterEntropyError):
    """Automatic detection found several equally plausible candidates."""


class WaterModelError(WaterEntropyError):
    """The water residues are inconsistent or of an unsupported geometry."""


class TrajectoryError(WaterEntropyError):
    """The trajectory is unusable for hydration-site analysis."""
