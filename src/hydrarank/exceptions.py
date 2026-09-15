"""Exception hierarchy for hydrarank."""


class HydraRankError(Exception):
    """Base class for all errors raised by hydrarank."""


class EmptySelectionError(HydraRankError):
    """A selection that must not be empty matched zero atoms."""


class AmbiguousSelectionError(HydraRankError):
    """Automatic detection found several equally plausible candidates."""


class WaterModelError(HydraRankError):
    """The water residues are inconsistent or of an unsupported geometry."""


class TrajectoryError(HydraRankError):
    """The trajectory is unusable for hydration-site analysis."""
