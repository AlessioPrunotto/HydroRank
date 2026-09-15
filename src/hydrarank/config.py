"""Configuration objects for the preprocessing stage."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from hydrarank.exceptions import HydraRankError

#: Residue names used for water by the common force fields / MD engines.
#: Deliberately explicit: MDAnalysis' and MDTraj's built-in ``water`` keywords
#: do not know about OPC, which silently yields empty selections.
DEFAULT_WATER_RESNAMES: tuple[str, ...] = (
    "HOH",
    "WAT",
    "SOL",
    "H2O",
    "TIP",
    "TIP2",
    "TIP3",
    "TIP3P",
    "TIP4",
    "TIP4P",
    "TIP5",
    "TIP5P",
    "T3P",
    "T4P",
    "SPC",
    "SPCE",
    "OPC",
    "OPC3",
)

#: Monoatomic ions that must never be mistaken for a ligand.
DEFAULT_ION_RESNAMES: tuple[str, ...] = (
    "NA",
    "NA+",
    "SOD",
    "K",
    "K+",
    "POT",
    "LI",
    "LI+",
    "RB",
    "CS",
    "CL",
    "CL-",
    "CLA",
    "BR",
    "IOD",
    "F",
    "MG",
    "MG2+",
    "CA",
    "CAL",
    "ZN",
    "ZN2",
    "MN",
    "FE",
    "CU",
    "NI",
    "CO",
)


@dataclass(frozen=True)
class PreprocessConfig:
    """Everything needed to turn a raw MD trajectory into water observations.

    Distances are in angstrom (MDAnalysis convention) throughout.
    """

    topology: Path
    trajectory: list[Path] = field(default_factory=list)

    ligand_selection: str | None = None
    """MDAnalysis selection string for the ligand. ``None`` triggers auto-detection."""

    water_resnames: tuple[str, ...] = DEFAULT_WATER_RESNAMES
    ion_resnames: tuple[str, ...] = DEFAULT_ION_RESNAMES

    start: int = 0
    stop: int | None = None
    step: int = 1

    water_cutoff: float = 5.0
    """Max distance from any ligand heavy atom for a water oxygen to be retained."""

    pocket_cutoff: float = 8.0
    """Max distance from the ligand for a protein residue to be considered pocket."""

    min_ligand_heavy_atoms: int = 8
    """Auto-detection ignores non-protein residues smaller than this."""

    site_radius: float = 1.0
    """Radius of a hydration site, in angstrom."""

    density_factor: float = 2.0
    """A site must be at least this many times denser than bulk water to be kept."""

    temperature: float = 300.0
    """Temperature used to turn entropies into free energies, in kelvin."""

    hbond_distance: float = 3.5
    """Donor-acceptor distance cutoff for a hydrogen bond, in angstrom."""

    hbond_angle: float = 130.0
    """Minimum donor-H...acceptor angle for a hydrogen bond, in degrees."""

    enclosure_radius: float = 5.0
    """Radius used to count the solute heavy atoms around a water, in angstrom."""

    hbond_penalty: float = 1.0
    """Assumed free-energy cost of breaking one water-solute hydrogen bond, kcal/mol."""

    max_gap: int = 0
    """Frames of absence tolerated within a single residence episode."""

    output_dir: Path = Path("output")

    def __post_init__(self) -> None:
        if self.step < 1:
            raise HydraRankError(f"step must be >= 1, got {self.step}")
        if self.start < 0:
            raise HydraRankError(f"start must be >= 0, got {self.start}")
        if self.stop is not None and self.stop <= self.start:
            raise HydraRankError(
                f"stop must be greater than start (got start={self.start}, stop={self.stop})"
            )
        if self.min_ligand_heavy_atoms < 1:
            raise HydraRankError(
                f"min_ligand_heavy_atoms must be >= 1, got {self.min_ligand_heavy_atoms}"
            )
        if self.temperature <= 0:
            raise HydraRankError(f"temperature must be > 0, got {self.temperature}")
        if self.max_gap < 0:
            raise HydraRankError(f"max_gap must be >= 0, got {self.max_gap}")
        if self.hbond_distance <= 0:
            raise HydraRankError(f"hbond_distance must be > 0, got {self.hbond_distance}")
        if not 0.0 <= self.hbond_angle <= 180.0:
            raise HydraRankError(f"hbond_angle must be in [0, 180], got {self.hbond_angle}")
        if self.enclosure_radius <= 0:
            raise HydraRankError(f"enclosure_radius must be > 0, got {self.enclosure_radius}")
        if self.hbond_penalty < 0:
            raise HydraRankError(f"hbond_penalty must be >= 0, got {self.hbond_penalty}")
        if self.water_cutoff <= 0:
            raise HydraRankError(f"water_cutoff must be > 0, got {self.water_cutoff}")
        if self.site_radius <= 0:
            raise HydraRankError(f"site_radius must be > 0, got {self.site_radius}")
        if self.density_factor <= 0:
            raise HydraRankError(f"density_factor must be > 0, got {self.density_factor}")
        if self.pocket_cutoff < self.water_cutoff:
            raise HydraRankError(
                "pocket_cutoff must be >= water_cutoff so that every retained water "
                f"is surrounded by pocket atoms (got {self.pocket_cutoff} < {self.water_cutoff})"
            )

    @classmethod
    def from_yaml(cls, path: str | Path) -> PreprocessConfig:
        path = Path(path)
        with path.open() as handle:
            raw = yaml.safe_load(handle) or {}
        if not isinstance(raw, dict):
            raise HydraRankError(f"{path} must contain a YAML mapping")
        return cls.from_dict(raw, base_dir=path.parent)

    @classmethod
    def from_dict(cls, raw: dict[str, Any], base_dir: Path | None = None) -> PreprocessConfig:
        known = {f.name for f in dataclasses.fields(cls)}
        unknown = set(raw) - known
        if unknown:
            raise HydraRankError(f"unknown configuration keys: {sorted(unknown)}")

        data = dict(raw)
        base = Path(base_dir) if base_dir is not None else Path()

        if "topology" not in data:
            raise HydraRankError("configuration must define 'topology'")
        data["topology"] = _resolve(data["topology"], base)

        traj = data.get("trajectory") or []
        if isinstance(traj, str | Path):
            traj = [traj]
        data["trajectory"] = [_resolve(item, base) for item in traj]

        if "output_dir" in data:
            data["output_dir"] = _resolve(data["output_dir"], base)
        for key in ("water_resnames", "ion_resnames"):
            if key in data:
                data[key] = tuple(str(name).upper() for name in data[key])

        return cls(**data)

    def to_dict(self) -> dict[str, Any]:
        out = dataclasses.asdict(self)
        out["topology"] = str(self.topology)
        out["trajectory"] = [str(p) for p in self.trajectory]
        out["output_dir"] = str(self.output_dir)
        out["water_resnames"] = list(self.water_resnames)
        out["ion_resnames"] = list(self.ion_resnames)
        return out

    def to_yaml(self, path: str | Path) -> None:
        with Path(path).open("w") as handle:
            yaml.safe_dump(self.to_dict(), handle, sort_keys=False)


def _resolve(value: str | Path, base: Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path)
