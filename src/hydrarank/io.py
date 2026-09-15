"""Loading and validating an MD system, and reporting what was found."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from MDAnalysis import Universe

from hydrarank.config import PreprocessConfig
from hydrarank.exceptions import TrajectoryError
from hydrarank.selections import (
    WaterTopology,
    analyse_water_topology,
    classify_atoms,
    select_ligand,
    select_water,
)


def load_universe(config: PreprocessConfig) -> Universe:
    """Build a ``Universe`` from the configured files, with early validation."""
    topology = Path(config.topology)
    if not topology.is_file():
        raise FileNotFoundError(f"topology not found: {topology}")

    trajectories = [Path(p) for p in config.trajectory]
    for path in trajectories:
        if not path.is_file():
            raise FileNotFoundError(f"trajectory not found: {path}")

    universe = (
        Universe(str(topology), *[str(p) for p in trajectories])
        if trajectories
        else Universe(str(topology))
    )
    _validate(universe)
    return universe


def frame_slice(config: PreprocessConfig, n_frames: int) -> slice:
    """Translate the configured frame range into a slice, with bounds checking."""
    stop = n_frames if config.stop is None else min(config.stop, n_frames)
    if config.start >= stop:
        raise TrajectoryError(
            f"empty frame range: start={config.start}, stop={stop}, "
            f"trajectory has {n_frames} frames"
        )
    return slice(config.start, stop, config.step)


@dataclass
class SystemReport:
    """Human-readable summary of what the selections actually matched."""

    n_atoms: int
    n_frames: int
    dt_ps: float
    n_selected_frames: int
    box: tuple[float, float, float, float, float, float] | None
    is_orthorhombic: bool
    has_bonds: bool
    residue_counts: dict[str, int]
    n_waters: int
    water_model: str
    water_resname: str
    ligand_label: str
    n_ligand_atoms: int
    n_ligand_heavy_atoms: int
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_atoms": self.n_atoms,
            "n_frames": self.n_frames,
            "dt_ps": self.dt_ps,
            "n_selected_frames": self.n_selected_frames,
            "box": list(self.box) if self.box is not None else None,
            "is_orthorhombic": self.is_orthorhombic,
            "has_bonds": self.has_bonds,
            "residue_counts": self.residue_counts,
            "n_waters": self.n_waters,
            "water_model": self.water_model,
            "water_resname": self.water_resname,
            "ligand_label": self.ligand_label,
            "n_ligand_atoms": self.n_ligand_atoms,
            "n_ligand_heavy_atoms": self.n_ligand_heavy_atoms,
            "warnings": self.warnings,
        }

    def __str__(self) -> str:
        lines = [
            "System",
            f"  atoms              : {self.n_atoms}",
            f"  frames             : {self.n_frames} (dt = {self.dt_ps:g} ps)",
            f"  frames selected    : {self.n_selected_frames}",
            f"  box                : {self._box_str()}",
            f"  bonds in topology  : {'yes' if self.has_bonds else 'no'}",
            "Water",
            f"  molecules          : {self.n_waters} (resname {self.water_resname})",
            f"  model              : {self.water_model}",
            "Ligand",
            f"  selection          : {self.ligand_label}",
            f"  atoms              : {self.n_ligand_atoms} ({self.n_ligand_heavy_atoms} heavy)",
            "Residue composition",
        ]
        lines += [f"  {name:<8} {count}" for name, count in self.residue_counts.items()]
        if self.warnings:
            lines.append("Warnings")
            lines += [f"  - {message}" for message in self.warnings]
        return "\n".join(lines)

    def _box_str(self) -> str:
        if self.box is None:
            return "absent"
        a, b, c, alpha, beta, gamma = self.box
        shape = "orthorhombic" if self.is_orthorhombic else "triclinic"
        return f"{a:.2f} x {b:.2f} x {c:.2f} A, angles {alpha:.1f}/{beta:.1f}/{gamma:.1f} ({shape})"


def describe_system(
    universe: Universe, config: PreprocessConfig, top_n_residues: int = 12
) -> tuple[SystemReport, WaterTopology]:
    """Run the selections and summarise them; raises if a mandatory one is empty."""
    water = select_water(universe, config.water_resnames)
    water_topology = analyse_water_topology(water)
    ligand = select_ligand(
        universe,
        config.ligand_selection,
        water_resnames=config.water_resnames,
        ion_resnames=config.ion_resnames,
        min_heavy_atoms=config.min_ligand_heavy_atoms,
    )

    n_frames = universe.trajectory.n_frames
    frames = frame_slice(config, n_frames)
    n_selected = len(range(*frames.indices(n_frames)))

    dimensions = universe.dimensions
    has_box = dimensions is not None and bool(np.all(dimensions[:3] > 0))
    box = tuple(float(x) for x in dimensions) if has_box else None
    is_ortho = bool(box is not None and np.allclose(box[3:], 90.0))

    dt = float(getattr(universe.trajectory, "dt", 0.0) or 0.0)
    messages = _collect_warnings(universe, dt, config.step, n_selected)
    if not has_box:
        messages.insert(0, "no periodic box information; PBC-aware preprocessing is impossible")

    ligand_resnames = sorted(set(ligand.residues.resnames))
    label = config.ligand_selection or f"auto: resname {' '.join(ligand_resnames)}"

    report = SystemReport(
        n_atoms=universe.atoms.n_atoms,
        n_frames=n_frames,
        dt_ps=dt,
        n_selected_frames=n_selected,
        box=box,
        is_orthorhombic=is_ortho,
        has_bonds=_has_bonds(universe),
        residue_counts=_residue_counts(universe, top_n_residues),
        n_waters=water_topology.n_waters,
        water_model=water_topology.model,
        water_resname=" ".join(sorted(set(water.residues.resnames))),
        ligand_label=label,
        n_ligand_atoms=ligand.n_atoms,
        n_ligand_heavy_atoms=_n_heavy(ligand),
        warnings=messages,
    )
    return report, water_topology


def _validate(universe: Universe) -> None:
    if universe.atoms.n_atoms == 0:
        raise TrajectoryError("topology contains no atoms")
    try:
        trajectory = universe.trajectory
    except AttributeError:
        raise TrajectoryError(
            "no coordinates available: the topology file carries none, so a trajectory or "
            "coordinate file must be supplied"
        ) from None
    if trajectory.n_frames == 0:
        raise TrajectoryError("trajectory contains no frames")


def _collect_warnings(universe: Universe, dt: float, step: int, n_selected: int) -> list[str]:
    messages: list[str] = []
    effective_dt = dt * step
    if dt <= 0:
        messages.append("frame time step is unknown; residence times cannot be expressed in ps")
    elif effective_dt > 20.0:
        messages.append(
            f"effective frame spacing is {effective_dt:g} ps, which is longer than typical water "
            "residence times (1-100 ps); persistence and entropy estimates will be unreliable"
        )
    if n_selected < 500:
        messages.append(
            f"only {n_selected} frames selected; per-site statistics will be noisy "
            "(a few thousand frames is a reasonable target)"
        )
    if not _has_bonds(universe):
        messages.append(
            "topology carries no bond information; molecules cannot be made whole across "
            "periodic boundaries without guessing bonds"
        )
    return messages


def _has_bonds(universe: Universe) -> bool:
    return hasattr(universe.atoms, "bonds") and len(universe.atoms.bonds) > 0


def _residue_counts(universe: Universe, top_n: int) -> dict[str, int]:
    counts = Counter(universe.residues.resnames)
    return dict(counts.most_common(top_n))


def _n_heavy(atoms) -> int:
    roles = classify_atoms(atoms)
    return int(np.count_nonzero((roles != "H") & (roles != "M")))
