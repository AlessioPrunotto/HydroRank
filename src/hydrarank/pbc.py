"""Periodic-boundary treatment.

The transformations are applied lazily by MDAnalysis while iterating, so no
intermediate trajectory is written. The order matters and is not interchangeable:

1. ``unwrap`` the solute so protein and ligand are whole molecules again;
2. ``center_in_box`` on the ligand, so the pocket sits at the middle of the box;
3. ``wrap`` everything else *by residue*, which pulls the periodic images of the
   nearby waters next to the pocket while keeping each water molecule intact.

After this, plain Euclidean distances are valid inside the pocket and the rest of
the pipeline can ignore periodicity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.transformations import center_in_box, unwrap, wrap

from hydrarank.exceptions import HydraRankError, TrajectoryError

#: A whole molecule cannot contain a covalent bond longer than this (angstrom).
MAX_BOND_LENGTH = 3.0


@dataclass(frozen=True)
class PBCGroups:
    """The three atom groups the transformation stack acts on."""

    solute: AtomGroup
    """Molecules that must be made whole (protein + ligand)."""

    mobile: AtomGroup
    """Everything wrapped by residue (water, ions, lipids, ...)."""

    center: AtomGroup
    """Group placed at the centre of the box, normally the ligand."""


def build_pbc_groups(universe: Universe, ligand: AtomGroup) -> PBCGroups:
    solute = universe.select_atoms("protein or nucleic") | ligand
    mobile = universe.atoms - solute
    return PBCGroups(solute=solute, mobile=mobile, center=ligand)


def apply_pbc_transformations(
    universe: Universe, groups: PBCGroups, guess_bonds: bool = False
) -> Universe:
    """Attach the unwrap/centre/wrap stack to ``universe.trajectory``.

    Returns the same universe, mutated. Can only be called once per trajectory.
    """
    if universe.trajectory.transformations:
        raise HydraRankError("transformations have already been applied to this trajectory")
    _require_box(universe)
    _require_bonds(universe, groups.solute, guess_bonds)

    universe.trajectory.add_transformations(
        unwrap(groups.solute),
        center_in_box(groups.center, center="geometry"),
        wrap(groups.mobile, compound="residues"),
    )
    return universe


def validate_cutoff(universe: Universe, cutoff: float, name: str = "cutoff") -> None:
    """Ensure a selection cutoff is compatible with the minimum-image convention."""
    box = universe.dimensions
    if box is None:
        raise TrajectoryError("no box information; cannot validate the cutoff")
    half_box = float(np.min(box[:3])) / 2.0
    if cutoff >= half_box:
        raise HydraRankError(
            f"{name} {cutoff:g} A is not smaller than half the shortest box vector "
            f"({half_box:g} A); the minimum-image convention would be violated"
        )


def max_bond_length(atoms: AtomGroup) -> float:
    """Longest bond in ``atoms`` on the current frame, ignoring periodic images.

    A value above :data:`MAX_BOND_LENGTH` means the molecule is still split across
    the periodic boundary.
    """
    bonds = atoms.bonds
    if len(bonds) == 0:
        return 0.0
    return float(np.max(bonds.values()))


def assert_solute_whole(atoms: AtomGroup, tolerance: float = MAX_BOND_LENGTH) -> None:
    longest = max_bond_length(atoms)
    if longest > tolerance:
        raise HydraRankError(
            f"solute is not whole: longest bond is {longest:.2f} A (> {tolerance:.2f} A). "
            "The unwrap step failed, most likely because of missing or wrong bond information."
        )


def _require_box(universe: Universe) -> None:
    box = universe.dimensions
    if box is None or not np.all(box[:3] > 0):
        raise TrajectoryError(
            "no periodic box information; PBC treatment requires box dimensions in the trajectory"
        )


def _require_bonds(universe: Universe, solute: AtomGroup, guess_bonds: bool) -> None:
    if hasattr(universe.atoms, "bonds") and len(solute.bonds) > 0:
        return
    if not guess_bonds:
        raise TrajectoryError(
            "the topology carries no bonds, so molecules cannot be made whole. "
            "Use a topology that stores connectivity (PSF, TPR, prmtop) or set "
            "guess_bonds=True to infer them from distances."
        )
    solute.guess_bonds()
