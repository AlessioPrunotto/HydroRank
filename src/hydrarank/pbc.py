"""Periodic-boundary treatment.

The transformations are applied lazily by MDAnalysis while iterating, so no
intermediate trajectory is written. The order matters and is not interchangeable:

1. make the solute whole and track its atoms across periodic images;
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
from MDAnalysis.lib._cutil import make_whole
from MDAnalysis.lib.distances import apply_PBC, minimize_vectors
from MDAnalysis.transformations import center_in_box
from MDAnalysis.transformations.base import TransformationBase

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


class _FastUnwrap(TransformationBase):
    """Make molecules whole without repeating a bond-graph walk every frame.

    MDAnalysis' general-purpose ``unwrap`` transformation calls ``make_whole``
    for every molecular fragment on every frame.  For a protein this graph walk
    dominates the runtime.  Molecular-dynamics coordinates are continuous
    between saved frames, so after making the first frame whole we can recover
    subsequent images with one vectorised minimum-image operation.
    """

    parallelizable = False

    def __init__(self, atoms: AtomGroup):
        super().__init__(parallelizable=False)
        self.atoms = atoms
        self._indices = atoms.indices
        self._fragments = tuple(atoms.fragments)
        self._previous: np.ndarray | None = None
        self._previous_frame: int | None = None

    def _transform(self, ts):
        # A backwards seek starts a new trajectory pass.  Rebuild that frame
        # from bonds; forward iteration (including a configured stride) uses
        # the much cheaper no-jump update.
        if self._previous is None or (
            self._previous_frame is not None and ts.frame <= self._previous_frame
        ):
            for fragment in self._fragments:
                make_whole(fragment)
            current = ts.positions[self._indices].copy()
        else:
            raw = ts.positions[self._indices]
            current = self._previous + minimize_vectors(raw - self._previous, ts.dimensions)
            ts.positions[self._indices] = current

        self._previous = current
        self._previous_frame = ts.frame
        return ts


class _FastResidueWrap(TransformationBase):
    """Vectorised equivalent of ``AtomGroup.wrap(compound='residues')``."""

    def __init__(self, atoms: AtomGroup):
        super().__init__()
        self._indices = atoms.indices
        _, self._residue_inverse = np.unique(atoms.resindices, return_inverse=True)
        self._masses = atoms.masses.astype(np.float64, copy=False)
        self._total_masses = np.bincount(self._residue_inverse, weights=self._masses)
        if np.any(np.isclose(self._total_masses, 0.0)):
            raise HydraRankError("cannot wrap a residue whose total mass is zero")

    def _transform(self, ts):
        if self._indices.size == 0:
            return ts
        positions = ts.positions[self._indices]
        weighted = positions * self._masses[:, None]
        centers = np.column_stack(
            [np.bincount(self._residue_inverse, weights=weighted[:, axis]) for axis in range(3)]
        )
        centers /= self._total_masses[:, None]
        shifts = apply_PBC(centers.astype(np.float32), ts.dimensions) - centers
        ts.positions[self._indices] = positions + shifts[self._residue_inverse]
        return ts


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
        _FastUnwrap(groups.solute),
        center_in_box(groups.center, center="geometry"),
        _FastResidueWrap(groups.mobile),
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
