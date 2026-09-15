"""Superposition of every frame onto a common binding-site reference frame.

Hydration sites only exist in a frame where the pocket does not move, so all
coordinates are rigid-body fitted onto a reference. The fit uses the binding-site
protein backbone rather than the whole protein (domain motions would smear the
pocket) or the ligand alone (unstable for small or symmetric ligands, and it would
make the protein move instead).

The per-frame RMSD returned by the aligner is a diagnostic that must be inspected:
a poor fit inflates the apparent positional spread of the waters, which downstream
looks exactly like disorder and produces false "easy to displace" hits.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.analysis.align import rotation_matrix

from hydrarank.exceptions import EmptySelectionError


@dataclass(frozen=True)
class AlignmentReference:
    """Reference coordinates of the fitting group, centred on the origin."""

    indices: np.ndarray
    positions: np.ndarray
    centroid: np.ndarray
    frame: int

    @property
    def n_atoms(self) -> int:
        return int(self.positions.shape[0])


def select_alignment_group(universe: Universe, ligand: AtomGroup, cutoff: float = 8.0) -> AtomGroup:
    """Backbone alpha carbons of the protein residues lining the ligand."""
    pocket = f"protein and byres (around {cutoff} group ligand)"
    group = universe.select_atoms(f"name CA and {pocket}", ligand=ligand)
    if group.n_atoms == 0:
        group = universe.select_atoms(f"backbone and {pocket}", ligand=ligand)
    if group.n_atoms < 3:
        raise EmptySelectionError(
            f"only {group.n_atoms} backbone atoms found within {cutoff} A of the ligand; "
            "a rigid-body fit needs at least 3. Increase the pocket cutoff or check that "
            "the ligand is not a periodic image away from the protein."
        )
    return group


def build_reference(group: AtomGroup, frame: int = 0) -> AlignmentReference:
    """Snapshot ``group`` at ``frame`` as the reference, restoring the frame afterwards."""
    universe = group.universe
    original = universe.trajectory.frame
    try:
        universe.trajectory[frame]
        positions = group.positions.astype(np.float64, copy=True)
    finally:
        universe.trajectory[original]

    centroid = positions.mean(axis=0)
    return AlignmentReference(
        indices=group.indices.copy(),
        positions=positions - centroid,
        centroid=centroid,
        frame=frame,
    )


@dataclass(frozen=True)
class RigidMotion:
    """The rotation/translation that maps one frame onto the reference frame."""

    rotation: np.ndarray
    source_centroid: np.ndarray
    target_centroid: np.ndarray
    rmsd: float

    def apply(self, positions: np.ndarray) -> np.ndarray:
        moved = (np.asarray(positions, dtype=np.float64) - self.source_centroid) @ self.rotation.T
        return moved + self.target_centroid


class SiteAligner:
    """Rigid-body fit of arbitrary coordinates onto an :class:`AlignmentReference`."""

    def __init__(self, reference: AlignmentReference):
        self.reference = reference

    def fit(self, fit_positions: np.ndarray) -> RigidMotion:
        """Compute the motion that superposes ``fit_positions`` onto the reference."""
        fit_positions = np.asarray(fit_positions, dtype=np.float64)
        if fit_positions.shape[0] != self.reference.n_atoms:
            raise ValueError(
                f"fitting group has {fit_positions.shape[0]} atoms but the reference has "
                f"{self.reference.n_atoms}"
            )
        centroid = fit_positions.mean(axis=0)
        rotation, rmsd = rotation_matrix(fit_positions - centroid, self.reference.positions)
        return RigidMotion(
            rotation=rotation,
            source_centroid=centroid,
            target_centroid=self.reference.centroid,
            rmsd=float(rmsd),
        )

    def transform(
        self, fit_positions: np.ndarray, positions: np.ndarray
    ) -> tuple[np.ndarray, float]:
        """Fit ``fit_positions`` onto the reference and apply the same motion to ``positions``.

        Returns the transformed coordinates and the fit RMSD in angstrom.
        """
        motion = self.fit(fit_positions)
        return motion.apply(positions), motion.rmsd


def iter_alignment_rmsd(
    universe: Universe,
    group: AtomGroup,
    aligner: SiteAligner,
    frames: slice | None = None,
) -> Iterator[tuple[int, float]]:
    """Yield ``(frame_index, fit_rmsd)`` for a quick convergence check."""
    trajectory = universe.trajectory[frames] if frames is not None else universe.trajectory
    for ts in trajectory:
        _, rmsd = aligner.transform(group.positions, group.positions)
        yield ts.frame, rmsd
