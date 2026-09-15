"""Diagnostics for the preprocessing stage.

These numbers decide whether the downstream hydration-site analysis can be trusted
at all, so they are computed and reported explicitly rather than assumed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.lib.distances import capped_distance

from hydrarank.config import PreprocessConfig
from hydrarank.data import WaterObservations
from hydrarank.exceptions import HydraRankError
from hydrarank.pbc import MAX_BOND_LENGTH, max_bond_length
from hydrarank.preprocess import prepare_system


@dataclass
class PreprocessQC:
    """Per-frame diagnostics of the PBC + alignment stack."""

    frames: np.ndarray
    fit_rmsd: np.ndarray
    longest_bond: np.ndarray
    n_shell_waters: np.ndarray
    n_fit_atoms: int
    reference_frame: int
    water_cutoff: float

    @property
    def solute_always_whole(self) -> bool:
        return bool(np.max(self.longest_bond) <= MAX_BOND_LENGTH)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_frames": int(self.frames.size),
            "n_fit_atoms": self.n_fit_atoms,
            "reference_frame": self.reference_frame,
            "water_cutoff": self.water_cutoff,
            "fit_rmsd": {
                "mean": float(self.fit_rmsd.mean()),
                "max": float(self.fit_rmsd.max()),
                "last": float(self.fit_rmsd[-1]),
            },
            "longest_bond_max": float(self.longest_bond.max()),
            "solute_always_whole": self.solute_always_whole,
            "n_shell_waters": {
                "mean": float(self.n_shell_waters.mean()),
                "min": int(self.n_shell_waters.min()),
                "max": int(self.n_shell_waters.max()),
            },
        }

    def __str__(self) -> str:
        lines = [
            "Preprocessing QC",
            f"  frames analysed    : {self.frames.size}",
            f"  fit atoms          : {self.n_fit_atoms} (reference frame {self.reference_frame})",
            f"  fit RMSD           : mean {self.fit_rmsd.mean():.2f} A, "
            f"max {self.fit_rmsd.max():.2f} A, last {self.fit_rmsd[-1]:.2f} A",
            f"  longest bond       : {self.longest_bond.max():.2f} A "
            f"({'solute whole' if self.solute_always_whole else 'SOLUTE SPLIT'})",
            f"  waters within {self.water_cutoff:g} A: "
            f"mean {self.n_shell_waters.mean():.1f} "
            f"(min {self.n_shell_waters.min()}, max {self.n_shell_waters.max()})",
        ]
        if not self.solute_always_whole:
            lines.append("  ! unwrapping failed; check the bond information in the topology")
        if self.fit_rmsd.max() > 2.5:
            lines.append(
                "  ! large fit RMSD: the pocket moves relative to the reference, which blurs "
                "hydration sites and inflates apparent water disorder"
            )
        return "\n".join(lines)


def run_preprocess_qc(universe: Universe, config: PreprocessConfig) -> PreprocessQC:
    """Apply the PBC and alignment stack and collect per-frame diagnostics."""
    system = prepare_system(universe, config)

    indices, rmsds, bonds, counts = [], [], [], []
    for ts in system.universe.trajectory[system.frames]:
        motion = system.aligner.fit(system.fit_group.positions)
        pairs = capped_distance(
            motion.apply(system.oxygens.positions),
            system.ligand_reference,
            max_cutoff=config.water_cutoff,
            return_distances=False,
        )
        indices.append(ts.frame)
        rmsds.append(motion.rmsd)
        bonds.append(max_bond_length(system.groups.solute))
        counts.append(np.unique(pairs[:, 0]).size if len(pairs) else 0)

    return PreprocessQC(
        frames=np.asarray(indices),
        fit_rmsd=np.asarray(rmsds),
        longest_bond=np.asarray(bonds),
        n_shell_waters=np.asarray(counts),
        n_fit_atoms=system.fit_group.n_atoms,
        reference_frame=system.reference_frame,
        water_cutoff=config.water_cutoff,
    )


def qc_from_observations(observations: WaterObservations) -> PreprocessQC:
    """Reconstruct QC collected during preprocessing without rereading the trajectory."""
    metadata = observations.metadata
    try:
        longest_bond = np.asarray(metadata["qc_longest_bond"], dtype=float)
        n_shell_waters = np.asarray(metadata["qc_n_shell_waters"], dtype=int)
        n_fit_atoms = int(metadata["n_fit_atoms"])
        reference_frame = int(metadata["reference_frame"])
        water_cutoff = float(metadata["water_cutoff"])
    except (KeyError, TypeError, ValueError) as error:
        raise HydraRankError(
            "the observations cache does not contain preprocessing QC; regenerate it with "
            "'hydrarank analyse --force'"
        ) from error
    if longest_bond.shape != observations.frames.shape:
        raise HydraRankError("cached longest-bond QC does not match the selected frame count")
    if n_shell_waters.shape != observations.frames.shape:
        raise HydraRankError("cached water-shell QC does not match the selected frame count")
    return PreprocessQC(
        frames=observations.frames,
        fit_rmsd=observations.fit_rmsd,
        longest_bond=longest_bond,
        n_shell_waters=n_shell_waters,
        n_fit_atoms=n_fit_atoms,
        reference_frame=reference_frame,
        water_cutoff=water_cutoff,
    )
