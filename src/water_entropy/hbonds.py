"""Geometric hydrogen bonds between the shell waters and the solute.

A water can only be displaced cheaply if the ligand does not have to replace strong
interactions. Counting hydrogen bonds to protein and ligand, together with how
enclosed the water is, is the cheapest useful proxy for that.

The criterion is the usual one: heavy-atom donor-acceptor distance below 3.5 A and a
donor-H...acceptor angle above 130 degrees. It needs explicit hydrogens in the
topology, which PSF, TPR and prmtop all provide.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from MDAnalysis import AtomGroup
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.lib.distances import capped_distance

from water_entropy.exceptions import WaterEntropyError

DEFAULT_HBOND_DISTANCE = 3.5
DEFAULT_HBOND_ANGLE = 130.0
DEFAULT_ENCLOSURE_RADIUS = 5.0

_HYDROGEN_MAX_MASS = 4.5
_POLAR_MASS_RANGES = ((13.5, 14.6), (15.0, 17.0), (31.5, 32.5))  # N, O, S


@dataclass(frozen=True)
class PolarGroups:
    """Hydrogen-bonding partners of one solute selection, plus its heavy atoms."""

    acceptors: AtomGroup
    donors: AtomGroup
    """Heavy donor atoms, one entry per donor hydrogen."""

    donor_hydrogens: AtomGroup
    heavy: AtomGroup

    @property
    def n_donors(self) -> int:
        return int(self.donor_hydrogens.n_atoms)


def find_polar_groups(atoms: AtomGroup) -> PolarGroups:
    """Locate acceptors, donor/hydrogen pairs and heavy atoms in a solute selection."""
    universe = atoms.universe
    try:
        masses = universe.atoms.masses
    except (NoDataError, AttributeError):
        raise WaterEntropyError(
            "hydrogen-bond analysis requires atom masses to identify hydrogens and polar atoms; "
            "use a topology with masses (such as PSF, TPR or prmtop)"
        ) from None
    if masses.size != universe.atoms.n_atoms or not np.any(masses > 0):
        raise WaterEntropyError(
            "hydrogen-bond analysis requires valid atom masses; the topology contains none"
        )
    polar = _is_polar(masses)
    hydrogen = masses < _HYDROGEN_MAX_MASS

    selected = np.zeros(universe.atoms.n_atoms, dtype=bool)
    selected[atoms.indices] = True

    donor_heavy, donor_h = [], []
    if hasattr(universe.atoms, "bonds"):
        bonds = atoms.bonds.to_indices()
        for left, right in ((bonds[:, 0], bonds[:, 1]), (bonds[:, 1], bonds[:, 0])):
            mask = polar[left] & hydrogen[right] & selected[left] & selected[right]
            donor_heavy.append(left[mask])
            donor_h.append(right[mask])

    return PolarGroups(
        acceptors=universe.atoms[np.flatnonzero(selected & polar)],
        donors=universe.atoms[np.concatenate(donor_heavy) if donor_heavy else np.empty(0, int)],
        donor_hydrogens=universe.atoms[np.concatenate(donor_h) if donor_h else np.empty(0, int)],
        heavy=universe.atoms[np.flatnonzero(selected & ~hydrogen & (masses > _HYDROGEN_MAX_MASS))],
    )


def count_hbonds(
    water_oxygen: np.ndarray,
    water_hydrogen: np.ndarray,
    polar: PolarGroups,
    distance: float = DEFAULT_HBOND_DISTANCE,
    angle: float = DEFAULT_HBOND_ANGLE,
) -> np.ndarray:
    """Hydrogen bonds between each water and the solute, counting both directions."""
    counts = np.zeros(water_oxygen.shape[0])
    if water_oxygen.shape[0] == 0:
        return counts

    acceptors = polar.acceptors.positions
    if acceptors.shape[0]:
        pairs = capped_distance(
            water_oxygen, acceptors, max_cutoff=distance, return_distances=False
        )
        for site in range(2):
            donated = _angle_ok(
                water_oxygen[pairs[:, 0]],
                water_hydrogen[pairs[:, 0], site],
                acceptors[pairs[:, 1]],
                angle,
            )
            np.add.at(counts, pairs[donated, 0], 1.0)

    if polar.n_donors:
        donors = polar.donors.positions
        hydrogens = polar.donor_hydrogens.positions
        pairs = capped_distance(donors, water_oxygen, max_cutoff=distance, return_distances=False)
        accepted = _angle_ok(
            donors[pairs[:, 0]],
            hydrogens[pairs[:, 0]],
            water_oxygen[pairs[:, 1]],
            angle,
        )
        np.add.at(counts, pairs[accepted, 1], 1.0)

    return counts


def count_neighbours(
    water_oxygen: np.ndarray,
    heavy: AtomGroup,
    radius: float = DEFAULT_ENCLOSURE_RADIUS,
) -> np.ndarray:
    """Solute heavy atoms around each water: a cheap measure of how buried it is."""
    counts = np.zeros(water_oxygen.shape[0])
    if water_oxygen.shape[0] == 0 or heavy.n_atoms == 0:
        return counts
    pairs = capped_distance(
        water_oxygen, heavy.positions, max_cutoff=radius, return_distances=False
    )
    if len(pairs):
        np.add.at(counts, pairs[:, 0], 1.0)
    return counts


def _angle_ok(
    donor: np.ndarray, hydrogen: np.ndarray, acceptor: np.ndarray, threshold: float
) -> np.ndarray:
    """True where the donor-H...acceptor angle exceeds ``threshold`` degrees."""
    if donor.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    to_donor = donor - hydrogen
    to_acceptor = acceptor - hydrogen
    norms = np.linalg.norm(to_donor, axis=1) * np.linalg.norm(to_acceptor, axis=1)
    cosine = np.einsum("ij,ij->i", to_donor, to_acceptor) / np.where(norms > 0, norms, 1.0)
    return cosine < np.cos(np.radians(threshold))


def _is_polar(masses: np.ndarray) -> np.ndarray:
    polar = np.zeros(masses.shape, dtype=bool)
    for low, high in _POLAR_MASS_RANGES:
        polar |= (masses >= low) & (masses <= high)
    return polar
