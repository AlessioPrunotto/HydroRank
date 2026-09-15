"""Robust selection of ligand, water and pocket atoms.

The selections here deliberately avoid the built-in ``water`` / ``protein``
shortcuts where those are known to be unreliable across force fields, and they
identify atom types by mass rather than by name, because water oxygen names vary
(``OW``, ``OH2``, ``O``, ...) and 4-/5-site models carry massless virtual sites.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.exceptions import NoDataError

from hydrarank.config import DEFAULT_ION_RESNAMES, DEFAULT_WATER_RESNAMES
from hydrarank.exceptions import (
    AmbiguousSelectionError,
    EmptySelectionError,
    WaterModelError,
)

_VIRTUAL_SITE_MAX_MASS = 0.1
_HYDROGEN_MAX_MASS = 4.5  # generous, to tolerate hydrogen mass repartitioning
_OXYGEN_MASS_RANGE = (15.0, 17.0)

_WATER_MODEL_BY_N_SITES = {3: "3-site (SPC/TIP3P-like)", 4: "4-site (TIP4P/OPC-like)", 5: "5-site"}


@dataclass(frozen=True)
class WaterTopology:
    """Per-residue atom indices of the water molecules, split by role.

    All index arrays refer to absolute atom indices in the parent ``Universe``.
    """

    resids: np.ndarray  # (n_waters,)
    resindices: np.ndarray  # (n_waters,), globally unique within the topology
    oxygen_ix: np.ndarray  # (n_waters,)
    hydrogen_ix: np.ndarray  # (n_waters, 2)
    virtual_ix: np.ndarray  # (n_waters, n_virtual), possibly with 0 columns
    n_sites: int

    @property
    def n_waters(self) -> int:
        return int(self.resids.size)

    @property
    def model(self) -> str:
        return _WATER_MODEL_BY_N_SITES.get(self.n_sites, f"{self.n_sites}-site (unknown)")

    def real_atom_ix(self) -> np.ndarray:
        """Indices of all non-virtual water atoms, sorted."""
        stacked = np.concatenate([self.oxygen_ix[:, None], self.hydrogen_ix], axis=1)
        return np.sort(stacked.ravel())


def select_water(
    universe: Universe, resnames: tuple[str, ...] = DEFAULT_WATER_RESNAMES
) -> AtomGroup:
    """Select all water atoms (including virtual sites) by residue name."""
    water = universe.select_atoms("resname " + " ".join(resnames))
    if water.n_atoms == 0:
        present = sorted(set(universe.residues.resnames))
        raise EmptySelectionError(
            "no water found. Residue names searched: "
            f"{list(resnames)}. Residue names present in the system: {present}. "
            "Pass the correct name via 'water_resnames' in the configuration."
        )
    return water


def analyse_water_topology(water: AtomGroup) -> WaterTopology:
    """Split a water selection into oxygen / hydrogen / virtual-site indices.

    Raises if the water residues are not all of the same geometry, which usually
    means the water resname list has caught something that is not water.
    """
    if water.n_atoms == 0:
        raise EmptySelectionError("empty water selection")

    resindices = water.resindices
    if np.any(np.diff(resindices) < 0):
        raise WaterModelError("water atoms are not grouped by residue in the topology")

    _, counts = np.unique(resindices, return_counts=True)
    n_sites = int(counts[0])
    if not np.all(counts == n_sites):
        sizes = sorted(set(counts.tolist()))
        raise WaterModelError(
            f"water residues have inconsistent sizes {sizes}; the water resname list "
            "probably matches non-water residues"
        )
    if n_sites not in _WATER_MODEL_BY_N_SITES:
        raise WaterModelError(f"unsupported water model with {n_sites} sites per molecule")

    roles = classify_atoms(water).reshape(-1, n_sites)
    if not np.all(roles == roles[0]):
        raise WaterModelError("water residues do not all have the same atom ordering")

    pattern = roles[0]
    indices = water.indices.reshape(-1, n_sites)
    oxygen_cols = np.flatnonzero(pattern == "O")
    hydrogen_cols = np.flatnonzero(pattern == "H")
    virtual_cols = np.flatnonzero(pattern == "M")

    if oxygen_cols.size != 1 or hydrogen_cols.size != 2:
        raise WaterModelError(
            f"expected 1 oxygen and 2 hydrogens per water, found {oxygen_cols.size} and "
            f"{hydrogen_cols.size} (atom names: {list(water[:n_sites].names)})"
        )

    return WaterTopology(
        resids=water.residues.resids.copy(),
        resindices=water.residues.resindices.copy(),
        oxygen_ix=indices[:, oxygen_cols[0]],
        hydrogen_ix=indices[:, hydrogen_cols],
        virtual_ix=indices[:, virtual_cols],
        n_sites=n_sites,
    )


def select_ligand(
    universe: Universe,
    selection: str | None = None,
    water_resnames: tuple[str, ...] = DEFAULT_WATER_RESNAMES,
    ion_resnames: tuple[str, ...] = DEFAULT_ION_RESNAMES,
    min_heavy_atoms: int = 8,
) -> AtomGroup:
    """Select the ligand, either from an explicit selection or by auto-detection.

    Auto-detection is a convenience for interactive use only: it refuses to guess
    when several candidate residues are found.
    """
    if selection is not None:
        ligand = universe.select_atoms(selection)
        if ligand.n_atoms == 0:
            raise EmptySelectionError(f"ligand selection {selection!r} matched no atoms")
        return ligand

    candidates = _ligand_candidates(universe, water_resnames, ion_resnames, min_heavy_atoms)
    if not candidates:
        raise EmptySelectionError(
            "could not auto-detect a ligand: no non-protein residue with at least "
            f"{min_heavy_atoms} heavy atoms. Set 'ligand_selection' explicitly."
        )
    if len(candidates) > 1:
        labels = [f"resname {r.resname} resid {r.resid}" for r in candidates]
        raise AmbiguousSelectionError(
            f"ligand auto-detection found {len(candidates)} candidates: {labels}. "
            "Set 'ligand_selection' explicitly."
        )
    return candidates[0].atoms


def select_pocket(universe: Universe, ligand: AtomGroup, cutoff: float = 8.0) -> AtomGroup:
    """Whole protein residues with any atom within ``cutoff`` of the ligand.

    Evaluated on the current frame only, and without periodic images: call this
    after the PBC/centering transformations are in place.
    """
    pocket = universe.select_atoms(
        f"protein and byres (around {cutoff} group ligand)", ligand=ligand
    )
    if pocket.n_atoms == 0:
        raise EmptySelectionError(
            f"no protein residue within {cutoff} A of the ligand on the current frame; "
            "the ligand may be a periodic image away from the protein"
        )
    return pocket


def _ligand_candidates(
    universe: Universe,
    water_resnames: tuple[str, ...],
    ion_resnames: tuple[str, ...],
    min_heavy_atoms: int,
) -> list:
    excluded = " ".join({name.upper() for name in (*water_resnames, *ion_resnames)})
    others = universe.select_atoms(f"not protein and not nucleic and not resname {excluded}")
    if others.n_atoms == 0:
        return []

    roles = classify_atoms(others)
    is_heavy = (roles != "H") & (roles != "M")
    residues = others.residues
    resindex_of_atom = np.searchsorted(residues.resindices, others.resindices)
    heavy_per_residue = np.bincount(
        resindex_of_atom, weights=is_heavy, minlength=residues.n_residues
    )
    return [residues[i] for i in np.flatnonzero(heavy_per_residue >= min_heavy_atoms)]


def classify_atoms(atoms: AtomGroup) -> np.ndarray:
    """Return an array of ``'O'``/``'H'``/``'M'``/``'X'`` roles, one per atom.

    ``'M'`` marks massless virtual sites, ``'X'`` any other real element.
    """
    try:
        masses = np.asarray(atoms.masses, dtype=float)
    except (NoDataError, AttributeError):
        masses = None

    if masses is None or masses.size == 0 or not np.any(masses > 0):
        return _classify_by_name(atoms)

    roles = np.full(atoms.n_atoms, "X", dtype="<U1")
    roles[masses <= _VIRTUAL_SITE_MAX_MASS] = "M"
    roles[(masses > _VIRTUAL_SITE_MAX_MASS) & (masses < _HYDROGEN_MAX_MASS)] = "H"
    roles[(masses >= _OXYGEN_MASS_RANGE[0]) & (masses <= _OXYGEN_MASS_RANGE[1])] = "O"
    return roles


def _classify_by_name(atoms: AtomGroup) -> np.ndarray:
    names = np.asarray([str(name).upper().lstrip("0123456789") for name in atoms.names])
    roles = np.full(names.size, "X", dtype="<U1")
    for i, name in enumerate(names):
        if not name:
            continue
        if name.startswith(("MW", "EP", "LP", "DW")):
            roles[i] = "M"
        elif name.startswith("H"):
            roles[i] = "H"
        elif name.startswith("O"):
            roles[i] = "O"
    return roles
