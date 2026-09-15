"""Synthetic in-memory systems, so the tests need no trajectory files."""

from __future__ import annotations

import numpy as np
import pytest
from MDAnalysis import Universe


def _build(
    residue_specs,
    n_frames: int = 3,
    box=(50.0, 50.0, 50.0, 90.0, 90.0, 90.0),
    positions=None,
    bonds=None,
):
    """Build a Universe from ``[(resname, [(name, mass), ...]), ...]``."""
    names, masses, resnames, atom_resindex = [], [], [], []
    for resindex, (resname, atoms) in enumerate(residue_specs):
        resnames.append(resname)
        for name, mass in atoms:
            names.append(name)
            masses.append(mass)
            atom_resindex.append(resindex)

    universe = Universe.empty(
        n_atoms=len(names),
        n_residues=len(resnames),
        atom_resindex=np.array(atom_resindex),
        trajectory=True,
    )
    universe.add_TopologyAttr("name", names)
    universe.add_TopologyAttr("mass", masses)
    universe.add_TopologyAttr("resname", resnames)
    universe.add_TopologyAttr("resid", list(range(1, len(resnames) + 1)))
    if bonds is not None:
        universe.add_bonds(bonds)

    if positions is None:
        rng = np.random.default_rng(0)
        coordinates = rng.uniform(0, 20, size=(n_frames, len(names), 3))
    else:
        coordinates = np.asarray(positions, dtype=float)
        if coordinates.ndim == 2:
            coordinates = np.repeat(coordinates[None], n_frames, axis=0)
    universe.load_new(coordinates.astype(np.float32), order="fac")
    for ts in universe.trajectory:
        ts.dimensions = np.array(box, dtype=np.float32)
    universe.trajectory[0]
    return universe


OPC_WATER = ("OPC", [("OW", 15.999), ("HW1", 1.008), ("HW2", 1.008), ("MW", 0.0)])
TIP3_WATER = ("TIP3", [("OH2", 15.999), ("H1", 1.008), ("H2", 1.008)])
LIGAND = ("547", [(f"C{i}", 12.011) for i in range(10)] + [("H1", 1.008)])
SODIUM = ("SOD", [("SOD", 22.99)])


@pytest.fixture
def opc_system() -> Universe:
    return _build([LIGAND, SODIUM, *[OPC_WATER] * 5])


@pytest.fixture
def tip3_system() -> Universe:
    return _build([LIGAND, *[TIP3_WATER] * 4])


@pytest.fixture
def build_universe():
    return _build


ALA = [("N", 14.007), ("CA", 12.011), ("C", 12.011), ("O", 15.999), ("CB", 12.011)]
POCKET_LIGAND = (
    "547",
    [("C0", 12.011), ("N1", 14.007), ("S1", 32.06), ("O1", 15.999)],
)
N_POCKET_RESIDUES = 4


def _pocket_positions():
    protein = [
        [20.0 + 2.0 * r + 0.4 * i, 25.0, 25.0] for r in range(N_POCKET_RESIDUES) for i in range(5)
    ]
    ligand = [[24.0, 27.0, 25.0], [25.0, 27.0, 25.0], [24.0, 28.0, 25.0], [25.0, 28.0, 25.0]]
    waters = []
    for offset in (0.0, 2.0):
        origin = np.array([24.5, 29.5 + offset, 25.0])
        waters += [
            (origin + shift).tolist()
            for shift in ([0, 0, 0], [0.9, 0, 0], [-0.3, 0.85, 0], [0.1, 0.1, 0])
        ]
    return np.array(protein + ligand + waters)


def _pocket_bonds():
    bonds = []
    for r in range(N_POCKET_RESIDUES):
        base = 5 * r
        bonds += [
            (base, base + 1),
            (base + 1, base + 2),
            (base + 2, base + 3),
            (base + 1, base + 4),
        ]
    base = 5 * N_POCKET_RESIDUES
    bonds += [(base, base + 1), (base + 1, base + 2), (base + 2, base + 3)]
    return bonds


@pytest.fixture
def pocket_system() -> Universe:
    """Four alanines, a bonded ligand and two waters, all inside the box."""
    residues = [("ALA", ALA) for _ in range(N_POCKET_RESIDUES)] + [
        POCKET_LIGAND,
        OPC_WATER,
        OPC_WATER,
    ]
    return _build(residues, n_frames=4, positions=_pocket_positions(), bonds=_pocket_bonds())
