import numpy as np
import pytest
from MDAnalysis import Universe

from hydrarank.exceptions import HydraRankError
from hydrarank.hbonds import count_hbonds, count_neighbours, find_polar_groups


def test_water_donates_hydrogen_bond_to_solute_acceptor(build_universe):
    universe = build_universe(
        [("SOL", [("O", 15.999)])],
        positions=np.array([[2.8, 0.0, 0.0]]),
    )
    polar = find_polar_groups(universe.atoms)
    water_oxygen = np.array([[0.0, 0.0, 0.0]])
    water_hydrogen = np.array([[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])

    assert count_hbonds(water_oxygen, water_hydrogen, polar)[0] == 1


def test_solute_donates_hydrogen_bond_to_water(build_universe):
    universe = build_universe(
        [("SOL", [("N", 14.007), ("H", 1.008)])],
        positions=np.array([[2.8, 0.0, 0.0], [1.8, 0.0, 0.0]]),
        bonds=[(0, 1)],
    )
    polar = find_polar_groups(universe.atoms)
    water_oxygen = np.array([[0.0, 0.0, 0.0]])
    water_hydrogen = np.array([[[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]]])

    assert count_hbonds(water_oxygen, water_hydrogen, polar)[0] == 1


def test_enclosure_counts_nearby_solute_heavy_atoms(build_universe):
    universe = build_universe(
        [("SOL", [("C1", 12.011), ("C2", 12.011), ("H", 1.008)])],
        positions=np.array([[1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
    )
    polar = find_polar_groups(universe.atoms)

    assert count_neighbours(np.array([[0.0, 0.0, 0.0]]), polar.heavy, radius=2.0)[0] == 1


def test_missing_masses_raise_an_actionable_error():
    universe = Universe.empty(1, trajectory=True)
    with pytest.raises(HydraRankError, match="requires atom masses"):
        find_polar_groups(universe.atoms)
