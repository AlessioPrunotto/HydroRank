import numpy as np
import pytest

from water_entropy.exceptions import TrajectoryError, WaterEntropyError
from water_entropy.pbc import (
    apply_pbc_transformations,
    assert_solute_whole,
    build_pbc_groups,
    max_bond_length,
    validate_cutoff,
)

BOX = (50.0, 50.0, 50.0, 90.0, 90.0, 90.0)
LIGAND = ("547", [(f"C{i}", 12.011) for i in range(4)])
WATER = ("OPC", [("OW", 15.999), ("HW1", 1.008), ("HW2", 1.008), ("MW", 0.0)])

# ligand straddling the x boundary, water sitting next to its periodic image
LIGAND_XYZ = [
    [1.0, 25.0, 25.0],
    [0.5, 25.0, 25.0],
    [49.5, 25.0, 25.0],
    [49.0, 25.0, 25.0],
]
WATER_XYZ = [
    [47.0, 25.0, 25.0],
    [47.9, 25.0, 25.0],
    [46.7, 25.9, 25.0],
    [47.1, 25.1, 25.0],
]
LIGAND_BONDS = [(0, 1), (1, 2), (2, 3)]


@pytest.fixture
def straddling_system(build_universe):
    return build_universe(
        [LIGAND, WATER],
        n_frames=2,
        box=BOX,
        positions=np.array(LIGAND_XYZ + WATER_XYZ),
        bonds=LIGAND_BONDS,
    )


def test_solute_is_split_before_transformation(straddling_system):
    ligand = straddling_system.select_atoms("resname 547")
    assert max_bond_length(ligand) > 3.0
    with pytest.raises(WaterEntropyError, match="not whole"):
        assert_solute_whole(ligand)


def test_transformations_make_solute_whole_and_centre_it(straddling_system):
    ligand = straddling_system.select_atoms("resname 547")
    groups = build_pbc_groups(straddling_system, ligand)
    assert groups.mobile.n_atoms == 4  # the water only
    apply_pbc_transformations(straddling_system, groups)

    for _ in straddling_system.trajectory:
        assert_solute_whole(ligand)
        assert np.allclose(ligand.center_of_geometry(), [25.0, 25.0, 25.0], atol=1e-3)


def test_water_is_wrapped_next_to_the_ligand_without_being_split(straddling_system):
    ligand = straddling_system.select_atoms("resname 547")
    water = straddling_system.select_atoms("resname OPC")
    apply_pbc_transformations(straddling_system, build_pbc_groups(straddling_system, ligand))

    straddling_system.trajectory[0]
    oxygen = water.positions[0]
    spread = np.linalg.norm(water.positions - oxygen, axis=1)
    assert spread.max() < 2.0  # molecule intact
    assert np.linalg.norm(oxygen - ligand.center_of_geometry()) < 5.0


def test_transformations_applied_only_once(straddling_system):
    ligand = straddling_system.select_atoms("resname 547")
    groups = build_pbc_groups(straddling_system, ligand)
    apply_pbc_transformations(straddling_system, groups)
    with pytest.raises(WaterEntropyError, match="already been applied"):
        apply_pbc_transformations(straddling_system, groups)


def test_missing_bonds_raise_with_actionable_message(build_universe):
    universe = build_universe([LIGAND, WATER], box=BOX)
    ligand = universe.select_atoms("resname 547")
    with pytest.raises(TrajectoryError, match="no bonds"):
        apply_pbc_transformations(universe, build_pbc_groups(universe, ligand))


def test_guessed_bonds_allow_unwrapping(build_universe):
    universe = build_universe([LIGAND, WATER], box=BOX, positions=np.array(LIGAND_XYZ + WATER_XYZ))
    universe.add_TopologyAttr("type", ["C"] * 4 + ["O", "H", "H", "M"])
    ligand = universe.select_atoms("resname 547")
    apply_pbc_transformations(universe, build_pbc_groups(universe, ligand), guess_bonds=True)
    assert universe.trajectory.transformations


def test_validate_cutoff(straddling_system):
    validate_cutoff(straddling_system, 5.0)
    with pytest.raises(WaterEntropyError, match="pocket_cutoff.*minimum-image"):
        validate_cutoff(straddling_system, 25.0, name="pocket_cutoff")
