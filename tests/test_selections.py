import numpy as np
import pytest

from water_entropy.exceptions import (
    AmbiguousSelectionError,
    EmptySelectionError,
    WaterModelError,
)
from water_entropy.selections import (
    analyse_water_topology,
    classify_atoms,
    select_ligand,
    select_water,
)


def test_opc_water_is_found_and_virtual_site_separated(opc_system):
    water = select_water(opc_system)
    assert water.n_atoms == 20

    topology = analyse_water_topology(water)
    assert topology.n_waters == 5
    assert topology.n_sites == 4
    assert topology.virtual_ix.shape == (5, 1)
    assert topology.hydrogen_ix.shape == (5, 2)
    assert "4-site" in topology.model

    # the virtual site must never leak into the real-atom set
    assert topology.real_atom_ix().size == 15
    assert not np.intersect1d(topology.real_atom_ix(), topology.virtual_ix.ravel()).size


def test_three_site_water(tip3_system):
    topology = analyse_water_topology(select_water(tip3_system))
    assert topology.n_sites == 3
    assert topology.virtual_ix.shape == (4, 0)
    assert topology.real_atom_ix().size == 12


def test_unknown_water_resname_raises_with_helpful_message(opc_system):
    with pytest.raises(EmptySelectionError, match="OPC"):
        select_water(opc_system, resnames=("HOH", "TIP3"))


def test_mixed_water_sizes_rejected(build_universe):
    universe = build_universe(
        [
            ("OPC", [("OW", 15.999), ("HW1", 1.008), ("HW2", 1.008), ("MW", 0.0)]),
            ("OPC", [("OW", 15.999), ("HW1", 1.008), ("HW2", 1.008)]),
        ]
    )
    with pytest.raises(WaterModelError, match="inconsistent sizes"):
        analyse_water_topology(select_water(universe))


def test_ligand_autodetected_ignoring_water_and_ions(opc_system):
    ligand = select_ligand(opc_system)
    assert set(ligand.residues.resnames) == {"547"}
    assert ligand.n_atoms == 11


def test_ligand_autodetection_refuses_when_ambiguous(build_universe):
    ligand = ("547", [(f"C{i}", 12.011) for i in range(10)])
    universe = build_universe([ligand, ("ABC", [(f"C{i}", 12.011) for i in range(10)])])
    with pytest.raises(AmbiguousSelectionError, match="2 candidates"):
        select_ligand(universe)


def test_explicit_ligand_selection_wins(opc_system):
    assert select_ligand(opc_system, "resname 547").n_atoms == 11


def test_empty_explicit_ligand_selection_raises(opc_system):
    with pytest.raises(EmptySelectionError):
        select_ligand(opc_system, "resname NOPE")


def test_classify_atoms_uses_mass_not_name(build_universe):
    universe = build_universe([("OPC", [("O00", 15.999), ("H", 1.008), ("H", 1.008), ("X", 0.0)])])
    assert list(classify_atoms(universe.atoms)) == ["O", "H", "H", "M"]
