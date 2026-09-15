import numpy as np
import pytest
from MDAnalysis.lib.transformations import rotation_matrix as rotation_about_axis

from water_entropy.alignment import (
    AlignmentReference,
    SiteAligner,
    build_reference,
    iter_alignment_rmsd,
    select_alignment_group,
)
from water_entropy.exceptions import EmptySelectionError

ALA = [("N", 14.007), ("CA", 12.011), ("C", 12.011), ("O", 15.999), ("CB", 12.011)]
LIGAND = ("547", [(f"C{i}", 12.011) for i in range(4)])


def _protein_ligand(build_universe, n_frames=3):
    residues = [("ALA", ALA) for _ in range(4)] + [LIGAND]
    x_offsets = [0.0, 5.0, 12.0, 30.0]
    positions = [[x + i, 0.0, 0.0] for x in x_offsets for i in range(5)]
    positions += [[2.0, 1.0, 0.0], [3.0, 1.0, 0.0], [2.0, 2.0, 0.0], [3.0, 2.0, 0.0]]
    return build_universe(residues, n_frames=n_frames, positions=np.array(positions))


def test_alignment_group_is_the_pocket_backbone(build_universe):
    universe = _protein_ligand(build_universe)
    ligand = universe.select_atoms("resname 547")
    group = select_alignment_group(universe, ligand, cutoff=20.0)
    assert group.n_atoms == 4
    assert set(group.names) == {"CA"}


def test_alignment_group_too_small_raises(build_universe):
    universe = _protein_ligand(build_universe)
    ligand = universe.select_atoms("resname 547")
    with pytest.raises(EmptySelectionError, match="at least 3"):
        select_alignment_group(universe, ligand, cutoff=8.0)


def test_aligner_recovers_a_rigid_body_motion():
    rng = np.random.default_rng(1)
    reference_positions = rng.uniform(-10, 10, size=(12, 3))
    centroid = reference_positions.mean(axis=0)
    reference = AlignmentReference(
        indices=np.arange(12),
        positions=reference_positions - centroid,
        centroid=centroid,
        frame=0,
    )

    rotation = rotation_about_axis(0.7, [0.3, -0.5, 0.8])[:3, :3]
    moved = reference_positions @ rotation.T + np.array([13.0, -4.0, 7.5])

    aligned, rmsd = SiteAligner(reference).transform(moved, moved)
    assert rmsd < 1e-6
    assert np.allclose(aligned, reference_positions, atol=1e-6)


def test_aligner_applies_the_fit_to_other_coordinates():
    reference_positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    centroid = reference_positions.mean(axis=0)
    reference = AlignmentReference(
        indices=np.arange(3),
        positions=reference_positions - centroid,
        centroid=centroid,
        frame=0,
    )
    shift = np.array([100.0, 0.0, 0.0])
    water = np.array([[0.5, 0.5, 3.0]]) + shift

    aligned, rmsd = SiteAligner(reference).transform(reference_positions + shift, water)
    assert rmsd < 1e-6
    assert np.allclose(aligned, [[0.5, 0.5, 3.0]], atol=1e-6)


def test_aligner_rejects_mismatched_group_size():
    reference = AlignmentReference(
        indices=np.arange(3), positions=np.zeros((3, 3)), centroid=np.zeros(3), frame=0
    )
    with pytest.raises(ValueError, match="4 atoms but the reference has 3"):
        SiteAligner(reference).transform(np.zeros((4, 3)), np.zeros((4, 3)))


def test_build_reference_restores_the_current_frame(build_universe):
    universe = _protein_ligand(build_universe)
    ligand = universe.select_atoms("resname 547")
    group = select_alignment_group(universe, ligand, cutoff=20.0)
    universe.trajectory[2]

    reference = build_reference(group, frame=1)
    assert reference.frame == 1
    assert reference.n_atoms == 4
    assert np.allclose(reference.positions.mean(axis=0), 0.0, atol=1e-6)
    assert universe.trajectory.frame == 2


def test_rmsd_series_is_zero_for_a_rigid_system(build_universe):
    universe = _protein_ligand(build_universe)
    ligand = universe.select_atoms("resname 547")
    group = select_alignment_group(universe, ligand, cutoff=20.0)
    aligner = SiteAligner(build_reference(group))

    values = dict(iter_alignment_rmsd(universe, group, aligner))
    assert len(values) == 3
    assert max(values.values()) < 1e-5
