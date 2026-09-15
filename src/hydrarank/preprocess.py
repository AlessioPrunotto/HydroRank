"""Turn a raw trajectory into aligned water observations around the ligand.

The shell is defined against the *reference* ligand pose in the aligned frame, not
against the moving ligand: a region that follows the ligand's wobble would make site
occupancies depend on ligand motion rather than on water behaviour.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.lib.distances import capped_distance

from hydrarank.alignment import SiteAligner, build_reference, select_alignment_group
from hydrarank.config import PreprocessConfig
from hydrarank.data import WaterObservations
from hydrarank.hbonds import (
    PolarGroups,
    count_hbonds,
    count_neighbours,
    find_polar_groups,
)
from hydrarank.io import frame_slice
from hydrarank.pbc import (
    PBCGroups,
    apply_pbc_transformations,
    build_pbc_groups,
    max_bond_length,
    validate_cutoff,
)
from hydrarank.selections import (
    WaterTopology,
    analyse_water_topology,
    classify_atoms,
    select_ligand,
    select_pocket,
    select_water,
)


@dataclass
class PreparedSystem:
    """A universe with the PBC and alignment machinery in place."""

    universe: Universe
    ligand: AtomGroup
    ligand_heavy: AtomGroup
    water: WaterTopology
    oxygens: AtomGroup
    groups: PBCGroups
    fit_group: AtomGroup
    aligner: SiteAligner
    frames: slice
    reference_frame: int
    ligand_reference: np.ndarray
    polar: PolarGroups

    @property
    def n_selected_frames(self) -> int:
        return len(range(*self.frames.indices(self.universe.trajectory.n_frames)))


def prepare_system(universe: Universe, config: PreprocessConfig) -> PreparedSystem:
    """Run the selections, attach the PBC transformations and build the fit reference."""
    ligand = select_ligand(
        universe,
        config.ligand_selection,
        water_resnames=config.water_resnames,
        ion_resnames=config.ion_resnames,
        min_heavy_atoms=config.min_ligand_heavy_atoms,
    )
    ligand_roles = classify_atoms(ligand)
    ligand_heavy = ligand[(ligand_roles != "H") & (ligand_roles != "M")]
    water = analyse_water_topology(select_water(universe, config.water_resnames))

    for name, cutoff in (
        ("water_cutoff", config.water_cutoff),
        ("pocket_cutoff", config.pocket_cutoff),
        ("hbond_distance", config.hbond_distance),
        ("enclosure_radius", config.enclosure_radius),
    ):
        validate_cutoff(universe, cutoff, name=name)
    groups = build_pbc_groups(universe, ligand)
    apply_pbc_transformations(universe, groups)

    frames = frame_slice(config, universe.trajectory.n_frames)
    reference_frame = frames.start
    universe.trajectory[reference_frame]
    fit_group = select_alignment_group(universe, ligand, config.pocket_cutoff)
    aligner = SiteAligner(build_reference(fit_group, frame=reference_frame))
    solute = select_pocket(universe, ligand, config.pocket_cutoff) | ligand

    return PreparedSystem(
        universe=universe,
        ligand=ligand,
        ligand_heavy=ligand_heavy,
        water=water,
        oxygens=universe.atoms[water.oxygen_ix],
        groups=groups,
        fit_group=fit_group,
        aligner=aligner,
        frames=frames,
        reference_frame=reference_frame,
        ligand_reference=ligand_heavy.positions.astype(np.float64, copy=True),
        polar=find_polar_groups(solute),
    )


def run_preprocess(
    universe: Universe,
    config: PreprocessConfig,
    progress: Callable[[int, int], None] | None = None,
    collect_qc: bool = False,
) -> WaterObservations:
    """Extract aligned first-shell waters, optionally collecting QC in the same pass."""
    system = prepare_system(universe, config)
    hydrogens = system.universe.atoms[system.water.hydrogen_ix.ravel()]

    oxygen_chunks, hydrogen_chunks, frame_chunks, resid_chunks, water_id_chunks = (
        [],
        [],
        [],
        [],
        [],
    )
    hbond_chunks, enclosure_chunks = [], []
    source_frames, times, rmsds = [], [], []
    longest_bonds, shell_counts = [], []

    for index, ts in enumerate(system.universe.trajectory[system.frames]):
        motion = system.aligner.fit(system.fit_group.positions)
        aligned_oxygen = motion.apply(system.oxygens.positions)

        pairs = capped_distance(
            aligned_oxygen,
            system.ligand_reference,
            max_cutoff=config.water_cutoff,
            return_distances=False,
        )
        selected = np.unique(pairs[:, 0]) if len(pairs) else np.empty(0, dtype=int)

        source_frames.append(ts.frame)
        times.append(float(ts.time))
        rmsds.append(motion.rmsd)
        if collect_qc:
            longest_bonds.append(max_bond_length(system.groups.solute))
            shell_counts.append(int(selected.size))
        if progress is not None:
            progress(index + 1, system.n_selected_frames)
        if selected.size == 0:
            continue

        raw_oxygen = system.oxygens.positions[selected]
        raw_hydrogen = hydrogens.positions.reshape(-1, 2, 3)[selected]
        hbond_chunks.append(
            count_hbonds(
                raw_oxygen,
                raw_hydrogen,
                system.polar,
                distance=config.hbond_distance,
                angle=config.hbond_angle,
            )
        )
        enclosure_chunks.append(
            count_neighbours(raw_oxygen, system.polar.heavy, radius=config.enclosure_radius)
        )

        oxygen_chunks.append(aligned_oxygen[selected])
        hydrogen_chunks.append(motion.apply(raw_hydrogen.reshape(-1, 3)).reshape(-1, 2, 3))
        frame_chunks.append(np.full(selected.size, index, dtype=np.int32))
        resid_chunks.append(system.water.resids[selected])
        water_id_chunks.append(system.water.resindices[selected])

    return WaterObservations(
        oxygen=_concat(oxygen_chunks, (0, 3)),
        hydrogen=_concat(hydrogen_chunks, (0, 2, 3)),
        frame=_concat(frame_chunks, (0,), dtype=np.int32),
        resid=_concat(resid_chunks, (0,), dtype=np.int64),
        water_id=_concat(water_id_chunks, (0,), dtype=np.int64),
        frames=np.asarray(source_frames, dtype=np.int64),
        times=np.asarray(times, dtype=np.float64),
        fit_rmsd=np.asarray(rmsds, dtype=np.float64),
        ligand_reference=system.ligand_reference,
        hb_solute=_concat(hbond_chunks, (0,)),
        enclosure=_concat(enclosure_chunks, (0,)),
        metadata={
            "water_cutoff": config.water_cutoff,
            "pocket_cutoff": config.pocket_cutoff,
            "hbond_distance": config.hbond_distance,
            "hbond_angle": config.hbond_angle,
            "enclosure_radius": config.enclosure_radius,
            "reference_frame": system.reference_frame,
            "n_fit_atoms": int(system.fit_group.n_atoms),
            "ligand_selection": config.ligand_selection,
            "topology": str(config.topology),
            "trajectory": [str(p) for p in config.trajectory],
            **(
                {
                    "qc_longest_bond": longest_bonds,
                    "qc_n_shell_waters": shell_counts,
                }
                if collect_qc
                else {}
            ),
        },
    )


def _concat(chunks: list, empty_shape: tuple, dtype=np.float64) -> np.ndarray:
    if not chunks:
        return np.empty(empty_shape, dtype=dtype)
    return np.concatenate(chunks).astype(dtype, copy=False)
