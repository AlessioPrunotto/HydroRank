"""Reproduce the apo/holo HSP90 validation reported in ``summary.md``.

The analysis deliberately uses the cached HydraRank observations.  It never loads the
multi-gigabyte trajectories after the first-frame reference systems have been prepared,
so peak memory is set by one compressed observations cache rather than by a trajectory.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.align import rotation_matrix
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from hydrarank.clustering import bulk_equivalent_count, cluster_positions
from hydrarank.config import PreprocessConfig
from hydrarank.data import WaterObservations
from hydrarank.entropy import (
    minus_t_delta_s,
    orientational_entropy,
    translational_entropy,
    water_orientations,
)
from hydrarank.io import load_universe
from hydrarank.preprocess import prepare_system

matplotlib.use("Agg", force=True)


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = ROOT / "output" / "hsp90"
SYSTEMS = ("apo", "3rlp", "3rlq", "3rlr")
HOLO = ("3rlp", "3rlq", "3rlr")
LIGANDS = {"3rlp": "3RP", "3rlq": "3RQ", "3rlr": "3RR"}
TARGET_WATERS = (("3rlp", 3, "W3"), ("3rlq", 249, "W249"), ("3rlq", 286, "W286"))
TEMPERATURE = 303.15
SITE_RADIUS = 1.0
CRYSTAL_MATCH_CUTOFF = 1.4


@dataclass(frozen=True)
class Motion:
    rotation: np.ndarray
    source_centroid: np.ndarray
    target_centroid: np.ndarray
    rmsd: float

    def apply(self, coordinates: np.ndarray) -> np.ndarray:
        coordinates = np.asarray(coordinates, dtype=float)
        return (coordinates - self.source_centroid) @ self.rotation.T + self.target_centroid


def fit_motion(source: np.ndarray, target: np.ndarray) -> Motion:
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    source_centroid = source.mean(axis=0)
    target_centroid = target.mean(axis=0)
    rotation, rmsd = rotation_matrix(
        source - source_centroid,
        target - target_centroid,
    )
    return Motion(rotation, source_centroid, target_centroid, float(rmsd))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def site_rows(system: str) -> list[dict[str, Any]]:
    rows = read_json(OUTPUT / system / "results.json")["ranking"]["sites"]
    return sorted(rows, key=lambda row: row["site"])


def site_centers(system: str) -> np.ndarray:
    return np.asarray([[row[key] for key in ("x", "y", "z")] for row in site_rows(system)])


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0])
    fieldnames.extend(key for row in rows[1:] for key in row if key not in fieldnames)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def prepared_systems():
    prepared = {}
    for system in SYSTEMS:
        config = PreprocessConfig.from_yaml(HERE / f"hsp90_{system}.yaml")
        prepared[system] = prepare_system(load_universe(config), config)
    return prepared


def ca_coordinates(prepared, system: str, resids: list[int]) -> np.ndarray:
    atoms = prepared[system].universe
    return np.asarray(
        [
            atoms.select_atoms(f"protein and name CA and resid {resid}").positions[0]
            for resid in resids
        ]
    )


def common_frame(prepared) -> tuple[dict[str, Motion], list[int]]:
    reference_resids = sorted({int(value) for value in prepared["3rlp"].fit_group.resids})
    common_resids = [
        resid
        for resid in reference_resids
        if all(
            prepared[name].universe.select_atoms(f"protein and name CA and resid {resid}").n_atoms
            == 1
            for name in SYSTEMS
        )
    ]
    reference = ca_coordinates(prepared, "3rlp", common_resids)
    motions = {
        system: fit_motion(ca_coordinates(prepared, system, common_resids), reference)
        for system in SYSTEMS
    }
    return motions, common_resids


def charmm_gui_text(system: str, suffix: str) -> str:
    system_dir = ROOT / "sample_traj" / system
    direct = list(system_dir.glob(f"charmm-gui-*/*{suffix}"))
    if direct:
        return direct[0].read_text()
    archive = system_dir / "charmm-gui.tgz"
    with tarfile.open(archive) as bundle:
        member = next(item for item in bundle.getmembers() if item.name.endswith(suffix))
        extracted = bundle.extractfile(member)
        if extracted is None:
            raise RuntimeError(f"could not read {member.name}")
        return extracted.read().decode()


def parse_pdb(text: str) -> list[dict[str, Any]]:
    atoms = []
    for line in io.StringIO(text):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms.append(
            {
                "record": line[:6].strip(),
                "name": line[12:16].strip(),
                "resname": line[17:21].strip(),
                "resid": int(line[22:26]),
                "coordinate": np.asarray(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                ),
                "bfactor": float(line[60:66]),
            }
        )
    return atoms


def parse_cif(system: str) -> list[dict[str, Any]]:
    atoms = []
    for line in charmm_gui_text(system, f"/{system}.cif").splitlines():
        fields = line.split()
        if len(fields) < 21 or fields[0] not in {"ATOM", "HETATM"}:
            continue
        atoms.append(
            {
                "record": fields[0],
                "element": fields[2].upper(),
                "name": fields[3],
                "resname": fields[5],
                "coordinate": np.asarray(fields[10:13], dtype=float),
                "bfactor": float(fields[14]),
                "resid": int(fields[16]),
                "chain": fields[18],
            }
        )
    return atoms


def crystal_to_simulation_motion(prepared, system: str, common_resids: list[int]) -> Motion:
    protein = parse_pdb(charmm_gui_text(system, f"/{system}_prob.pdb"))
    lookup = {
        (atom["resid"] - 8, atom["name"]): atom["coordinate"]
        for atom in protein
        if atom["name"] == "CA"
    }
    resids = [resid for resid in common_resids if (resid, "CA") in lookup]
    crystal = np.asarray([lookup[(resid, "CA")] for resid in resids])
    simulation = ca_coordinates(prepared, system, resids)
    return fit_motion(crystal, simulation)


def cif_chain_to_simulation_motion(
    prepared, system: str, common_resids: list[int], chain: str = "A"
) -> Motion:
    atoms = parse_cif(system)
    lookup = {
        (atom["resid"] - 8, atom["name"]): atom["coordinate"]
        for atom in atoms
        if atom["record"] == "ATOM" and atom["chain"] == chain and atom["name"] == "CA"
    }
    resids = [resid for resid in common_resids if (resid, "CA") in lookup]
    crystal = np.asarray([lookup[(resid, "CA")] for resid in resids])
    simulation = ca_coordinates(prepared, system, resids)
    return fit_motion(crystal, simulation)


def transform_crystal_system(
    prepared,
    common_motions: dict[str, Motion],
    common_resids: list[int],
    system: str,
) -> dict[str, Any]:
    crystal_motion = crystal_to_simulation_motion(prepared, system, common_resids)
    waters = parse_pdb(charmm_gui_text(system, f"/{system}_watb.pdb"))
    ligand = parse_pdb(charmm_gui_text(system, f"/{system}_hetc.pdb"))

    def transform(coordinates: np.ndarray) -> np.ndarray:
        return common_motions[system].apply(crystal_motion.apply(coordinates))

    ligand_coordinates = np.asarray([atom["coordinate"] for atom in ligand])
    water_coordinates = np.asarray([atom["coordinate"] for atom in waters])
    distances = np.linalg.norm(
        water_coordinates[:, None, :] - ligand_coordinates[None, :, :], axis=2
    ).min(axis=1)
    nearby = distances <= 5.0
    return {
        "ligand": transform(ligand_coordinates),
        "waters": transform(water_coordinates[nearby]),
        "water_atoms": [atom for atom, keep in zip(waters, nearby, strict=True) if keep],
        "fit_rmsd": crystal_motion.rmsd,
    }


def paper_target_coordinates(
    prepared,
    common_motions: dict[str, Motion],
    common_resids: list[int],
) -> list[dict[str, Any]]:
    output = []
    for system, resid, label in TARGET_WATERS:
        atoms = parse_cif(system)
        water = next(
            atom
            for atom in atoms
            if atom["record"] == "HETATM"
            and atom["resname"] == "HOH"
            and atom["chain"] == "A"
            and atom["resid"] == resid
        )
        crystal_motion = cif_chain_to_simulation_motion(prepared, system, common_resids, chain="A")
        coordinate = common_motions[system].apply(crystal_motion.apply(water["coordinate"]))
        output.append(
            {
                "label": label,
                "source_structure": system.upper(),
                "resid": resid,
                "coordinate": coordinate,
                "bfactor": water["bfactor"],
                "fit_rmsd": crystal_motion.rmsd,
            }
        )
    return output


def assign_fixed_sites(oxygen: np.ndarray, centers: np.ndarray, radius: float = 1.0):
    distance, label = cKDTree(centers).query(oxygen)
    return np.where(distance <= radius, label, -1), distance


def fixed_site_metrics(
    observations: WaterObservations,
    oxygen: np.ndarray,
    centers: np.ndarray,
    frame_mask: np.ndarray | None = None,
) -> list[dict[str, Any]]:
    labels, _ = assign_fixed_sites(oxygen, centers)
    if frame_mask is None:
        frame_mask = np.ones(observations.n_frames, dtype=bool)
    selected_frames = np.flatnonzero(frame_mask)
    frame_lookup = np.full(observations.n_frames, -1, dtype=int)
    frame_lookup[selected_frames] = np.arange(selected_frames.size)
    observation_mask = frame_mask[observations.frame]
    orientations = water_orientations(observations.oxygen, observations.hydrogen)
    rows = []
    for site in range(centers.shape[0]):
        members = observation_mask & (labels == site)
        occupied = np.unique(frame_lookup[observations.frame[members]])
        s_trans = translational_entropy(oxygen[members])
        s_orient = orientational_entropy(orientations[members])
        hb = float(np.mean(observations.hb_solute[members])) if np.any(members) else np.nan
        rows.append(
            {
                "site": site + 1,
                "n_observations": int(np.count_nonzero(members)),
                "occupancy": float(occupied.size / selected_frames.size),
                "s_trans": float(s_trans),
                "s_orient": float(s_orient),
                "minus_t_delta_s": float(minus_t_delta_s(s_trans + s_orient, TEMPERATURE)),
                "hb_solute": hb,
                "score": float(minus_t_delta_s(s_trans + s_orient, TEMPERATURE) - hb),
            }
        )
    return rows


def point_occupancy(
    observations: WaterObservations,
    oxygen: np.ndarray,
    point: np.ndarray,
    radius: float,
    frame_mask: np.ndarray | None = None,
) -> float:
    if frame_mask is None:
        frame_mask = np.ones(observations.n_frames, dtype=bool)
    members = np.linalg.norm(oxygen - point, axis=1) <= radius
    members &= frame_mask[observations.frame]
    occupied = np.unique(observations.frame[members])
    return float(occupied.size / np.count_nonzero(frame_mask))


def one_to_one_matches(
    reference: np.ndarray, predicted: np.ndarray, cutoff: float
) -> tuple[list[tuple[int, int, float]], int, int]:
    if reference.size == 0 or predicted.size == 0:
        return [], int(reference.shape[0]), int(predicted.shape[0])
    distances = np.linalg.norm(reference[:, None] - predicted[None, :], axis=2)
    ref_index, pred_index = linear_sum_assignment(distances)
    matches = [
        (int(i), int(j), float(distances[i, j]))
        for i, j in zip(ref_index, pred_index, strict=True)
        if distances[i, j] <= cutoff
    ]
    return matches, reference.shape[0] - len(matches), predicted.shape[0] - len(matches)


def crystal_validation(
    crystal: dict[str, dict[str, Any]],
    common_site_centers: dict[str, np.ndarray],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    match_rows = []
    metric_rows = []
    comparisons = [(name, name) for name in HOLO] + [("3rlp", "apo")]
    for crystal_name, prediction_name in comparisons:
        reference = crystal[crystal_name]["waters"]
        predicted = common_site_centers[prediction_name]
        for cutoff in (1.0, 1.4, 2.0):
            matches, n_missed, n_extra = one_to_one_matches(reference, predicted, cutoff)
            metric_rows.append(
                {
                    "crystal": crystal_name.upper(),
                    "prediction": prediction_name,
                    "cutoff_A": cutoff,
                    "n_crystal_waters": len(reference),
                    "n_sites": len(predicted),
                    "n_matched": len(matches),
                    "recall": len(matches) / len(reference),
                    "precision": len(matches) / len(predicted),
                    "n_missed": n_missed,
                    "n_extra": n_extra,
                }
            )
        matches, _, _ = one_to_one_matches(reference, predicted, CRYSTAL_MATCH_CUTOFF)
        for reference_index, prediction_index, distance in matches:
            atom = crystal[crystal_name]["water_atoms"][reference_index]
            match_rows.append(
                {
                    "crystal": crystal_name.upper(),
                    "prediction": prediction_name,
                    "water_resid": atom["resid"],
                    "water_bfactor_A2": atom["bfactor"],
                    "site": prediction_index + 1,
                    "distance_A": distance,
                }
            )
    return match_rows, metric_rows


def ligand_polar_contacts() -> list[dict[str, Any]]:
    rows = []
    target_resids = {51, 52, 93}
    for system in HOLO:
        atoms = parse_cif(system)
        ligand = [
            atom
            for atom in atoms
            if atom["record"] == "HETATM"
            and atom["resname"] == LIGANDS[system]
            and atom["chain"] == "A"
            and atom["element"] in {"N", "O", "S"}
        ]
        protein = [
            atom
            for atom in atoms
            if atom["record"] == "ATOM"
            and atom["chain"] == "A"
            and atom["resid"] in target_resids
            and atom["element"] in {"N", "O", "S"}
        ]
        for ligand_atom in ligand:
            for protein_atom in protein:
                distance = float(
                    np.linalg.norm(ligand_atom["coordinate"] - protein_atom["coordinate"])
                )
                if distance <= 3.5:
                    rows.append(
                        {
                            "structure": system.upper(),
                            "ligand": LIGANDS[system],
                            "ligand_atom": ligand_atom["name"],
                            "protein_residue": f"{protein_atom['resname']}{protein_atom['resid']}",
                            "protein_atom": protein_atom["name"],
                            "distance_A": distance,
                            "interpretation": (
                                "polar heavy-atom contact; hydrogen-bond geometry not assigned"
                            ),
                        }
                    )
    return rows


def block_convergence(
    observations: dict[str, WaterObservations],
    common_oxygen: dict[str, np.ndarray],
    apo_centers: np.ndarray,
    targets: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    site_rows_output = []
    target_rows_output = []
    for system in SYSTEMS:
        data = observations[system]
        blocks = np.array_split(np.arange(data.n_frames), 5)
        for block_index, frames in enumerate(blocks, start=1):
            frame_mask = np.zeros(data.n_frames, dtype=bool)
            frame_mask[frames] = True
            metrics = fixed_site_metrics(
                data, common_oxygen[system], apo_centers, frame_mask=frame_mask
            )
            for row in metrics:
                site_rows_output.append(
                    {
                        "system": system,
                        "block": block_index,
                        "start_ns": data.times[frames[0]] / 1000.0,
                        "stop_ns": data.times[frames[-1]] / 1000.0,
                        **row,
                    }
                )
            for target in targets:
                for radius in (1.0, 1.4):
                    target_rows_output.append(
                        {
                            "system": system,
                            "block": block_index,
                            "start_ns": data.times[frames[0]] / 1000.0,
                            "stop_ns": data.times[frames[-1]] / 1000.0,
                            "target": target["label"],
                            "radius_A": radius,
                            "occupancy": point_occupancy(
                                data,
                                common_oxygen[system],
                                target["coordinate"],
                                radius,
                                frame_mask=frame_mask,
                            ),
                        }
                    )
    return site_rows_output, target_rows_output


def sensitivity_analysis(
    apo_observations: WaterObservations,
    apo_oxygen: np.ndarray,
    apo_rows: list[dict[str, Any]],
    crystal_waters: np.ndarray,
    targets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows = []
    settings = [
        ("site_radius", 0.8, 2.0, 0.0),
        ("site_radius", 1.0, 2.0, 0.0),
        ("site_radius", 1.2, 2.0, 0.0),
        ("density_factor", 1.0, 1.5, 0.0),
        ("density_factor", 1.0, 2.5, 0.0),
        ("discard_ns", 1.0, 2.0, 5.0),
        ("discard_ns", 1.0, 2.0, 10.0),
    ]
    for parameter, radius, density, discard_ns in settings:
        keep_frames = apo_observations.times >= discard_ns * 1000.0
        keep_observations = keep_frames[apo_observations.frame]
        n_frames = int(np.count_nonzero(keep_frames))
        min_count = bulk_equivalent_count(n_frames, radius, density)
        centers, _ = cluster_positions(
            apo_oxygen[keep_observations], radius=radius, min_count=min_count
        )
        matches, _, _ = one_to_one_matches(crystal_waters, centers, CRYSTAL_MATCH_CUTOFF)
        target_distances = [
            float(np.linalg.norm(centers - target["coordinate"], axis=1).min())
            for target in targets
        ]
        rows.append(
            {
                "analysis": "clustering",
                "parameter": parameter,
                "value": (
                    radius
                    if parameter == "site_radius"
                    else density
                    if parameter == "density_factor"
                    else discard_ns
                ),
                "n_sites": len(centers),
                "crystal_matches_at_1.4_A": len(matches),
                **{
                    f"{target['label']}_nearest_site_A": distance
                    for target, distance in zip(targets, target_distances, strict=True)
                },
            }
        )

    baseline_scores = np.asarray([row["score"] for row in apo_rows])
    baseline_rank = np.argsort(np.argsort(-baseline_scores))
    for penalty in (0.5, 1.0, 1.5):
        scores = np.asarray(
            [row["minus_t_delta_s"] - penalty * row["hb_solute"] for row in apo_rows]
        )
        ranks = np.argsort(np.argsort(-scores))
        rows.append(
            {
                "analysis": "ranking",
                "parameter": "hbond_penalty",
                "value": penalty,
                "spearman_vs_default": float(spearmanr(baseline_rank, ranks).statistic),
                "top_site": int(np.argmax(scores) + 1),
            }
        )
    return rows


def make_figures(
    common_rows: list[dict[str, Any]],
    target_rows: list[dict[str, Any]],
    block_targets: list[dict[str, Any]],
    sensitivity: list[dict[str, Any]],
) -> None:
    figure_dir = HERE / "figures"
    figure_dir.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    target_labels = [target[2] for target in TARGET_WATERS]
    x = np.arange(len(target_labels))
    width = 0.19
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for index, system in enumerate(SYSTEMS):
        values = [
            next(
                row["occupancy"]
                for row in target_rows
                if row["system"] == system and row["target"] == target and row["radius_A"] == 1.4
            )
            for target in target_labels
        ]
        ax.bar(x + (index - 1.5) * width, values, width, label=system.upper())
    ax.set_xticks(x, target_labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fixed crystallographic-site occupancy")
    ax.legend(ncol=4, frameon=False)
    ax.set_title("Occupancy of crystallographic waters across HSP90 systems", pad=22)
    ax.text(
        0.5,
        1.01,
        "Fraction of frames with a water oxygen within 1.4 Å of the experimental position",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8,
        color="0.35",
    )
    fig.text(
        0.5,
        0.01,
        "Bars represent separate 50 ns simulations; absent bars correspond to zero occupancy.",
        ha="center",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(figure_dir / "target_water_occupancy.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    for axis, target in zip(axes, target_labels, strict=True):
        for system in SYSTEMS:
            rows = [
                row
                for row in block_targets
                if row["system"] == system and row["target"] == target and row["radius_A"] == 1.4
            ]
            axis.plot(
                [row["block"] for row in rows],
                [row["occupancy"] for row in rows],
                "o-",
                label=system.upper(),
            )
        axis.set_title(target)
        axis.set_xlabel("10 ns block")
        axis.set_ylim(-0.03, 1.03)
    axes[0].set_ylabel("Occupancy")
    axes[-1].legend(frameon=False, fontsize=8)
    fig.suptitle("Crystallographic-water occupancy in consecutive 10 ns blocks", y=0.99)
    fig.text(
        0.5,
        0.92,
        "Occupancy uses the same fixed 1.4 Å sphere in every block and system",
        ha="center",
        fontsize=8,
        color="0.35",
    )
    fig.text(
        0.5,
        0.01,
        "Block-to-block variation assesses within-trajectory stability, "
        "not replica reproducibility.",
        ha="center",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.89))
    fig.savefig(figure_dir / "target_water_block_convergence.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    systems = list(SYSTEMS)
    highlighted_sites = {2, 4, 9}
    other_sites_labelled = False
    for site in range(1, 20):
        if site in highlighted_sites:
            continue
        values = [
            next(
                row["occupancy"]
                for row in common_rows
                if row["system"] == system and row["site"] == site
            )
            for system in systems
        ]
        ax.plot(
            systems,
            values,
            color="0.75",
            linewidth=0.8,
            label="Other apo sites" if not other_sites_labelled else None,
            zorder=1,
        )
        other_sites_labelled = True
    highlighted = (
        (9, "#d62728", "W3 ↔ apo site 9 (0.32 Å)", "-", "#d62728"),
        (2, "#1f77b4", "W249 ↔ apo site 2 (0.46 Å)", "-", "#1f77b4"),
        (
            4,
            "#2ca02c",
            "Nearest apo site to W286: site 4 (1.50 Å; outside 1.4 Å cutoff)",
            "--",
            "white",
        ),
    )
    for site, color, label, linestyle, marker_face in highlighted:
        values = [
            next(
                row["occupancy"]
                for row in common_rows
                if row["system"] == system and row["site"] == site
            )
            for system in systems
        ]
        ax.plot(
            systems,
            values,
            marker="o",
            linestyle=linestyle,
            linewidth=2.3,
            color=color,
            markerfacecolor=marker_face,
            markeredgecolor=color,
            markeredgewidth=1.5,
            label=label,
            zorder=2,
        )
    ax.set_title("Apo-derived hydration-site occupancy across HSP90 systems", pad=22)
    ax.text(
        0.5,
        1.01,
        "Highlighted labels report the distance from each crystallographic water "
        "to the apo-site center",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8,
        color="0.35",
    )
    ax.set_ylabel("Fixed apo-site occupancy")
    ax.set_ylim(-0.03, 1.03)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.text(
        0.5,
        0.01,
        "Occupancy is measured within 1.0 Å of fixed apo-derived centers; "
        "x-axis categories are separate simulations, not time.",
        ha="center",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(figure_dir / "apo_site_occupancy_series.png", dpi=220)
    plt.close(fig)

    clustering = [row for row in sensitivity if row["analysis"] == "clustering"]
    fig, ax = plt.subplots(figsize=(9.4, 5.0))
    positions = np.asarray([0.0, 1.0, 2.0, 3.6, 4.6, 6.2, 7.2])
    labels = [
        "Radius\n0.8 Å",
        "Radius\n1.0 Å\n(default)",
        "Radius\n1.2 Å",
        "Density factor\n1.5",
        "Density factor\n2.5",
        "Discard\n5 ns",
        "Discard\n10 ns",
    ]
    colors = ["#4c78a8"] * 3 + ["#f58518"] * 2 + ["#54a24b"] * 2
    counts = [row["n_sites"] for row in clustering]
    bars = ax.bar(positions, counts, color=colors, width=0.82)
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.axhline(
        19,
        color="0.25",
        linestyle="--",
        linewidth=1,
        label="Default result: 19 sites",
    )
    ax.set_xticks(positions, labels)
    ax.set_ylabel("Apo hydration sites")
    ax.set_ylim(0, max(counts) + 4)
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    ax.set_title("Sensitivity of the apo hydration-site count", pad=22)
    ax.text(
        0.5,
        1.01,
        "One parameter varied at a time; defaults: radius 1.0 Å, density factor 2.0, discard 0 ns",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8,
        color="0.35",
    )
    fig.text(
        0.5,
        0.01,
        "Bars report the number of detected sites only—not accuracy, occupancy, "
        "or ranking quality.",
        ha="center",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(figure_dir / "sensitivity_site_count.png", dpi=220)
    plt.close(fig)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    prepared = prepared_systems()
    common_motions, common_resids = common_frame(prepared)
    observations = {
        system: WaterObservations.load(OUTPUT / system / "observations.npz") for system in SYSTEMS
    }
    common_oxygen = {
        system: common_motions[system].apply(data.oxygen) for system, data in observations.items()
    }
    common_centers = {
        system: common_motions[system].apply(site_centers(system)) for system in SYSTEMS
    }
    common_ligands = {
        system: common_motions[system].apply(prepared[system].ligand_reference)
        for system in SYSTEMS
    }
    apo_centers = common_centers["apo"]
    apo_rows = site_rows("apo")

    crystal = {
        system: transform_crystal_system(prepared, common_motions, common_resids, system)
        for system in HOLO
    }
    targets = paper_target_coordinates(prepared, common_motions, common_resids)

    fixed_metrics = {
        system: fixed_site_metrics(observations[system], common_oxygen[system], apo_centers)
        for system in SYSTEMS
    }
    common_rows = []
    for site_index, apo_row in enumerate(apo_rows):
        for system in SYSTEMS:
            ligand_distance = float(
                np.linalg.norm(common_ligands[system] - apo_centers[site_index], axis=1).min()
            )
            metric = fixed_metrics[system][site_index]
            ratio = metric["occupancy"] / max(fixed_metrics["apo"][site_index]["occupancy"], 1e-12)
            if system == "apo":
                state = "apo reference"
            elif ligand_distance < 2.2:
                state = "sterically occupied"
            elif (
                ratio < 0.25
                and fixed_metrics["apo"][site_index]["occupancy"] - metric["occupancy"] > 0.2
            ):
                state = "depleted"
            elif ratio >= 0.65:
                state = "retained"
            else:
                state = "perturbed"
            common_rows.append(
                {
                    "site": site_index + 1,
                    "system": system,
                    "apo_rank": apo_row["rank"],
                    "apo_category": apo_row["category"],
                    "apo_score": apo_row["score"],
                    "occupancy": metric["occupancy"],
                    "occupancy_ratio_to_apo": ratio,
                    "ligand_distance_A": ligand_distance,
                    "full_shell_coverage": ligand_distance <= 4.0,
                    "state": state,
                }
            )

    target_rows = []
    for target in targets:
        target_distances = np.linalg.norm(apo_centers - target["coordinate"], axis=1)
        nearest_site = int(np.argmin(target_distances))
        for system in SYSTEMS:
            ligand_distance = float(
                np.linalg.norm(common_ligands[system] - target["coordinate"], axis=1).min()
            )
            for radius in (1.0, 1.4):
                target_rows.append(
                    {
                        "target": target["label"],
                        "source_structure": target["source_structure"],
                        "water_bfactor_A2": target["bfactor"],
                        "nearest_apo_site": nearest_site + 1,
                        "nearest_apo_site_distance_A": float(target_distances[nearest_site]),
                        "system": system,
                        "radius_A": radius,
                        "occupancy": point_occupancy(
                            observations[system],
                            common_oxygen[system],
                            target["coordinate"],
                            radius,
                        ),
                        "ligand_distance_A": ligand_distance,
                    }
                )

    crystal_matches, crystal_metrics = crystal_validation(crystal, common_centers)
    block_sites, block_targets = block_convergence(
        observations, common_oxygen, apo_centers, targets
    )
    sensitivity = sensitivity_analysis(
        observations["apo"],
        common_oxygen["apo"],
        apo_rows,
        crystal["3rlp"]["waters"],
        targets,
    )
    contacts = ligand_polar_contacts()

    write_csv(HERE / "common_apo_sites.csv", common_rows)
    write_csv(HERE / "paper_target_waters.csv", target_rows)
    write_csv(HERE / "crystal_water_matches.csv", crystal_matches)
    write_csv(HERE / "crystal_recovery_metrics.csv", crystal_metrics)
    write_csv(HERE / "block_site_metrics.csv", block_sites)
    write_csv(HERE / "block_target_occupancy.csv", block_targets)
    write_csv(HERE / "sensitivity.csv", sensitivity)
    write_csv(HERE / "ligand_polar_contacts.csv", contacts)
    make_figures(common_rows, target_rows, block_targets, sensitivity)

    incremental = []
    for system in HOLO:
        occupied = [
            row
            for row in common_rows
            if row["system"] == system and row["state"] == "sterically occupied"
        ]
        incremental.append(
            {
                "system": system.upper(),
                "n_sterically_occupied_apo_sites": len(occupied),
                "sum_apo_scores": float(sum(row["apo_score"] for row in occupied)),
                "sites": [row["site"] for row in occupied],
            }
        )

    metrics = {
        "temperature_K": TEMPERATURE,
        "common_frame": {
            "reference": "3RLP HydraRank frame",
            "n_ca_atoms": len(common_resids),
            "fit_rmsd_A": {system: common_motions[system].rmsd for system in SYSTEMS},
        },
        "crystal_to_simulation_fit_rmsd_A": {
            system: crystal[system]["fit_rmsd"] for system in HOLO
        },
        "paper_targets": [
            {
                key: value.tolist() if isinstance(value, np.ndarray) else value
                for key, value in target.items()
            }
            for target in targets
        ],
        "incremental_steric_coverage": incremental,
        "limitations": [
            "One 50 ns trajectory per state; no independent replicas.",
            "The apo trajectory retained crystallographic waters at initialization.",
            "Crystal-water occupancy is not a thermodynamic ground truth.",
            "W286 is not robustly recovered by the conventional 3RLQ trajectory.",
            "Three complexes support a qualitative series comparison, not a fitted affinity model.",
        ],
    }
    (HERE / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    tracked = [
        HERE / "run_validation.py",
        *[HERE / f"hsp90_{system}.yaml" for system in SYSTEMS],
        *[OUTPUT / system / "results.json" for system in SYSTEMS],
    ]
    manifest = {
        "files": [
            {
                "path": str(path.relative_to(ROOT)),
                "size": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in tracked
        ]
    }
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
