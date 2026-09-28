#!/usr/bin/env python3
# ruff: noqa: E501
"""Reproducible HydraRank/SSTMap comparison for the harmonized 3RLP run."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
import MDAnalysis as mda

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.align import rotation_matrix
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr, spearmanr

from hydrarank.data import WaterObservations
from hydrarank.selections import classify_atoms

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
FIG = OUT / "figures"
HYDRA_JSON = ROOT / "output/hsp90/3rlp/results.json"
HYDRA_OBS = ROOT / "output/hsp90/3rlp/observations.npz"
SST_LOCAL = (
    ROOT / "sstmap_analysis/3rlp/full_hydrarank61/SSTMap_HSA/3rlp_hydrarank61_hsa_summary.txt"
)
SST_GLOBAL = (
    ROOT / "sstmap_analysis/3rlp/full_opc_neighbors/SSTMap_HSA/3rlp_opc_fixed_hsa_summary.txt"
)
SST_LIGAND = ROOT / "sstmap_analysis/3rlp/input/3RP_aligned.pdb"
MATCH_CUTOFF = 1.0
R_KCAL = 1.987204259e-3
TEMPERATURE = 303.15
SSTMAP_TEMPERATURE = 300.0

SST_NUMERIC_COLUMNS = [
    "index",
    "x",
    "y",
    "z",
    "nwat",
    "occupancy",
    "Esw",
    "EswLJ",
    "EswElec",
    "Eww",
    "EwwLJ",
    "EwwElec",
    "Etot",
    "Ewwnbr",
    "TSsw_trans",
    "TSsw_orient",
    "TStot",
    "Nnbrs",
    "Nhbww",
    "Nhbsw",
    "Nhbtot",
    "f_hb_ww",
    "f_enc",
    "Acc_ww",
    "Don_ww",
    "Acc_sw",
    "Don_sw",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def load_sst(path: Path) -> dict[str, np.ndarray]:
    values = np.loadtxt(path, skiprows=1, usecols=range(len(SST_NUMERIC_COLUMNS)), ndmin=2)
    return {name: values[:, i] for i, name in enumerate(SST_NUMERIC_COLUMNS)}


def ligand_transform() -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return SST-to-Hydra rotation/centroids from the shared ligand heavy atoms."""
    target = WaterObservations.load(HYDRA_OBS).ligand_reference.astype(float)
    universe = mda.Universe(str(SST_LIGAND))
    ligand = universe.select_atoms("resname 3RP")
    roles = classify_atoms(ligand)
    source = ligand.positions[(roles != "H") & (roles != "M")].astype(float)
    if source.shape != target.shape:
        raise RuntimeError(f"ligand heavy-atom mismatch: SST {source.shape}, Hydra {target.shape}")
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    rotation, _ = rotation_matrix(source - source_center, target - target_center)
    transformed = (source - source_center) @ rotation.T + target_center
    rmsd = float(np.sqrt(np.mean(np.sum((transformed - target) ** 2, axis=1))))
    return rotation, source_center, target_center, rmsd


def transform(
    coords: np.ndarray, params: tuple[np.ndarray, np.ndarray, np.ndarray, float]
) -> np.ndarray:
    rotation, source_center, target_center, _ = params
    return (coords - source_center) @ rotation.T + target_center


def match_sites(
    a: np.ndarray, b: np.ndarray, cutoff: float
) -> tuple[list[tuple[int, int, float]], list[int], list[int]]:
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    rows, cols = linear_sum_assignment(distances)
    pairs = [
        (int(i), int(j), float(distances[i, j]))
        for i, j in zip(rows, cols, strict=True)
        if distances[i, j] <= cutoff
    ]
    matched_a = {i for i, _, _ in pairs}
    matched_b = {j for _, j, _ in pairs}
    return pairs, sorted(set(range(len(a))) - matched_a), sorted(set(range(len(b))) - matched_b)


def bootstrap_mean_ci(values: np.ndarray, seed: int = 20260918) -> list[float]:
    rng = np.random.default_rng(seed)
    samples = rng.choice(values, size=(10_000, values.size), replace=True).mean(axis=1)
    return [float(x) for x in np.percentile(samples, [2.5, 97.5])]


def comparison_stats(a: np.ndarray, b: np.ndarray) -> dict[str, float | list[float]]:
    delta = a - b
    return {
        "n": int(a.size),
        "pearson_r": float(pearsonr(a, b).statistic),
        "spearman_rho": float(spearmanr(a, b).statistic),
        "mean_hydrarank_minus_sstmap": float(delta.mean()),
        "mean_difference_across_site_bootstrap_95pct_ci": bootstrap_mean_ci(delta),
        "mae": float(np.mean(np.abs(delta))),
        "rmse": float(np.sqrt(np.mean(delta**2))),
    }


def scatter(ax, x, y, xlabel, ylabel, title):
    ax.scatter(x, y, s=45, color="#2266aa", edgecolor="white", linewidth=0.6, zorder=3)
    low = min(float(np.min(x)), float(np.min(y)))
    high = max(float(np.max(x)), float(np.max(y)))
    margin = max((high - low) * 0.08, 0.03)
    ax.plot([low - margin, high + margin], [low - margin, high + margin], "--", color="0.45", lw=1)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(alpha=0.2)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    hydra_doc = json.loads(HYDRA_JSON.read_text())
    hydra = sorted(hydra_doc["ranking"]["sites"], key=lambda row: int(row["site"]))
    hxyz = np.array([[row["x"], row["y"], row["z"]] for row in hydra], dtype=float)
    local = load_sst(SST_LOCAL)
    global_ = load_sst(SST_GLOBAL)
    tfm = ligand_transform()
    lxyz = transform(np.column_stack([local["x"], local["y"], local["z"]]), tfm)
    gxyz = transform(np.column_stack([global_["x"], global_["y"], global_["z"]]), tfm)

    pairs, unmatched_h, unmatched_l = match_sites(hxyz, lxyz, MATCH_CUTOFF)
    if len(pairs) != len(hydra):
        raise RuntimeError("not every HydraRank site matched the harmonized SSTMap run")

    hidx = np.array([i for i, _, _ in pairs])
    lidx = np.array([j for _, j, _ in pairs])
    distances = np.array([d for _, _, d in pairs])
    h_occ = np.array([hydra[i]["occupancy"] for i in hidx])
    h_trans = -R_KCAL * TEMPERATURE * np.array([hydra[i]["s_trans"] for i in hidx])
    h_orient = -R_KCAL * TEMPERATURE * np.array([hydra[i]["s_orient"] for i in hidx])
    h_total = np.array([hydra[i]["minus_t_delta_s"] for i in hidx])
    h_hbond = np.array([hydra[i]["hb_solute"] for i in hidx])
    s_occ = local["occupancy"][lidx]
    temperature_scale = TEMPERATURE / SSTMAP_TEMPERATURE
    s_trans = -local["TSsw_trans"][lidx] * temperature_scale
    s_orient = -local["TSsw_orient"][lidx] * temperature_scale
    s_total = -local["TStot"][lidx] * temperature_scale
    s_hbond = local["Nhbsw"][lidx]
    # SSTMap has no ranking. This proxy deliberately applies HydraRank's published
    # score form (entropy cost - 1 kcal/mol per solute H-bond) to SSTMap observables.
    sst_score_proxy = s_total - s_hbond

    rows = []
    for k, (hi, si, distance) in enumerate(pairs):
        h = hydra[hi]
        rows.append(
            {
                "hydrarank_site": int(h["site"]),
                "hydrarank_rank": int(h["rank"]),
                "sstmap_site": int(local["index"][si]),
                "center_distance_A": distance,
                "hydrarank_occupancy": h_occ[k],
                "sstmap_occupancy": s_occ[k],
                "hydrarank_minus_TdS_trans_kcal_mol": h_trans[k],
                "sstmap_minus_TdS_trans_kcal_mol": s_trans[k],
                "hydrarank_minus_TdS_orient_kcal_mol": h_orient[k],
                "sstmap_minus_TdS_orient_kcal_mol": s_orient[k],
                "hydrarank_minus_TdS_total_kcal_mol": h_total[k],
                "sstmap_minus_TdS_total_kcal_mol": s_total[k],
                "hydrarank_solute_hbonds": h_hbond[k],
                "sstmap_solute_hbonds": s_hbond[k],
                "hydrarank_residence_ps": h["mean_residence_ps"],
                "hydrarank_enclosure": h["enclosure"],
                "hydrarank_score": h["score"],
                "sstmap_hydrarank_like_score_proxy": sst_score_proxy[k],
                "hydrarank_category": h["category"],
                "sstmap_Esw_kcal_mol": local["Esw"][si],
                "sstmap_Eww_kcal_mol": local["Eww"][si],
                "sstmap_Etot_kcal_mol": local["Etot"][si],
                "sstmap_water_neighbors": local["Nnbrs"][si],
                "sstmap_water_water_hbonds": local["Nhbww"][si],
            }
        )
    with (OUT / "matched_sites.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    unmatched_rows = []
    for i in unmatched_h:
        nearest = float(np.min(np.linalg.norm(lxyz - hxyz[i], axis=1)))
        unmatched_rows.append(
            {"method": "HydraRank", "site": int(hydra[i]["site"]), "nearest_other_A": nearest}
        )
    for j in unmatched_l:
        nearest = float(np.min(np.linalg.norm(hxyz - lxyz[j], axis=1)))
        unmatched_rows.append(
            {"method": "SSTMap", "site": int(local["index"][j]), "nearest_other_A": nearest}
        )
    with (OUT / "unmatched_sites.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["method", "site", "nearest_other_A"])
        writer.writeheader()
        writer.writerows(unmatched_rows)

    # Sensitivity of SSTMap results to the alignment choice.
    lg_pairs, _, lg_unmatched = match_sites(lxyz, gxyz, MATCH_CUTOFF)
    local_i = np.array([i for i, _, _ in lg_pairs])
    global_i = np.array([j for _, j, _ in lg_pairs])
    lg_dist = np.array([d for _, _, d in lg_pairs])
    alignment_rows = []
    for i, j, d in lg_pairs:
        alignment_rows.append(
            {
                "local_sstmap_site": int(local["index"][i]),
                "global_sstmap_site": int(global_["index"][j]),
                "center_shift_A": d,
                "local_occupancy": local["occupancy"][i],
                "global_occupancy": global_["occupancy"][j],
                "local_minus_TdS_total_kcal_mol": -local["TStot"][i] * temperature_scale,
                "global_minus_TdS_total_kcal_mol": -global_["TStot"][j] * temperature_scale,
            }
        )
    with (OUT / "alignment_sensitivity.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(alignment_rows[0]))
        writer.writeheader()
        writer.writerows(alignment_rows)

    fields = {
        "occupancy": (h_occ, s_occ),
        "minus_T_delta_S_trans_kcal_mol": (h_trans, s_trans),
        "minus_T_delta_S_orient_kcal_mol": (h_orient, s_orient),
        "minus_T_delta_S_total_kcal_mol": (h_total, s_total),
        "solute_hbonds_native_definitions": (h_hbond, s_hbond),
        "hydrarank_score_vs_sstmap_derived_score_proxy": (
            np.array([hydra[i]["score"] for i in hidx]),
            sst_score_proxy,
        ),
    }
    hydra_scores = np.array([hydra[i]["score"] for i in hidx])
    ranking_overlap = {}
    for k in (3, 5, 10):
        hydra_top = set(np.argsort(-hydra_scores)[:k])
        proxy_top = set(np.argsort(-sst_score_proxy)[:k])
        ranking_overlap[f"top_{k}"] = {
            "overlap": len(hydra_top & proxy_top),
            "possible": k,
            "jaccard": len(hydra_top & proxy_top) / len(hydra_top | proxy_top),
        }
    metrics = {
        "scope": "single-system controlled method comparison; not general validation",
        "settings": {
            "temperature_K": TEMPERATURE,
            "sstmap_native_entropy_temperature_K": SSTMAP_TEMPERATURE,
            "bulk_density_A^-3": 0.0333,
            "frames": 25001,
            "frame_spacing_ps": 2.0,
            "site_radius_A": 1.0,
            "matching_cutoff_A": MATCH_CUTOFF,
        },
        "site_detection": {
            "hydrarank_sites": len(hydra),
            "sstmap_sites": len(local["index"]),
            "matched_pairs": len(pairs),
            "hydrarank_unmatched": len(unmatched_h),
            "sstmap_unmatched": len(unmatched_l),
            "matched_fraction_of_union": len(pairs)
            / (len(pairs) + len(unmatched_h) + len(unmatched_l)),
            "center_distance_A": {
                "mean": float(distances.mean()),
                "median": float(np.median(distances)),
                "max": float(distances.max()),
            },
            "matches_at_cutoffs_A": {
                str(c): len(match_sites(hxyz, lxyz, c)[0]) for c in (0.25, 0.5, 0.75, 1.0)
            },
        },
        "matched_site_metrics": {
            name: comparison_stats(*arrays) for name, arrays in fields.items()
        },
        "ranking_comparison": {
            "definition": "SSTMap proxy = SSTMap -TDeltaS total - 1 kcal/mol * SSTMap solute H-bonds; this is not a native SSTMap score",
            "top_k_overlap": ranking_overlap,
        },
        "reference_frame_transform": {"ligand_heavy_atoms": 18, "ligand_fit_rmsd_A": tfm[3]},
        "sstmap_alignment_sensitivity": {
            "local_sites": len(local["index"]),
            "global_sites": len(global_["index"]),
            "matched_pairs_within_1A": len(lg_pairs),
            "local_sites_unmatched": len(lg_unmatched),
            "center_shift_A": {
                "mean": float(lg_dist.mean()),
                "median": float(np.median(lg_dist)),
                "max": float(lg_dist.max()),
            },
            "occupancy": comparison_stats(
                local["occupancy"][local_i], global_["occupancy"][global_i]
            ),
            "minus_T_delta_S_total_kcal_mol": comparison_stats(
                -local["TStot"][local_i] * temperature_scale,
                -global_["TStot"][global_i] * temperature_scale,
            ),
        },
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (HYDRA_JSON, HYDRA_OBS, SST_LOCAL, SST_GLOBAL, SST_LIGAND)
        },
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    # Figures.
    plt.style.use("seaborn-v0_8-whitegrid")
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(*hxyz.T, s=55, label="HydraRank", color="#2266aa")
    ax.scatter(*lxyz.T, s=45, label="SSTMap", marker="^", color="#dd7722")
    for hi, si, _ in pairs:
        ax.plot(*np.vstack([hxyz[hi], lxyz[si]]).T, color="0.65", lw=0.8)
    ax.set(xlabel="x (Å)", ylabel="y (Å)", zlabel="z (Å)", title="3RLP hydration-site centers")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "site_overlay.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5.5))
    scatter(ax, s_occ, h_occ, "SSTMap occupancy", "HydraRank occupancy", "Matched-site occupancy")
    fig.tight_layout()
    fig.savefig(FIG / "occupancy_comparison.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    for ax, x, y, title in zip(
        axes,
        (s_trans, s_orient, s_total),
        (h_trans, h_orient, h_total),
        ("Translational", "Orientational", "Total"),
        strict=True,
    ):
        scatter(ax, x, y, "SSTMap −TΔS (kcal/mol)", "HydraRank −TΔS (kcal/mol)", title)
    fig.tight_layout()
    fig.savefig(FIG / "entropy_comparison.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5.5))
    scatter(
        ax,
        s_hbond,
        h_hbond,
        "SSTMap solute H-bonds",
        "HydraRank solute H-bonds",
        "Native H-bond definitions",
    )
    fig.tight_layout()
    fig.savefig(FIG / "hbond_comparison.png", dpi=220)
    plt.close(fig)

    scores = hydra_scores
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    for ax, y, title in zip(
        axes,
        (sst_score_proxy, local["Esw"][lidx], local["Etot"][lidx]),
        ("Derived SSTMap score proxy", "SSTMap solute–water energy", "SSTMap total energy"),
        strict=True,
    ):
        ax.scatter(scores, y, color="#5b4b9a", s=45)
        ax.set(xlabel="HydraRank displacement score", ylabel="kcal/mol", title=title)
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIG / "ranking_context.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    scatter(
        axes[0],
        global_["occupancy"][global_i],
        local["occupancy"][local_i],
        "Global-fit occupancy",
        "Local-fit occupancy",
        "SSTMap occupancy",
    )
    scatter(
        axes[1],
        -global_["TStot"][global_i] * temperature_scale,
        -local["TStot"][local_i] * temperature_scale,
        "Global-fit −TΔS",
        "Local-fit −TΔS",
        "SSTMap entropy (kcal/mol)",
    )
    axes[2].hist(lg_dist, bins=8, color="#338866", edgecolor="white")
    axes[2].set(
        xlabel="Center shift (Å)", ylabel="Matched sites", title="Alignment-induced site shift"
    )
    fig.tight_layout()
    fig.savefig(FIG / "alignment_sensitivity.png", dpi=220)
    plt.close(fig)

    m = metrics["matched_site_metrics"]
    sd = metrics["site_detection"]
    sens = metrics["sstmap_alignment_sensitivity"]
    rank_metric = m["hydrarank_score_vs_sstmap_derived_score_proxy"]
    report = f"""# HydraRank–SSTMap comparison: HSP90–3RLP

## Executive result

Under harmonized preprocessing, HydraRank and SSTMap identify nearly the same core hydration-site geometry: all **{len(hydra)} HydraRank sites** have unique SSTMap matches within 1.0 Å, while SSTMap identifies **{len(local["index"])} sites**, including two additional lower-occupancy sites. Matched-center separation is {sd["center_distance_A"]["mean"]:.3f} Å on average and {sd["center_distance_A"]["max"]:.3f} Å at maximum.

This is strong evidence of implementation-level spatial concordance for 3RLP. It is not yet scientific validation: it covers one system and the two methods share related inhomogeneous-solvation/nearest-neighbor ideas.

## Quantitative comparison

| Quantity (19 matched sites) | Pearson r | Spearman ρ | HydraRank − SSTMap mean | MAE | RMSE |
|---|---:|---:|---:|---:|---:|
| Occupancy | {m["occupancy"]["pearson_r"]:.3f} | {m["occupancy"]["spearman_rho"]:.3f} | {m["occupancy"]["mean_hydrarank_minus_sstmap"]:.3f} | {m["occupancy"]["mae"]:.3f} | {m["occupancy"]["rmse"]:.3f} |
| −TΔS translational (kcal/mol) | {m["minus_T_delta_S_trans_kcal_mol"]["pearson_r"]:.3f} | {m["minus_T_delta_S_trans_kcal_mol"]["spearman_rho"]:.3f} | {m["minus_T_delta_S_trans_kcal_mol"]["mean_hydrarank_minus_sstmap"]:.3f} | {m["minus_T_delta_S_trans_kcal_mol"]["mae"]:.3f} | {m["minus_T_delta_S_trans_kcal_mol"]["rmse"]:.3f} |
| −TΔS orientational (kcal/mol) | {m["minus_T_delta_S_orient_kcal_mol"]["pearson_r"]:.3f} | {m["minus_T_delta_S_orient_kcal_mol"]["spearman_rho"]:.3f} | {m["minus_T_delta_S_orient_kcal_mol"]["mean_hydrarank_minus_sstmap"]:.3f} | {m["minus_T_delta_S_orient_kcal_mol"]["mae"]:.3f} | {m["minus_T_delta_S_orient_kcal_mol"]["rmse"]:.3f} |
| −TΔS total (kcal/mol) | {m["minus_T_delta_S_total_kcal_mol"]["pearson_r"]:.3f} | {m["minus_T_delta_S_total_kcal_mol"]["spearman_rho"]:.3f} | {m["minus_T_delta_S_total_kcal_mol"]["mean_hydrarank_minus_sstmap"]:.3f} | {m["minus_T_delta_S_total_kcal_mol"]["mae"]:.3f} | {m["minus_T_delta_S_total_kcal_mol"]["rmse"]:.3f} |
| Solute H-bonds (native definitions) | {m["solute_hbonds_native_definitions"]["pearson_r"]:.3f} | {m["solute_hbonds_native_definitions"]["spearman_rho"]:.3f} | {m["solute_hbonds_native_definitions"]["mean_hydrarank_minus_sstmap"]:.3f} | {m["solute_hbonds_native_definitions"]["mae"]:.3f} | {m["solute_hbonds_native_definitions"]["rmse"]:.3f} |

Positive entropy values above are the unfavorable cost **−TΔS**. HydraRank's dimensionless entropies were converted with −RTΔS/R at 303.15 K. SSTMap's reported 300 K TΔS values were sign-inverted and multiplied by 303.15/300; this is an exact temperature normalization for SSTMap's linear T·S conversion and does not require rerunning its trajectory analysis. The bootstrap intervals in `metrics.json` describe across-site variation only and are **not trajectory-sampling confidence intervals**.

## What can and cannot be compared

- Occupancy, site positions, and translational/orientational/total entropy are direct matched-site comparisons after unit/sign normalization.
- Solute hydrogen bonds are informative but not estimator-equivalent: HydraRank uses a >130° donor–H–acceptor threshold, while SSTMap uses the stricter equivalent of >150°.
- SSTMap provides water–solute, water–water, and total interaction energies; HydraRank has no energy analogue.
- HydraRank provides residence time, enclosure, category, and a displacement ranking; SSTMap has no direct analogues. Correlations between the HydraRank score and SSTMap energies are therefore exploratory, not a validation of the rank.
- SSTMap's `f_enc` is undefined in this HSA implementation, so it cannot validate HydraRank enclosure.

## Ranking comparison

SSTMap does not provide a displacement ranking, so a direct native-rank comparison is impossible. As a diagnostic, applying HydraRank's score form to SSTMap outputs—SSTMap −TΔS minus 1 kcal/mol per SSTMap solute H-bond—gives Pearson r={rank_metric["pearson_r"]:.3f} and Spearman ρ={rank_metric["spearman_rho"]:.3f} against the HydraRank score. The top-3, top-5, and top-10 overlaps are {ranking_overlap["top_3"]["overlap"]}/3, {ranking_overlap["top_5"]["overlap"]}/5, and {ranking_overlap["top_10"]["overlap"]}/10. This is meaningful moderate agreement, but the proxy is not an independent SSTMap ranking and its H-bond definition differs.

## Alignment sensitivity

The locally aligned definitive SSTMap calculation gives {len(local["index"])} sites; the earlier whole-protein-aligned control gives {len(global_["index"])}. There are {len(lg_pairs)} matched sites within 1 Å. Their center shift is {sens["center_shift_A"]["mean"]:.3f} Å on average (maximum {sens["center_shift_A"]["max"]:.3f} Å); total entropy has Pearson r={sens["minus_T_delta_S_total_kcal_mol"]["pearson_r"]:.3f} between alignments. This confirms that local-pocket alignment is not a cosmetic choice and that only the local-fit run should be used for the primary controlled comparison.

## Interpretation and next validation steps

For 3RLP, the tools agree exceptionally well on where the principal waters are, but quantitative thermodynamic agreement must be judged from the statistics above rather than inferred from spatial overlap. Any systematic entropy offset is scientifically plausible because the tools do not have identical estimators and implementation details.

The next defensible validation step is to repeat the controlled analysis for 3RLQ and 3RLR, then measure convergence by trajectory blocks or independent replicas. Experimental displacement data or ligand-series affinity changes would be needed to test whether HydraRank's product-level ranking is predictive; SSTMap alone is a computational comparator, not ground truth.

## Reproducibility

Both methods used 25,001 frames at 2 ps spacing, the same 61-atom local-pocket fit, a 5 Å hydration region, 1 Å sites, 0.0333 Å⁻³ bulk density, and a 2× density threshold. Entropy costs are reported at the simulation temperature of 303.15 K. Exact input checksums are recorded in `metrics.json`. Recreate every table, metric, and figure with:

```bash
uv run python comparison/3rlp/compare.py
```
"""
    (OUT / "report.md").write_text(report)
    print(
        json.dumps(
            {
                "matched_sites": len(pairs),
                "sstmap_extra_sites": len(unmatched_l),
                "mean_center_distance_A": distances.mean(),
                "total_entropy_pearson_r": m["minus_T_delta_S_total_kcal_mol"]["pearson_r"],
                "occupancy_pearson_r": m["occupancy"]["pearson_r"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
