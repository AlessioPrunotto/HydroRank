"""Check shared atom identities, coordinate preparation, and pilot outputs."""

import argparse
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from common import ROOT, dump
from scipy.spatial.distance import cdist

from hydrarank.selections import analyse_water_topology, select_water


def validate(campaign):
    raw = mda.Universe(str(ROOT / "sample_traj/3rlp/step5.tpr"))
    sst = mda.Universe(str(ROOT / "sstmap_analysis/3rlp/input/3rlp_sstmap_topology.gro"))
    assert len(raw.atoms) == len(sst.atoms) == 62822
    for attr in ("names", "resnames", "resindices"):
        assert np.array_equal(getattr(raw.atoms, attr), getattr(sst.atoms, attr)), attr
    water = analyse_water_topology(select_water(raw))
    assert np.array_equal(water.oxygen_ix, sst.select_atoms("resname OPC and name O").indices)
    assert len(water.oxygen_ix) == 14799
    generated = mda.Universe(
        str(ROOT / "sample_traj/3rlp/step5.tpr"), str(campaign / "preparation/aligned.xtc")
    )
    ligand = mda.Universe(str(campaign / "preparation/ligand.pdb"))
    generated.trajectory[0]
    assert (
        np.max(
            np.abs(
                generated.select_atoms("resname 3RP and mass 2:100").positions
                - ligand.atoms.positions
            )
        )
        < 0.02
    )
    checks = dict(
        atom_order="exact names/resnames/resindices match",
        atom_count=62822,
        opc_oxygen_indices="exact match",
        opc_oxygen_count=14799,
        preparation=json.loads((campaign / "preparation/preparation.json").read_text()),
    )
    # Independent reference membership remains available even if GIST crashes.
    spacing = 0.5
    lower = np.floor((ligand.atoms.positions.min(axis=0) - 5) / spacing) * spacing
    upper = np.ceil((ligand.atoms.positions.max(axis=0) + 5) / spacing) * spacing
    oxygens = generated.atoms[water.oxygen_ix]
    counts = []
    for _ in generated.trajectory:
        xyz = oxygens.positions
        counts.append(int(np.count_nonzero(np.all((xyz >= lower) & (xyz < upper), axis=1))))
    checks["gist_reference_box"] = dict(
        lower_A=lower.tolist(),
        upper_A=upper.tolist(),
        spacing_A=spacing,
        dimensions=np.rint((upper - lower) / spacing).astype(int).tolist(),
        volume_A3=float(np.prod(upper - lower)),
        frames=len(counts),
        reference_mean_waters=float(np.mean(counts)),
        reference_min_waters=min(counts),
        reference_max_waters=max(counts),
        boundary="lower inclusive, upper exclusive; independent Cartesian count, not GIST output",
    )
    if (campaign / "preflight_hsa_full/adapter.json").exists():
        full = np.loadtxt(
            campaign / "preflight_hsa_full/benchmark_hsa_summary.txt",
            skiprows=1,
            usecols=range(27),
            ndmin=2,
        )
        reduced = np.loadtxt(
            campaign / "preflight_hsa_reduced/benchmark_hsa_summary.txt",
            skiprows=1,
            usecols=range(27),
            ndmin=2,
        )
        np.testing.assert_allclose(full[:, :6], reduced[:, :6], atol=1e-5)
        np.testing.assert_allclose(full[:, 14:17], reduced[:, 14:17], atol=1e-5)
        assert np.any(np.abs(full[:, 6:13]) > 0)
        assert np.allclose(reduced[:, 6:13], 0, equal_nan=False)
        hydra = json.loads((campaign / "preflight_hydrarank/results/results.json").read_text())
        hrows = hydra["ranking"]["sites"]
        positions = np.array([[r["x"], r["y"], r["z"]] for r in hrows])
        checks["preflight_hsa_full_and_reduced"] = (
            "identical site positions, occupancies and entropy; reduced energies zero"
        )
        checks["hydra_to_hsa_nearest_center_A"] = (
            cdist(positions, full[:, 1:4]).min(axis=1).tolist()
        )
        checks["hydra_short_site_count"] = len(hrows)
        checks["sst_short_site_count"] = len(full)
    cached = list(campaign.glob("cached_hydrarank_*/adapter.json"))
    checks["cache_reuse"] = {
        str(p.parent.name): json.loads(p.read_text())["result"]["cache_reused"] for p in cached
    }
    assert all(checks["cache_reuse"].values())
    dump(campaign / "validation.json", checks)
    return checks


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("campaign", type=Path)
    args = p.parse_args()
    validate(args.campaign.resolve())
    print("Validation passed")
