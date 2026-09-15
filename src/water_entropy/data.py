"""The data contract between preprocessing and everything downstream.

One row per (frame, water molecule) pair, in the aligned binding-site frame.
Angstrom throughout.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from water_entropy.exceptions import WaterEntropyError

_FORMAT_VERSION = 3
_SUPPORTED_FORMAT_VERSIONS = (2, 3)


@dataclass
class WaterObservations:
    """Water molecules seen near the ligand, superposed onto a common frame."""

    oxygen: np.ndarray  # (n_obs, 3)
    hydrogen: np.ndarray  # (n_obs, 2, 3)
    frame: np.ndarray  # (n_obs,) index into `frames`
    resid: np.ndarray  # (n_obs,)
    frames: np.ndarray  # (n_frames,) source trajectory frame numbers
    times: np.ndarray  # (n_frames,) picoseconds
    fit_rmsd: np.ndarray  # (n_frames,)
    ligand_reference: np.ndarray  # (n_ligand_heavy, 3)
    metadata: dict[str, Any] = field(default_factory=dict)
    hb_solute: np.ndarray | None = None  # (n_obs,) hydrogen bonds to protein + ligand
    enclosure: np.ndarray | None = None  # (n_obs,) solute heavy atoms nearby
    water_id: np.ndarray | None = None  # (n_obs,) globally unique topology residue index

    def __post_init__(self) -> None:
        n_obs = self.oxygen.shape[0]
        if self.hb_solute is None:
            self.hb_solute = np.zeros(n_obs)
        if self.enclosure is None:
            self.enclosure = np.zeros(n_obs)
        if self.water_id is None:
            # Compatibility for programmatic callers and format-2 caches. New
            # preprocessing always supplies the globally unique topology resindex.
            self.water_id = self.resid.copy()
        if self.oxygen.shape != (n_obs, 3):
            raise WaterEntropyError(f"oxygen has shape {self.oxygen.shape}, expected {(n_obs, 3)}")
        for name, expected in (
            ("hydrogen", (n_obs, 2, 3)),
            ("frame", (n_obs,)),
            ("resid", (n_obs,)),
            ("hb_solute", (n_obs,)),
            ("enclosure", (n_obs,)),
            ("water_id", (n_obs,)),
        ):
            actual = getattr(self, name).shape
            if actual != expected:
                raise WaterEntropyError(f"{name} has shape {actual}, expected {expected}")
        if self.frames.shape != self.times.shape or self.frames.shape != self.fit_rmsd.shape:
            raise WaterEntropyError("frames, times and fit_rmsd must have the same length")
        if self.frames.ndim != 1:
            raise WaterEntropyError("frames, times and fit_rmsd must be one-dimensional")
        if self.n_frames < 1:
            raise WaterEntropyError("frames must contain at least one trajectory frame")
        for name in ("frame", "resid", "water_id", "frames"):
            if not np.issubdtype(getattr(self, name).dtype, np.integer):
                raise WaterEntropyError(f"{name} must contain integers")
        if self.ligand_reference.ndim != 2 or self.ligand_reference.shape[1] != 3:
            raise WaterEntropyError("ligand_reference must have shape (n_ligand_heavy_atoms, 3)")
        if n_obs and (np.min(self.frame) < 0 or np.max(self.frame) >= self.n_frames):
            raise WaterEntropyError("observation frame indices fall outside the frames array")
        for name in (
            "oxygen",
            "hydrogen",
            "times",
            "fit_rmsd",
            "ligand_reference",
            "hb_solute",
            "enclosure",
        ):
            if not np.all(np.isfinite(getattr(self, name))):
                raise WaterEntropyError(f"{name} contains non-finite values")
        if np.any(self.hb_solute < 0) or np.any(self.enclosure < 0):
            raise WaterEntropyError("hb_solute and enclosure counts cannot be negative")
        if self.frames.size > 1 and np.any(np.diff(self.frames) <= 0):
            raise WaterEntropyError("source trajectory frame numbers must be strictly increasing")
        if self.times.size > 1:
            intervals = np.diff(self.times)
            if np.any(intervals <= 0):
                raise WaterEntropyError("trajectory times must be strictly increasing")
            if not np.allclose(intervals, intervals[0], rtol=1e-5, atol=1e-8):
                raise WaterEntropyError(
                    "trajectory times are irregular; residence times require uniformly "
                    "sampled frames"
                )

    @property
    def n_observations(self) -> int:
        return int(self.oxygen.shape[0])

    @property
    def n_frames(self) -> int:
        return int(self.frames.size)

    @property
    def dipole(self) -> np.ndarray:
        """Unit vector from the oxygen towards the midpoint of the two hydrogens."""
        vectors = self.hydrogen.mean(axis=1) - self.oxygen
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        return vectors / np.where(norms > 0, norms, 1.0)

    @property
    def dt_ps(self) -> float:
        if self.times.size < 2:
            return 0.0
        return float(self.times[1] - self.times[0])

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path,
            oxygen=self.oxygen,
            hydrogen=self.hydrogen,
            frame=self.frame,
            resid=self.resid,
            frames=self.frames,
            times=self.times,
            fit_rmsd=self.fit_rmsd,
            ligand_reference=self.ligand_reference,
            hb_solute=self.hb_solute,
            enclosure=self.enclosure,
            water_id=self.water_id,
            metadata=json.dumps({**self.metadata, "format_version": _FORMAT_VERSION}),
        )
        return path

    @classmethod
    def load(cls, path: str | Path) -> WaterObservations:
        with np.load(Path(path), allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata"]))
            version = metadata.pop("format_version", None)
            if version not in _SUPPORTED_FORMAT_VERSIONS:
                raise WaterEntropyError(
                    f"unsupported observations format version {version}, expected one of "
                    f"{_SUPPORTED_FORMAT_VERSIONS}"
                )
            return cls(
                oxygen=data["oxygen"],
                hydrogen=data["hydrogen"],
                frame=data["frame"],
                resid=data["resid"],
                frames=data["frames"],
                times=data["times"],
                fit_rmsd=data["fit_rmsd"],
                ligand_reference=data["ligand_reference"],
                hb_solute=data["hb_solute"],
                enclosure=data["enclosure"],
                water_id=data["water_id"] if "water_id" in data.files else data["resid"],
                metadata=metadata,
            )
