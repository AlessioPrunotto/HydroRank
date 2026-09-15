"""Clustering of water oxygen positions into hydration sites.

Uses the density-peak scheme that WaterMap/SSTMap popularised: repeatedly take the
position with the most neighbours inside a small sphere, call it a site, remove the
waters it claims, and continue until no remaining peak is denser than bulk water.
It needs nothing beyond a KD-tree, is deterministic, and produces sites of a fixed,
physically meaningful radius -- unlike k-means (which needs the number of sites up
front) or DBSCAN (whose clusters can grow into elongated, unphysical blobs).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from hydrarank.data import WaterObservations
from hydrarank.entropy import BULK_WATER_DENSITY
from hydrarank.exceptions import HydraRankError

#: Radius of a hydration site: roughly half the O-O distance of two hydrogen-bonded
#: waters, so that two sites cannot describe the same water.
DEFAULT_SITE_RADIUS = 1.0


@dataclass
class HydrationSites:
    """Hydration sites and the assignment of every water observation to them."""

    centers: np.ndarray  # (n_sites, 3)
    labels: np.ndarray  # (n_obs,) site index, -1 when unassigned
    radius: float
    n_frames: int
    min_count: int

    def __post_init__(self) -> None:
        if self.centers.ndim != 2 or self.centers.shape[1] != 3:
            raise HydraRankError(f"centers must have shape (n_sites, 3), got {self.centers.shape}")
        if self.labels.ndim != 1:
            raise HydraRankError("labels must be one-dimensional")
        if not np.issubdtype(self.labels.dtype, np.integer):
            raise HydraRankError("labels must contain integers")
        if self.radius <= 0:
            raise HydraRankError(f"radius must be > 0, got {self.radius}")
        if self.n_frames < 1:
            raise HydraRankError(f"n_frames must be >= 1, got {self.n_frames}")
        if self.min_count < 1:
            raise HydraRankError(f"min_count must be >= 1, got {self.min_count}")
        if not np.all(np.isfinite(self.centers)):
            raise HydraRankError("centers contain non-finite values")
        if self.labels.size and (np.min(self.labels) < -1 or np.max(self.labels) >= self.n_sites):
            raise HydraRankError("labels contain an invalid site index")

    @property
    def n_sites(self) -> int:
        return int(self.centers.shape[0])

    @property
    def counts(self) -> np.ndarray:
        """Number of water observations claimed by each site."""
        return np.bincount(self.labels[self.labels >= 0], minlength=self.n_sites)

    def frame_occupancy(self, frame: np.ndarray) -> np.ndarray:
        """Fraction of frames in which each site holds at least one water."""
        occupancy = np.zeros(self.n_sites)
        for site in range(self.n_sites):
            occupied = np.unique(frame[self.labels == site])
            occupancy[site] = occupied.size / self.n_frames
        return occupancy

    def mean_waters(self) -> np.ndarray:
        """Average number of waters per frame in each site; can exceed 1."""
        return self.counts / self.n_frames

    def spread(self, positions: np.ndarray) -> np.ndarray:
        """RMS distance of each site's members from its centre."""
        values = np.zeros(self.n_sites)
        for site in range(self.n_sites):
            members = positions[self.labels == site]
            values[site] = np.sqrt(np.mean(np.sum((members - self.centers[site]) ** 2, axis=1)))
        return values

    def summary(self, observations: WaterObservations) -> list[dict[str, Any]]:
        occupancy = self.frame_occupancy(observations.frame)
        mean_waters = self.mean_waters()
        spread = self.spread(observations.oxygen)
        counts = self.counts
        return [
            {
                "site": index + 1,
                "x": float(self.centers[index, 0]),
                "y": float(self.centers[index, 1]),
                "z": float(self.centers[index, 2]),
                "n_waters": int(counts[index]),
                "occupancy": float(occupancy[index]),
                "mean_waters": float(mean_waters[index]),
                "spread": float(spread[index]),
            }
            for index in range(self.n_sites)
        ]

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_sites": self.n_sites,
            "radius": self.radius,
            "n_frames": self.n_frames,
            "min_count": self.min_count,
            "centers": self.centers.tolist(),
            "counts": self.counts.tolist(),
        }


def bulk_equivalent_count(n_frames: int, radius: float, density_factor: float = 2.0) -> int:
    """Occupancy a sphere of ``radius`` would reach at ``density_factor`` times bulk density.

    Used as the stopping criterion: a site that is no denser than bulk water carries no
    information about the pocket.
    """
    if n_frames < 1:
        raise HydraRankError(f"n_frames must be >= 1, got {n_frames}")
    if radius <= 0:
        raise HydraRankError(f"radius must be > 0, got {radius}")
    if density_factor <= 0:
        raise HydraRankError(f"density_factor must be > 0, got {density_factor}")
    volume = 4.0 / 3.0 * np.pi * radius**3
    return max(1, int(np.ceil(density_factor * BULK_WATER_DENSITY * volume * n_frames)))


def cluster_positions(
    positions: np.ndarray,
    radius: float = DEFAULT_SITE_RADIUS,
    min_count: int = 1,
    max_sites: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Density-peak clustering of a point cloud.

    Returns the site centres (refined to the centroid of their members) and the
    per-point site labels, with ``-1`` for points that belong to no site.
    """
    positions = np.asarray(positions, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise HydraRankError(f"expected an (n, 3) array of positions, got {positions.shape}")
    if radius <= 0:
        raise HydraRankError(f"radius must be > 0, got {radius}")
    if min_count < 1:
        raise HydraRankError(f"min_count must be >= 1, got {min_count}")
    if max_sites is not None and max_sites < 1:
        raise HydraRankError(f"max_sites must be >= 1, got {max_sites}")

    n_points = positions.shape[0]
    labels = np.full(n_points, -1, dtype=np.int64)
    if n_points == 0:
        return np.empty((0, 3)), labels

    neighbours = cKDTree(positions).query_ball_point(positions, r=radius)
    counts = np.fromiter((len(item) for item in neighbours), dtype=np.int64, count=n_points)
    available = np.ones(n_points, dtype=bool)

    centers = []
    while True:
        if max_sites is not None and len(centers) >= max_sites:
            break
        peak = int(np.argmax(np.where(available, counts, -1)))
        if counts[peak] < min_count or not available[peak]:
            break

        members = np.array([i for i in neighbours[peak] if available[i]], dtype=np.int64)
        labels[members] = len(centers)
        centers.append(positions[members].mean(axis=0))

        available[members] = False
        for member in members:
            for neighbour in neighbours[member]:
                counts[neighbour] -= 1

    return np.asarray(centers).reshape(-1, 3), labels


def cluster_hydration_sites(
    observations: WaterObservations,
    radius: float = DEFAULT_SITE_RADIUS,
    density_factor: float = 2.0,
    min_count: int | None = None,
    max_sites: int | None = None,
) -> HydrationSites:
    """Cluster the oxygen positions of ``observations`` into hydration sites."""
    threshold = (
        min_count
        if min_count is not None
        else bulk_equivalent_count(observations.n_frames, radius, density_factor)
    )
    centers, labels = cluster_positions(
        observations.oxygen, radius=radius, min_count=threshold, max_sites=max_sites
    )
    return HydrationSites(
        centers=centers,
        labels=labels,
        radius=radius,
        n_frames=observations.n_frames,
        min_count=threshold,
    )
