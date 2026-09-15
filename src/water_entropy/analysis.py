"""The hydration-site table: occupancy, persistence and entropy, site by site."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from water_entropy.clustering import HydrationSites
from water_entropy.data import WaterObservations
from water_entropy.entropy import (
    MIN_SAMPLES,
    minus_t_delta_s,
    orientational_entropy,
    translational_entropy,
    water_orientations,
)
from water_entropy.exceptions import WaterEntropyError


@dataclass(frozen=True)
class ResidenceStats:
    """How long individual water molecules stay in a site before being replaced."""

    mean_frames: float
    max_frames: int
    n_episodes: int
    n_censored: int

    def mean_ps(self, dt_ps: float) -> float:
        return self.mean_frames * dt_ps


@dataclass
class SiteAnalysis:
    """Per-site metrics, one row per hydration site."""

    sites: HydrationSites
    observations: WaterObservations
    occupancy: np.ndarray
    mean_waters: np.ndarray
    spread: np.ndarray
    residence: list[ResidenceStats]
    s_trans: np.ndarray
    s_orient: np.ndarray
    hb_solute: np.ndarray
    enclosure: np.ndarray
    temperature: float

    @property
    def n_sites(self) -> int:
        return self.sites.n_sites

    @property
    def s_total(self) -> np.ndarray:
        return self.s_trans + self.s_orient

    @property
    def minus_t_delta_s(self) -> np.ndarray:
        """Free-energy cost of the ordering, in kcal/mol; larger means more to gain."""
        return minus_t_delta_s(self.s_total, self.temperature)

    def rows(self) -> list[dict[str, Any]]:
        dt = self.observations.dt_ps
        return [
            {
                "site": index + 1,
                "x": float(self.sites.centers[index, 0]),
                "y": float(self.sites.centers[index, 1]),
                "z": float(self.sites.centers[index, 2]),
                "n_waters": int(self.sites.counts[index]),
                "occupancy": float(self.occupancy[index]),
                "mean_waters": float(self.mean_waters[index]),
                "spread": float(self.spread[index]),
                "mean_residence_ps": self.residence[index].mean_ps(dt),
                "n_episodes": self.residence[index].n_episodes,
                "hb_solute": float(self.hb_solute[index]),
                "enclosure": float(self.enclosure[index]),
                "s_trans": float(self.s_trans[index]),
                "s_orient": float(self.s_orient[index]),
                "s_total": float(self.s_total[index]),
                "minus_t_delta_s": float(self.minus_t_delta_s[index]),
            }
            for index in range(self.n_sites)
        ]

    def to_dict(self) -> dict[str, Any]:
        return {
            "temperature": self.temperature,
            "site_radius": self.sites.radius,
            "n_frames": self.sites.n_frames,
            "dt_ps": self.observations.dt_ps,
            "sites": self.rows(),
        }


def residence_stats(
    frames: np.ndarray, resids: np.ndarray, n_frames: int, max_gap: int = 0
) -> ResidenceStats:
    """Length of the uninterrupted visits of individual waters to one site.

    ``max_gap`` frames of absence are tolerated within a single visit, which avoids
    counting a water that briefly wobbles out of the site as a new arrival.
    """
    if frames.size == 0:
        return ResidenceStats(mean_frames=0.0, max_frames=0, n_episodes=0, n_censored=0)

    lengths, censored = [], 0
    for resid in np.unique(resids):
        visits = np.unique(frames[resids == resid])
        breaks = np.flatnonzero(np.diff(visits) > max_gap + 1) + 1
        for episode in np.split(visits, breaks):
            lengths.append(int(episode[-1] - episode[0] + 1))
            if episode[0] == 0 or episode[-1] == n_frames - 1:
                censored += 1

    return ResidenceStats(
        mean_frames=float(np.mean(lengths)),
        max_frames=int(np.max(lengths)),
        n_episodes=len(lengths),
        n_censored=censored,
    )


def analyse_sites(
    observations: WaterObservations,
    sites: HydrationSites,
    temperature: float = 300.0,
    max_gap: int = 0,
    min_samples: int = MIN_SAMPLES,
) -> SiteAnalysis:
    """Compute occupancy, persistence and entropy proxies for every hydration site."""
    if sites.labels.shape != (observations.n_observations,):
        raise WaterEntropyError(
            "site labels must contain one entry per water observation "
            f"({sites.labels.size} != {observations.n_observations})"
        )
    if sites.n_frames != observations.n_frames:
        raise WaterEntropyError(
            f"sites and observations disagree on frame count ({sites.n_frames} != "
            f"{observations.n_frames})"
        )
    orientations = water_orientations(observations.oxygen, observations.hydrogen)

    residence, s_trans, s_orient = [], [], []
    for index in range(sites.n_sites):
        members = np.flatnonzero(sites.labels == index)
        residence.append(
            residence_stats(
                observations.frame[members],
                observations.water_id[members],
                sites.n_frames,
                max_gap=max_gap,
            )
        )
        s_trans.append(translational_entropy(observations.oxygen[members], min_samples=min_samples))
        s_orient.append(orientational_entropy(orientations[members], min_samples=min_samples))

    return SiteAnalysis(
        sites=sites,
        observations=observations,
        occupancy=sites.frame_occupancy(observations.frame),
        mean_waters=sites.mean_waters(),
        spread=sites.spread(observations.oxygen),
        residence=residence,
        s_trans=np.asarray(s_trans, dtype=float),
        s_orient=np.asarray(s_orient, dtype=float),
        hb_solute=_site_means(observations.hb_solute, sites),
        enclosure=_site_means(observations.enclosure, sites),
        temperature=temperature,
    )


def _site_means(values: np.ndarray, sites: HydrationSites) -> np.ndarray:
    return np.asarray(
        [
            (
                float(np.mean(values[sites.labels == index]))
                if np.any(sites.labels == index)
                else 0.0
            )
            for index in range(sites.n_sites)
        ]
    )


def format_analysis(analysis: SiteAnalysis) -> str:
    """Render the hydration-site table."""
    observations = analysis.observations
    sites = analysis.sites
    assigned = int(np.count_nonzero(sites.labels >= 0))

    header = [
        f"Hydration sites (radius {sites.radius:g} A, {sites.n_frames} frames at "
        f"{observations.dt_ps:g} ps, T = {analysis.temperature:g} K)",
        f"  {sites.n_sites} sites, {assigned} of {observations.n_observations} "
        "water observations assigned",
        "",
        f"{'site':>5} {'x':>7} {'y':>7} {'z':>7} {'wat':>5} {'occup':>6} {'res/ps':>7} "
        f"{'encl':>5} {'hb':>5} {'S_tr':>7} {'S_or':>7} {'-TdS':>7}",
    ]
    rows = [
        f"{row['site']:>5} {row['x']:>7.2f} {row['y']:>7.2f} {row['z']:>7.2f} "
        f"{row['n_waters']:>5} {row['occupancy']:>6.2f} {row['mean_residence_ps']:>7.0f} "
        f"{row['enclosure']:>5.1f} {row['hb_solute']:>5.2f} "
        f"{row['s_trans']:>7.2f} {row['s_orient']:>7.2f} {row['minus_t_delta_s']:>7.2f}"
        for row in analysis.rows()
    ]
    footer = [
        "",
        "S_tr, S_or: excess translational / orientational entropy per water relative to",
        "bulk, in units of the gas constant (negative = more ordered than bulk).",
        "-TdS: free-energy cost of that ordering at the given temperature, in kcal/mol.",
        "encl: mean nearby solute heavy atoms; hb: mean solute hydrogen bonds per water.",
    ]
    if observations.dt_ps > 20:
        footer.append(
            f"! frames are {observations.dt_ps:g} ps apart, so residence times are "
            "unresolved and only bound from below"
        )
    return "\n".join([*header, *rows, *footer])
