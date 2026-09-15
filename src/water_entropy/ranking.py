"""Ranking of hydration sites by how much there is to gain from displacing them.

Two quantities decide the answer and they pull in opposite directions:

* the entropy released when the water is freed, ``-T dS``, which is large for sites
  that hold a water in a fixed position and orientation;
* the hydrogen bonds that water makes to protein and ligand, which the new ligand
  substituent has to replace, and which are pure cost if it cannot.

The score below simply subtracts the second from the first. The conversion factor
``hbond_penalty`` is a blunt instrument -- real hydrogen-bond strengths vary with
geometry and environment -- so the category labels matter more than the number: they
say *why* a site is or is not a target, which is what a medicinal chemist acts on.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from water_entropy.analysis import SiteAnalysis
from water_entropy.exceptions import WaterEntropyError

#: Ordered and weakly bound: displacing it should pay off directly.
DISPLACEABLE = "displaceable"

#: Ordered but hydrogen bonded to the solute: worth targeting only if the ligand can
#: reproduce those hydrogen bonds.
REPLACE_HBONDS = "replace-hbonds"

#: Already as disordered as bulk water, or too rarely occupied to tell.
BULK_LIKE = "bulk-like"

DEFAULT_ENTROPY_THRESHOLD = 2.0
DEFAULT_HBOND_THRESHOLD = 1.0


@dataclass
class SiteRanking:
    """Hydration sites ordered by displacement score."""

    analysis: SiteAnalysis
    score: np.ndarray
    category: list[str]
    order: np.ndarray
    hbond_penalty: float
    entropy_threshold: float
    hbond_threshold: float

    def rows(self) -> list[dict[str, Any]]:
        source = self.analysis.rows()
        return [
            {
                **source[site],
                "rank": position + 1,
                "score": float(self.score[site]),
                "category": self.category[site],
            }
            for position, site in enumerate(self.order)
        ]

    def to_dict(self, top: int | None = None) -> dict[str, Any]:
        rows = self.rows()
        return {
            "hbond_penalty": self.hbond_penalty,
            "entropy_threshold": self.entropy_threshold,
            "hbond_threshold": self.hbond_threshold,
            "temperature": self.analysis.temperature,
            "sites": rows if top is None else rows[:top],
        }


def rank_sites(
    analysis: SiteAnalysis,
    hbond_penalty: float = 1.0,
    entropy_threshold: float = DEFAULT_ENTROPY_THRESHOLD,
    hbond_threshold: float = DEFAULT_HBOND_THRESHOLD,
) -> SiteRanking:
    """Score and classify every hydration site as a displacement target."""
    if hbond_penalty < 0:
        raise WaterEntropyError(f"hbond_penalty must be >= 0, got {hbond_penalty}")
    if entropy_threshold < 0:
        raise WaterEntropyError(f"entropy_threshold must be >= 0, got {entropy_threshold}")
    if hbond_threshold < 0:
        raise WaterEntropyError(f"hbond_threshold must be >= 0, got {hbond_threshold}")
    entropy_gain = analysis.minus_t_delta_s
    score = entropy_gain - hbond_penalty * analysis.hb_solute

    category = []
    for index in range(analysis.n_sites):
        if not np.isfinite(entropy_gain[index]) or entropy_gain[index] < entropy_threshold:
            category.append(BULK_LIKE)
        elif analysis.hb_solute[index] >= hbond_threshold:
            category.append(REPLACE_HBONDS)
        else:
            category.append(DISPLACEABLE)

    order = np.argsort(np.where(np.isfinite(score), -score, np.inf), kind="stable")
    return SiteRanking(
        analysis=analysis,
        score=score,
        category=category,
        order=order,
        hbond_penalty=hbond_penalty,
        entropy_threshold=entropy_threshold,
        hbond_threshold=hbond_threshold,
    )


def format_ranking(ranking: SiteRanking, top: int | None = None) -> str:
    """Render the ranked table of displacement targets."""
    rows = ranking.rows()
    shown = rows if top is None else rows[:top]

    header = [
        f"Displacement ranking ({len(rows)} sites, T = {ranking.analysis.temperature:g} K, "
        f"{ranking.hbond_penalty:g} kcal/mol per solute hydrogen bond)",
        "",
        f"{'rank':>4} {'site':>4} {'occup':>6} {'encl':>5} {'hb':>5} "
        f"{'-TdS':>6} {'score':>6}  category",
    ]
    body = [
        f"{row['rank']:>4} {row['site']:>4} {row['occupancy']:>6.2f} {row['enclosure']:>5.1f} "
        f"{row['hb_solute']:>5.2f} {row['minus_t_delta_s']:>6.2f} {row['score']:>6.2f}  "
        f"{row['category']}"
        for row in shown
    ]
    counts = {name: ranking.category.count(name) for name in (DISPLACEABLE, REPLACE_HBONDS)}
    footer = [
        "",
        "encl: solute heavy atoms around the water (how buried it is).",
        "hb:   mean hydrogen bonds to protein and ligand per frame.",
        "score = -TdS - penalty * hb; the penalty is a heuristic, so treat the",
        "category as the primary signal and the score only as a tie-breaker.",
        "",
        f"{counts[DISPLACEABLE]} sites are ordered and weakly bound ({DISPLACEABLE}); "
        f"{counts[REPLACE_HBONDS]} are ordered but hydrogen bonded ({REPLACE_HBONDS}).",
    ]
    return "\n".join([*header, *body, *footer])
