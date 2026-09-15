from types import SimpleNamespace

import numpy as np
import pytest

from hydrarank.exceptions import HydraRankError
from hydrarank.ranking import (
    BULK_LIKE,
    DISPLACEABLE,
    REPLACE_HBONDS,
    format_ranking,
    rank_sites,
)


def _analysis():
    rows = [
        {"site": 1, "occupancy": 0.9, "enclosure": 8.0, "hb_solute": 0.2, "minus_t_delta_s": 4.0},
        {"site": 2, "occupancy": 0.8, "enclosure": 6.0, "hb_solute": 1.5, "minus_t_delta_s": 4.5},
        {"site": 3, "occupancy": 0.4, "enclosure": 2.0, "hb_solute": 0.0, "minus_t_delta_s": 1.0},
    ]
    return SimpleNamespace(
        n_sites=3,
        temperature=300.0,
        minus_t_delta_s=np.array([4.0, 4.5, 1.0]),
        hb_solute=np.array([0.2, 1.5, 0.0]),
        rows=lambda: rows,
    )


def test_rank_sites_scores_orders_and_classifies_sites():
    ranking = rank_sites(_analysis(), hbond_penalty=1.0)

    assert ranking.score == pytest.approx([3.8, 3.0, 1.0])
    assert list(ranking.order) == [0, 1, 2]
    assert ranking.category == [DISPLACEABLE, REPLACE_HBONDS, BULK_LIKE]


def test_ranking_output_can_be_limited():
    ranking = rank_sites(_analysis())

    assert len(ranking.to_dict(top=2)["sites"]) == 2
    table = format_ranking(ranking, top=1)
    assert "Displacement ranking" in table
    rows = [line.split() for line in table.splitlines()]
    assert sum(len(row) == 8 and row[0].isdigit() and row[1].isdigit() for row in rows) == 1


def test_ties_are_stable_and_nan_scores_sort_last():
    analysis = _analysis()
    analysis.minus_t_delta_s = np.array([3.0, 3.0, np.nan])
    analysis.hb_solute = np.zeros(3)

    ranking = rank_sites(analysis)

    assert list(ranking.order) == [0, 1, 2]
    assert ranking.category == [DISPLACEABLE, DISPLACEABLE, BULK_LIKE]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"hbond_penalty": -1.0},
        {"entropy_threshold": -1.0},
        {"hbond_threshold": -1.0},
    ],
)
def test_ranking_rejects_negative_parameters(kwargs):
    with pytest.raises(HydraRankError, match="must be >= 0"):
        rank_sites(_analysis(), **kwargs)
