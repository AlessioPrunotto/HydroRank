import csv

import pytest

from hydrarank.exceptions import HydraRankError
from hydrarank.export import write_csv, write_site_coordinates

ROWS = [
    {
        "site": 1,
        "x": 1.25,
        "y": 2.5,
        "z": 3.75,
        "occupancy": 0.8,
        "minus_t_delta_s": 2.4,
    }
]


def test_write_csv(tmp_path):
    path = write_csv(ROWS, tmp_path / "nested" / "sites.csv")
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["site"] == "1"
    assert rows[0]["occupancy"] == "0.8"


@pytest.mark.parametrize("suffix,marker", [(".pdb", "HETATM"), (".cif", "_atom_site.Cartn_x")])
def test_write_site_coordinates(tmp_path, suffix, marker):
    path = write_site_coordinates(ROWS, tmp_path / f"sites{suffix}")
    assert marker in path.read_text()


def test_export_rejects_empty_rows_and_unknown_coordinate_format(tmp_path):
    with pytest.raises(HydraRankError, match="no hydration sites"):
        write_csv([], tmp_path / "empty.csv")
    with pytest.raises(HydraRankError, match="must end"):
        write_site_coordinates(ROWS, tmp_path / "sites.xyz")
