import json

import numpy as np

from hydrarank.jsonio import dumps


def test_dumps_produces_strict_json_for_nonfinite_scientific_values():
    encoded = dumps({"finite": np.float64(1.5), "missing": np.nan, "infinite": np.inf})

    assert json.loads(encoded) == {"finite": 1.5, "missing": None, "infinite": None}
    assert "NaN" not in encoded and "Infinity" not in encoded
