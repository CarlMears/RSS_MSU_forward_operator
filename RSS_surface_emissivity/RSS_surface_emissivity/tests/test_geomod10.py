from pathlib import Path
from typing import Final

import pytest
import numpy as np


def test_sanity() -> None:
    import geomod10

    assert hasattr(geomod10, "wind_emiss")


# The Fortran implementation has paths under this prefix hard-coded
O_DRIVE: Final = Path("/mnt/oserver/o")


# Since the Fortran code exits forcefully if any of the input files don't exist,
# then don't even attempt to run this, just flag it in the pytest output
@pytest.mark.xfail(not O_DRIVE.exists(), reason="O: drive is unreachable", run=False)
def test_values() -> None:
    from geomod10 import wind_emiss

    freq = 19.0
    tht = 53.0
    sst = np.array([10.0, 10.0, 10.0, 10.0], dtype=np.float32)
    wind = np.array([5.0, 7.0, 9.0, 11.0], dtype=np.float32)
    phir = np.array([-999.0, -999.0, -999.0, -999.0], dtype=np.float32)

    emiss = wind_emiss(freq, tht, sst, wind, phir)

    expected_emiss = np.array(
        [
            [0.58449423, 0.28792182],
            [0.5853587, 0.29340607],
            [0.588001, 0.30069447],
            [0.59222317, 0.31007487],
        ],
        dtype=np.float32,
    )
    assert np.array_equal(emiss, expected_emiss)
