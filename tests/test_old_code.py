import numpy as np
import pytest
from astropy.constants import M_sun, M_jup
from hypothesis import settings, given, strategies as st

import rv
from utils.old import RV_from_params, companion_amplitude
from utils.parse import parse_paramfile


@pytest.fixture(params=["tests/test_params.txt", 'tests/test_params2.txt'])
def params(request):
    """Load Parameter file."""
    return parse_paramfile(request.param)


def test_RV_from_params_circular(params):
    """Maximum RV should be within gamma + k1, gamma + k2 for a circular orbit."""
    params["eccentricity"] = 0  # Circular orbit
    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)
    rvs = RV_from_params(time, params)
    min_val = params["mean_val"] - params["k1"]
    max_val = params["mean_val"] + params["k1"]
    assert np.all(rvs <= max_val)
    assert np.all(rvs >= min_val)


@settings(deadline=400)  # double deadline for this test
@given(st.floats(min_value=0, max_value=1e6, allow_nan=False, allow_infinity=False),
       st.floats(min_value=0.1, max_value=1e5, allow_nan=False, allow_infinity=False),
       st.floats(min_value=0.1, max_value=1e6, allow_nan=False, allow_infinity=False),
       st.floats(min_value=0, max_value=360, allow_nan=False, allow_infinity=False),
       st.floats(min_value=0, max_value=0.999, allow_nan=False, allow_infinity=False))
def test_RV_from_params(k1, period, tau, omega, ecc):
    """RV should be within theoretical limits."""
    params = {"mean_val": 0.0, "k1": k1, "period": period, "tau": tau, "omega": omega, "eccentricity": ecc}

    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)
    rvs = RV_from_params(time, params) - params["mean_val"]  # remove center

    a_1 = params["k1"] * (1 + params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    b_1 = params["k1"] * (1 - params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    # Round to avoid floating point errors
    rvs = np.around(rvs, decimals=8)
    max_val = np.around(a_1, decimals=8)
    min_val = np.around(-b_1, decimals=8)

    max_rv = np.max(rvs)
    min_rv = np.min(rvs)

    assert max_val >= max_rv
    assert min_val <= min_rv
    assert np.all(rvs <= max_val)
    assert np.all(rvs >= min_val)
    assert np.allclose(params["k1"], 0.5 * (a_1 + b_1))


def test_from_params_companion(params):
    """Maximum RV should be within gamma+k1, gamma+k2 for a circular orbit."""
    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)
    rvs = RV_from_params(time, params, ignore_mean=True, companion=True)

    a_2 = params["k2"] * (1 + params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    b_2 = params["k2"] * (1 - params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))

    max_val = -b_2
    min_val = a_2

    max_rv = np.max(rvs)
    min_rv = np.min(rvs)
    print("max_val", max_val, "min_val", min_val)
    print("max_rv", max_rv, "min_rv", min_rv)
    assert np.all(rvs <= max_val)
    assert np.all(rvs >= min_val)
    if params["eccentricity"] == 0:
        assert np.allclose(a_2, b_2)
    else:
        assert not np.allclose(a_2, b_2)
    assert np.allclose(params["k2"], 0.5 * (a_2 + b_2))


def test_RV_ignore_mean(params):
    """Maximum RV should be within gamma+k1, gamma+k2 for a circular orbit."""
    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)

    assert np.allclose(RV_from_params(time, params, ignore_mean=True),
                       RV_from_params(time, params, ignore_mean=False) - params["mean_val"])


@given(st.floats(allow_nan=False, allow_infinity=False),
       st.floats(min_value=0.01, allow_nan=False, allow_infinity=False),
       st.floats(min_value=0.001, allow_nan=False, allow_infinity=False))
def test_companion_amplitude_function(k1, m1, m2):
    m_ratio = (M_sun / M_jup).value
    assert np.allclose(companion_amplitude(k1, m1, m2), (-k1 * m1 * m_ratio / m2))


@pytest.mark.xfail(strict=True)
def test_mid_min_max():
    """Test error bars get added correctly.

     For parsing parameters with error bars."""
    param_1 = [10, -5, 2.5]
    assert np.allclose(rv.min_mid_max(param_1), [5, 10, 12.5])
    param_2 = [5.1, 2]
    assert np.allclose(rv.min_mid_max(param_2), [3.1, 5.1, 7.1])
    param_3 = None
    assert rv.min_mid_max(param_3) == [None, None, None]
    param_3 = 1
    assert rv.min_mid_max(param_3) == 1  # No change
    param_4 = [1]
    assert rv.min_mid_max(param_4) == [1]  # No change

    with pytest.raises(ValueError):
        rv.min_mid_max([1, 2, 3])
    with pytest.raises(ValueError):
        rv.min_mid_max([1, -2])
    with pytest.raises(ValueError):
        rv.min_mid_max([1, -2, -3])