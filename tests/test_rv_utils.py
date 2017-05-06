import ephem
import pytest
import datetime
import numpy as np
from hypothesis import strategies as st
from hypothesis import given, example, assume

from utils.rv_utils import mean_anomaly, true_anomaly
from utils.rv_utils import jd2datetime, datetime2jd


# issue with limits  0-pi only
@given(st.lists(st.floats(min_value=0, max_value=np.pi), min_size=1), st.floats(min_value=0.05, max_value=0.99))
def test_trueanomaly(ma, ecc):

    ma = np.asarray(ma)
    assume(np.all(np.abs(ma) > 0.0001))
    ta = true_anomaly(ma, ecc)
    # cos(E) = (e + cos(ta)) / (1 + e*cos(ta))
    E = np.arccos((ecc + np.cos(ta)) /(1 + ecc * np.cos(ta)))
    # ma = E - e*sin(E)
    print("VALUES", ma, (E - ecc * np.sin(E)))
    assert np.allclose(ma, E - ecc * np.sin(E), rtol=0.05)
    assert len(ta) == len(ma)


# issue with limits 0-pi only
@given(st.floats(min_value=0, max_value=np.pi), st.floats(min_value=0.01, max_value=0.99))
@example(2, 0.5)   # example with an integer
def test_trueanomaly_with_scalar(ma, ecc):
    assume(abs(ma) > 0.001)
    ta = true_anomaly(ma, ecc)
    # Contrast Result to E from ta
    # cos(E) = (e + cos(ta)) / (1 + e*cos(ta))
    E = np.arccos((ecc + np.cos(ta)) / (1 + ecc * np.cos(ta)))
    # ma = E - e*sin(E)
    MA = E - ecc * np.sin(E)
    assert np.allclose(ma, E - ecc * np.sin(E))
    assert len(ta) == 1


@given(st.floats(min_value=0.01, max_value=0.99))
def test_trueanomaly_errors(ecc):

    with pytest.raises(TypeError):
        true_anomaly([], ecc)

    with pytest.raises(ValueError):
        true_anomaly(np.array([]), ecc)


@given(st.lists(st.floats(), min_size=1), st.floats(), st.floats(min_value=0.01))
def test_mean_anomaly(t, t0, p):
    """Mean anomaly is an angle, doen't have a constraint value."""
    t = np.array(t)
    ma = mean_anomaly(t, t0, p)

    assert len(t) == len(ma)
    assert isinstance(t, np.ndarray)


@pytest.mark.xfail
def test_radial_velocity():
    assert False
# def radial_velocity(gamma, k, ta, omega, ecc):
#     # Calculate radial velocity of star
#     return gamma + k * (np.cos(ta + omega) + ecc * np.cos(omega))


@pytest.mark.xfail
def test_rv_curve():
    assert False
# def rv_curve_py(times, gamma, k, omega, ecc, t0, period):
#     ma = mean_anomaly(times, t0, period)
#     ta = true_anomaly(ma, ecc)
#     rv = radial_velocity(gamma, k, ta, omega, ecc)
#     return rv


@pytest.mark.parametrize("param_file", ["test/params.txt"])
def parameter_fixture(param_file):
    """Load Parameter file."""
    return parse_paramfile(param_file)


# def test_RV_from_params(parameter_fixture):
def test_RV_from_params_circular():
    """Maximum RV should be within gamma+k1, gamma+k2 for a circular orbit."""
    params = rv.parse_paramfile("tests/test_params.txt")
    params["eccentricity"] = 0  # Circular orbit
    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)
    # amp = gamma + K * sin()
    rvs = RV_from_params(time, params)
    min_val = params["mean_val"] - params["k1"]
    max_val = params["mean_val"] + params["k1"]
    assert np.all(rvs < max_val)
    assert np.all(rvs > min_val)


def test_RV_from_params():
    """RV should be within theroretical limits."""
    params = rv.parse_paramfile("tests/test_params.txt")
    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)
    rvs = RV_from_params(time, params) - params["mean_val"]  # remove center

    A1 = params["k1"] * (1 + params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    B1 = params["k1"] * (1 - params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))

    max_val = A1
    min_val = -B1

    max_rv = np.max(rvs)
    min_rv = np.min(rvs)
    print("max_val", max_val, "min_val", min_val)
    print("max_rv", max_rv, "min_rv", min_rv)
    assert np.all(rvs < max_val)
    assert np.all(rvs > min_val)
    if params["eccentricity"] == 0:
        assert np.allclose(A1, B1)
    else:
        assert not np.allclose(A1, B1)
    assert np.allclose(params["k1"], 0.5 * (A1 + B1))


def test_RV_ignore_mean():
    """Maximum RV should be within gamma+k1, gamma+k2 for a circular orbit."""
    params = rv.parse_paramfile("tests/test_params.txt")
    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)

    assert np.allclose(RV_from_params(time, params, ignore_mean=True),
                       RV_from_params(time, params, ignore_mean=False) - params["mean_val"])


# @pytest.mark.xfail
def test_from_params_companion():
    """Maximum RV should be within gamma+k1, gamma+k2 for a circular orbit."""
    params = rv.parse_paramfile("tests/test_params.txt")
    time = np.linspace(params["tau"], params["tau"] + params["period"], 1000)
    rvs = RV_from_params(time, params, ignore_mean=True, companion=True)

    A2 = params["k2"] * (1 + params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    B2 = params["k2"] * (1 - params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))

    max_val = -B2
    min_val = A2

    max_rv = np.max(rvs)
    min_rv = np.min(rvs)
    print("max_val", max_val, "min_val", min_val)
    print("max_rv", max_rv, "min_rv", min_rv)
    assert np.all(rvs < max_val)
    assert np.all(rvs > min_val)
    if params["eccentricity"] == 0:
        assert np.allclose(A2, B2)
    else:
        assert not np.allclose(A2, B2)
    assert np.allclose(params["k2"], 0.5 * (A2 + B2))


@pytest.mark.parametrize("jd, expected", [
    (2400000.5, (1858, 11, 17)),
    (2458130.1, (2018, 1, 11, 14, 24, 0))])
def test_jd2datetime(jd, expected):
    assert abs(jd2datetime(jd) - datetime.datetime(*expected)) < datetime.timedelta(seconds=1)


@pytest.mark.parametrize("date, expected", [
    ((2012, 2, 12, 11, 31, 10), 2455969.979977),
    ((1990, 9, 6, 20), 2448141.333333)])
def test_datetime2jd(date, expected):
    d = datetime.datetime(*date)
    assert np.allclose(datetime2jd(d), expected)
    assert np.allclose(datetime2jd(d), ephem.julian_date(d))


@given(st.floats(min_value=2200000, max_value=2600000))
def test_jdconversions(jd):
    """Test jd2datetime and datetime2jd are reversable"""

    assert np.allclose(datetime2jd(jd2datetime(jd)), jd)
    assert np.allclose(ephem.julian_date(jd2datetime(jd)), datetime2jd(jd2datetime(jd)))
