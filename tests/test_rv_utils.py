import ephem
import pytest
import datetime
import numpy as np
from hypothesis import strategies as st
from hypothesis import given, example, assume, settings

import rv
from utils.rv_utils import RV, JulianDate
from utils.rv_utils import RV_from_params


# issue with limits  0-pi only
@given(st.lists(st.floats(min_value=0, max_value=np.pi), min_size=1), st.floats(min_value=0.05, max_value=0.99))
def test_true_anomaly(ma, ecc):

    ma = np.asarray(ma)
    assume(np.all(np.abs(ma) > 0.0001))
    ta = RV.true_anomaly(ma, ecc)
    # cos(E) = (e + cos(ta)) / (1 + e*cos(ta))
    E = np.arccos((ecc + np.cos(ta)) /(1 + ecc * np.cos(ta)))
    # ma = E - e*sin(E)
    print("VALUES", ma, (E - ecc * np.sin(E)))
    assert np.allclose(ma, E - ecc * np.sin(E), rtol=0.05)
    assert len(ta) == len(ma)


# issue with limits 0-pi only
@settings(max_examples=50)
@given(st.floats(min_value=0, max_value=np.pi), st.floats(min_value=0.01, max_value=0.99))
@example(2, 0.5)   # example with an integer
def test_true_anomaly_with_scalar(ma, ecc):
    assume(abs(ma) > 0.001)
    ta = RV.true_anomaly(ma, ecc)
    # Contrast Result to E from ta
    # cos(E) = (e + cos(ta)) / (1 + e*cos(ta))
    E = np.arccos((ecc + np.cos(ta)) / (1 + ecc * np.cos(ta)))
    # ma = E - e*sin(E)
    MA = E - ecc * np.sin(E)
    assert np.allclose(ma, E - ecc * np.sin(E))
    assert len(ta) == 1


@given(st.floats(min_value=0.01, max_value=0.99))
def test_true_anomaly_errors(ecc):

    with pytest.raises(TypeError):
        RV.true_anomaly([], ecc)

    with pytest.raises(ValueError):
        RV.true_anomaly(np.array([]), ecc)


@given(st.lists(st.floats(), min_size=1), st.floats(), st.floats(min_value=0.01))
def test_mean_anomaly(t, t0, p):
    """Mean anomaly is an angle, doesn't have a constraint value."""
    t = np.array(t)
    ma = RV.mean_anomaly(t, t0, p)

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
#     ma = RV.mean_anomaly(times, t0, period)
#     ta = RV.true_anomaly(ma, ecc)
#     rv = radial_velocity(gamma, k, ta, omega, ecc)
#     return rv


@pytest.fixture(params=["tests/test_params.txt"])
def params(request):
    """Load Parameter file."""
    return rv.parse_paramfile(request.param)


def test_RV_from_params_circular(params):
    """Maximum RV should be within gamma + k1, gamma + k2 for a circular orbit."""
    params["eccentricity"] = 0  # Circular orbit
    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)
    rvs = RV_from_params(time, params)
    min_val = params["mean_val"] - params["k1"]
    max_val = params["mean_val"] + params["k1"]
    assert np.all(rvs <= max_val)
    assert np.all(rvs >= min_val)


# mean_val k1 period tau omega eccentricity
@settings(max_examples=100)
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


def test_RV_ignore_mean(params):
    """Maximum RV should be within gamma+k1, gamma+k2 for a circular orbit."""
    time = np.linspace(params["tau"], params["tau"] + params["period"], 200)

    assert np.allclose(RV_from_params(time, params, ignore_mean=True),
                       RV_from_params(time, params, ignore_mean=False) - params["mean_val"])


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


@pytest.mark.parametrize("jd, expected", [
    (2400000.5, (1858, 11, 17)),
    (2458130.1, (2018, 1, 11, 14, 24, 0))])
def test_JulianDate_to_datetime(jd, expected):
    jd = JulianDate(jd)
    assert abs(jd.to_datetime() - datetime.datetime(*expected)) < datetime.timedelta(seconds=1)
    assert isinstance(jd, JulianDate)
    assert isinstance(jd.to_datetime(), datetime.datetime)


@pytest.mark.parametrize("jd, expected", [
    (2400000.5, (1858, 11, 17)),
    (2458130.1, (2018, 1, 11, 14, 24, 0))])
def test_reduced_JulianDate_to_datetime(jd, expected):
    jd = JulianDate(jd)
    jd.reduce()
    assert abs(jd.to_datetime() - datetime.datetime(*expected)) < datetime.timedelta(seconds=1)
    assert isinstance(jd, JulianDate)
    assert isinstance(jd.to_datetime(), datetime.datetime)


@pytest.mark.parametrize("jd, expected", [
    (00000.5, (1858, 11, 17)),
    (58130.1, (2018, 1, 11, 14, 24, 0))])
def test_prereduced_JulianDate_to_datetime(jd, expected):
    jd = JulianDate(jd, reduced=True)
    jd.reduce()
    assert abs(jd.to_datetime() - datetime.datetime(*expected)) < datetime.timedelta(seconds=1)
    assert isinstance(jd, JulianDate)
    assert isinstance(jd.to_datetime(), datetime.datetime)


@pytest.mark.parametrize("date, expected", [
    ((2012, 2, 12, 11, 31, 10), 2455969.979977),
    ((1990, 9, 6, 20), 2448141.333333)])
def test_JulianDate_from_datetime(date, expected):
    d = datetime.datetime(*date)
    jd = JulianDate.from_datetime(d)
    assert np.allclose(jd.jd, expected)
    assert np.allclose(jd.jd, ephem.julian_date(d))


@given(st.floats(min_value=2200000, max_value=2600000))
def test_jd_conversions(jd):
    """Test JulianDate.to_datetime and JulianDate.from_datetime are reversible."""
    assert np.allclose(jd,
                       JulianDate(jd).jd)
    assert np.allclose(JulianDate(jd).jd,
                       JulianDate.from_datetime(JulianDate(jd).to_datetime()).jd)
    assert np.allclose(ephem.julian_date(JulianDate(jd).to_datetime()),
                       JulianDate.from_datetime(JulianDate(jd).to_datetime()).jd)


@pytest.mark.parametrize("date, expected", [
    ((2012, 2, 12, 11, 31, 10), 2455969.979977),
    ((1990, 9, 6, 20), 2448141.333333)])
def test_JulianDate_from_datetime(date, expected):
    d = datetime.datetime(*date)
    jd = JulianDate.from_datetime(d)
    assert np.allclose(jd.jd, expected)
    assert np.allclose(jd.jd, ephem.julian_date(d))


@pytest.mark.parametrize("date, expected", [
    ((2012, 2, 12, 11, 31, 10), 55969.979977),
    ((1990, 9, 6, 20), 48141.333333)])
def test_JulianDate_reduce_datetime(date, expected):
    d = datetime.datetime(*date)
    jd = JulianDate.from_datetime(d, reduced=True)

    assert jd.reduced
    assert np.allclose(jd.jd, expected)
    assert np.allclose(jd.jd, ephem.julian_date(d)-2400000)


@pytest.mark.parametrize("juliandate, expected", [
    (2455969.979977, 55969.979977),
    (2448141.333333, 48141.333333)])
def test_JulianDate_reduce_jd(juliandate, expected):
    jd = JulianDate(juliandate)
    assert jd.jd == juliandate
    assert not jd.reduced
    jd.reduce()
    assert jd.reduced
    assert np.allclose(jd.jd, expected)


@pytest.mark.parametrize("expected, juliandate", [
    ((2012, 2, 12, 11, 31, 10), 2455969.979977),
    ((1990, 9, 6, 20), 2448141.333333)])
def test_JulianDate_reduce_jd_to_datetime(expected, juliandate):
    jd = JulianDate(juliandate)
    jd.reduce()
    assert abs(jd.to_datetime() - datetime.datetime(*expected)) < datetime.timedelta(seconds=1)


@pytest.mark.parametrize("date, expected", [
    ("2012-2-12 11:31:10", 2455969.979977),
    ("1990-9-6 20:00:00", 2448141.333333)])
def test_JulianDate_from_str(date, expected):
    jd = JulianDate.from_str(date)
    assert np.allclose(jd.jd, expected)
    assert np.allclose(jd.jd, ephem.julian_date(date))


@pytest.mark.parametrize("datestr, dtime", [
    ("2012-2-12 11:31:10", (2012, 2, 12, 11, 31, 10)),
    ("1990-9-6 20:00:00", (1990, 9, 6, 20))])
def test_JulianDate_from_str_to_datetime(datestr, dtime):
    jd = JulianDate.from_str(datestr)
    assert abs(jd.to_datetime() - datetime.datetime(*dtime)) < datetime.timedelta(seconds=1)


from utils.rv_utils import strtimes2jd
@pytest.mark.parametrize("obstimes, expected", [(["2012-1-5", "2014-4-8"], [2455931.5, 2456755.5]), (["1999-4-8"], [2451276.5]), ("1999-4-8", 2451276.5)])
def test_str2jd_change(obstimes, expected):
    assert strtimes2jd(obstimes) == expected
