import datetime

import ephem
import numpy as np
import pytest
from astropy.constants import M_jup, M_sun
from hypothesis import given
from hypothesis import strategies as st

from utils.rv_utils import (RV, JulianDate, check_core_parameters, generate_companion_label, prepare_mass_params,
                            strtimes2jd)


@pytest.mark.xfail
def test_radial_velocity():
    assert False


@pytest.mark.xfail
def test_rv_curve():
    assert False


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
    assert np.allclose(jd.jd, ephem.julian_date(d) - 2400000)


@pytest.mark.parametrize("julian_date, expected", [
    (2455969.979977, 55969.979977),
    (2448141.333333, 48141.333333)])
def test_JulianDate_reduce_jd(julian_date, expected):
    jd = JulianDate(julian_date)
    assert jd.jd == julian_date
    assert not jd.reduced
    jd.reduce()
    assert jd.reduced
    assert np.allclose(jd.jd, expected)


@pytest.mark.parametrize("expected, julian_date", [
    ((2012, 2, 12, 11, 31, 10), 2455969.979977),
    ((1990, 9, 6, 20), 2448141.333333)])
def test_JulianDate_reduce_jd_to_datetime(expected, julian_date):
    jd = JulianDate(julian_date)
    jd.reduce()
    assert abs(jd.to_datetime() - datetime.datetime(*expected)) < datetime.timedelta(seconds=1)


@pytest.mark.parametrize("date, expected", [
    ("2012-02-12 11:31:10", 2455969.979977),
    ("1990-09-06 20:00:00", 2448141.333333)])
def test_JulianDate_from_str(date, expected):
    jd = JulianDate.from_str(date)
    assert np.allclose(jd.jd, expected)
    assert np.allclose(jd.jd, ephem.julian_date(date))


@pytest.mark.parametrize("datestr, dtime", [
    ("2012-02-12 11:31:10", (2012, 2, 12, 11, 31, 10)),
    ("1990-09-06 20:00:00", (1990, 9, 6, 20))])
def test_JulianDate_from_str_to_datetime(datestr, dtime):
    jd = JulianDate.from_str(datestr)
    assert abs(jd.to_datetime() - datetime.datetime(*dtime)) < datetime.timedelta(seconds=1)


@pytest.mark.parametrize("obstimes, expected",
                         [(["2012-01-05", "2014-04-08"], [2455931.5, 2456755.5]), (["1999-04-08"], [2451276.5]),
                          ("1999-04-08", 2451276.5)])
def test_strtimes2jd(obstimes, expected):
    assert strtimes2jd(obstimes, format="%Y-%m-%d") == expected


@pytest.mark.parametrize("obstimes, expected",
                         [(["2012-01-05", "2014-04-08"], [55931.5, 56755.5]),
                          (["1999-04-08"], [51276.5]),
                          ("1999-04-08", 51276.5)])
def test_strtimes2jd_with_reduce(obstimes, expected):
    assert strtimes2jd(obstimes, reduced=True, format="%Y-%m-%d") == expected


def test_strtimes2jd_with_None():
    assert strtimes2jd(None) == None


@pytest.mark.parametrize("input, expected", [
    (2456755.5, "2014-04-08"),
    (2455931.5, "2012-01-05"),
    (2451276.5, "1999-04-08")
])
def test_julian_date_to_str(input, expected):
    jd = JulianDate(input)
    format = "%Y-%m-%d"
    assert jd.to_str(format) == expected


@pytest.mark.parametrize("input", [
    2456755.5,
    2455931.5,
    2451276.5])
def test_julian_date_to_and_frm_str_starting_from_jd(input):
    assert JulianDate.from_str(JulianDate(input).to_str()).jd == input


@pytest.mark.parametrize("input, format", [
    ("2014-04-08", "%Y-%m-%d"),
    ("2012-01-05", "%Y-%m-%d"),
    ("1999-04-08 12:10:30", "%Y-%m-%d %H:%M:%S")
])
def test_julian_date_to_and_frm_str_starting_from_str(input, format):
    assert JulianDate.from_str(input, format).to_str(format) == input


@pytest.mark.parametrize("removed", [
    "name", "k1", "eccentricity", "omega", "tau", "period"])
def test_check_error_when_core_param_missing(removed):
    params = {"name": "test", "k1": 100, "eccentricity": 0.5, "omega": 20,
              "mean_val": 5, "tau": 5000, "period": 1}
    params.pop(removed)

    with pytest.raises(ValueError):
        check_core_parameters(params)


def test_missing_mean_val_is_set_zero_in_check_core_param():
    params_1 = {"name": "test", "k1": 100, "eccentricity": 0.5,
                "omega": 20, "tau": 5000, "period": 1}
    params_2 = params_1.copy()
    assert params_1.get("mean_val") is None
    params_1 = check_core_parameters(params_1)
    assert params_1["mean_val"] == 0.0

    params_2.update({"mean_val": ""})
    assert params_2["mean_val"] == ""
    params_2 = check_core_parameters(params_2)
    assert params_2["mean_val"] == 0.0


@pytest.mark.parametrize("msini_flag, k2_flag, ratio_flag, expected", [
    (True, True, True, "Mass ratio Companion"),
    (True, True, False, "Given k2 Companion"),
    (True, False, False, "M2sini Companion"),
    (False, False, False, "M2 Companion"),
    (False, True, False, "Given k2 Companion"),
    (False, True, True, "Mass ratio Companion"),
    (False, False, True, "Mass ratio Companion"),
])
def test_generate_companion_label(msini_flag, k2_flag, ratio_flag, expected):
    orbit = RV(msini_flag=msini_flag, k2_flag=k2_flag, ratio_flag=ratio_flag)
    assert generate_companion_label(orbit) == expected


def test_companion_label_with_no_flags():
    assert generate_companion_label(RV()) == "M2 Companion"


@pytest.mark.parametrize("only_msini", [True, False])
def test_prepare_mass_params_sets_msini_flag(only_msini):
    params = {"m_sini": 20, "m_sun": 0.8, "m_true": 35}
    assert prepare_mass_params(params, only_msini=only_msini)["msini_flag"] == only_msini


def test_prepare_mass_params_with_no_m1():
    params = {"m_sini": 20, "m_sun": 0.8, "m_true": 80}

    params = prepare_mass_params(params)
    assert params["m1"] == params["m_sun"] * M_sun / M_jup


@pytest.mark.parametrize("k2, expected", [
    (None, False),
    (1, True)])
def test_prepare_mass_params_with_k2_param(k2, expected):
    params = {"k2": k2, "m_sini": 20, "m_sun": 0.8, "m_true": 80}
    assert params.get("k2_flag") is None
    prepare_mass_params(params)
    assert params.get("k2_flag") == expected


@pytest.mark.parametrize("m2", [5, 10, 20])
def test_prepare_mass_params_with_m2(m2):
    """M2 does not change if given."""
    params = {"m2": m2, "m_sini": 20, "m_sun": 0.8, "m_true": 80}
    assert params["m2"] == m2

    params = prepare_mass_params(params)
    assert params["m2"] == m2


def test_prepare_mass_params_with_msni_m2():
    params = {"m_sini": 20, "m_sun": 0.8, "m_true": 80}
    params = prepare_mass_params(params, only_msini=True)
    assert params["m2"] == params["m_sini"]
    assert params["msini_flag"] is True


def test_prepare_mass_params_with_m_true_m2():
    params = {"m_sini": 20, "m_sun": 0.8, "m_true": 80}
    params = prepare_mass_params(params, only_msini=False)
    assert params["m2"] == params["m_true"]
    assert params["msini_flag"] is False


def test_prepare_mass_params_scales_m1_to_jup_mass():
    params = {"m1": 1., "m2": 0}
    params = prepare_mass_params(params)
    assert params["m1"] == M_sun / M_jup
