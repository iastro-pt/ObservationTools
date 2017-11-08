import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

import rv
import utils.rv_utils
from utils import parse


def test_parse_obslist():
    test_file = "tests/test_obstimes.txt"
    obs_list = parse.parse_obslist(test_file)
    assert isinstance(obs_list, list)
    print(obs_list)
    assert sorted(obs_list) == sorted(["2012-08-14 12:44:05", "2012-09-24 13:12:10"])


def test_parse_params():
    param_file = "tests/test_params.txt"
    params = parse.parse_paramfile(param_file)

    assert isinstance(params, dict)
    assert "name" in params.keys()
    assert params["name"] == "test"
    assert "period" in params.keys()
    assert isinstance(params["period"], float)


@given(st.lists(st.floats(allow_nan=False, allow_infinity=False), min_size=1))
def test_parse_list_string(in_list):
    """Test string of a list turned into list of the values.

    Transforms to floats if possible."""
    str_list = str(in_list)
    assert parse.parse_list_string(str_list) == in_list


# @example(['source1', 'source2 et. al. 2017']) # not equal.
# @example(in_list="source1, source2 et. al. 2017")
def test_parse_list_string_with_strings():
    """Test list parse with strings.

    Triggers since they have a comma separator.
    """
    in_list = "source1, source2 et. al. 2017"
    str_list = str(in_list)
    assert parse.parse_list_string(str_list) == ["source1", "source2 et. al. 2017"]


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


@pytest.mark.parametrize("times,obs_list,expected", [
    (None, None, None),
    (None, "tests/test_obstimes.txt", ["2012-08-14 12:44:05", "2012-09-24 13:12:10"]),
    (["2017-05-01", "2015-01-02", "2016-04-05 12:34:15"], None, ["2017-05-01", "2015-01-02", "2016-04-05 12:34:15"]),
    (["2017-05-01", "2015-01-02", "2016-04-05 12:34:15"], "tests/test_obstimes.txt",
     ["2017-05-01", "2015-01-02", "2016-04-05 12:34:15", "2012-08-14 12:44:05", "2012-09-24 13:12:10"]),
])
def test_join_times(times, obs_list, expected):
    assert utils.rv_utils.join_times(times, obs_list) == expected


@pytest.mark.parametrize("times, format, expected_jd", [
    ([], None, []),
    (["2017-05-01 00:00:00", "2015-01-02 00:00:00", "2016-04-05 12:34:15"], None,
     [2457874.5, 2457024.5, 2457484.023784722]),
    (["2017-05-01", "2015-01-02", "2016-04-05"], "%Y-%m-%d", [2457874.5, 2457024.5, 2457483.5]),
    (["2012-08-14 12:44:05", "2012-09-24 13:12:10"], "%Y-%m-%d %H:%M:%S", [2456154.030613426, 2456195.050115741]),
])
def test_strtimes2jd_with_formats(times, format, expected_jd):
    assert np.allclose(utils.rv_utils.strtimes2jd(times, reduced=False, format=format), expected_jd)


@pytest.mark.parametrize("times,expected_jd", [
    ([], []),
    (["2017-05-01 00:00:00", "2015-01-02 00:00:00", "2016-04-05 12:34:15"], [57874.5, 57024.5, 57484.023784722]),
    (["2012-08-14 12:44:05", "2012-09-24 13:12:10"], [56154.030613426, 56195.050115741]),
])
def test_reduced_strtimes2jd(times, expected_jd):
    jd = utils.rv_utils.strtimes2jd(times, reduced=True)
    assert np.allclose(jd, expected_jd)

#
# @pytest.mark.matlotlib_image_compare
# def test_default_rv_phase_plot():
#     return rv.main("data/HD30501_params.txt")
#
#
# @pytest.mark.matlotlib_image_compare
# def test_rv_phase_plot_with_obs():
#     return rv.main("data/HD30501_params.txt")
#
#
# @pytest.mark.matlotlib_image_compare
# def test_default_rv_time_plot():
#     return rv.main("data/HD30501_params.txt", mode="time")
#
#
# @pytest.mark.matlotlib_image_compare
# def test_rv_time_plot_with_obs():
#     return rv.main("data/HD30501_params.txt", mode="time")
