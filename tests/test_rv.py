import pytest
import numpy as np
#from ObservationTools import rv
import rv
from rv import obs_time_jd
from hypothesis import given, example
from hypothesis import strategies as st


def test_parse_obslist():
    test_file = "tests/test_obstimes.txt"
    obs_list = rv.parse_obslist(test_file)
    assert isinstance(obs_list, list)
    print(obs_list)
    assert sorted(obs_list) == sorted(["2012-08-14 12:44:05", "2012-09-24 13:12:10"])


def test_parse_params():
    param_file = "tests/test_params.txt"
    params = rv.parse_paramfile(param_file)

    assert isinstance(params, dict)
    assert "name" in params.keys()
    assert params["name"] == "hd30501"
    assert "period" in params.keys()
    assert isinstance(params["period"], float)


@given(st.lists(st.floats(allow_nan=False), min_size=1))
def test_parse_list_string(in_list):
    """Test string of a list turned into list of the values.

    Tranforms to floats if possible."""
    str_list = str(in_list)
    assert rv.parse_list_string(str_list) == in_list


#@example(['source1', 'source2 et. al. 2017']) # not equal.
#@example(in_list="source1, source2 et. al. 2017")
def test_parse_list_string_with_strings():
    """Test list parse with strings.

    Triggers since they have a comma separator.
    """
    in_list="source1, source2 et. al. 2017"
    str_list = str(in_list)
    assert rv.parse_list_string(str_list) == ["source1", "source2 et. al. 2017"]


def test_mid_min_max():
    """Test error bars get added correctly.

     For parsing paramerers with errorbars."""
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


def test_obs_times_jd():
    obs_times = ["2017-05-01", "2015-01-02", "2016-04-05 12:34:15"]
    jd = obs_time_jd(obs_times=obs_times)
    jd_expected_results = [2457874.5, 2457024.5, 2457484.023784722]
    assert np.allclose(jd, jd_expected_results)
    assert np.all([isinstance(j, float) for j in jd])


def test_obs_list_jd():
        obs_list = "tests/test_obstimes.txt"
        jd = obs_time_jd(obs_list=obs_list)
        jd_expected_results = [2456154.030613426, 2456195.050115741]
        assert np.allclose(jd, jd_expected_results)
        assert np.all([isinstance(j, float) for j in jd])


def test_both_obs_times_jd():
        obs_list = "tests/test_obstimes.txt"
        obs_times = ["2017-05-01", "2015-01-02", "2016-04-05 12:34:15"]
        jd = obs_time_jd(obs_times=obs_times, obs_list=obs_list)
        jd_expected_results = [2457874.5, 2457024.5, 2457484.023784722, 2456154.030613426, 2456195.050115741]
        assert np.allclose(jd, jd_expected_results)
        assert np.all([isinstance(j, float) for j in jd])

def test_neither_obs_times_jd():
        obs_list =  None
        obs_times = None
        jd = obs_time_jd(obs_times=obs_times, obs_list=obs_list)
        jd_expected_results = None
        assert jd is jd_expected_results
