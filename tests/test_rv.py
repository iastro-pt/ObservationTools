import pytest
import numpy as np
#from ObservationTools import rv
import rv
from rv import obs_time_jd


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
