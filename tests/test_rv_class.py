import numpy as np
import pytest
from hypothesis import given, strategies as st, example, assume

from utils.parse import parse_paramfile
from utils.rv_utils import RV


def test_rv_object_instance_of_rv_class():
    rv = RV()
    assert isinstance(rv, RV)


def test_initalize_rv_class_from_dict():
    params = {"k1": 1, "period": 2, "tau": 5000, "omega": 1, "eccentricity": 0.5, "mean_val": 5}
    rv = RV.from_dict(params)
    assert rv.semi_amp == params["k1"]
    assert rv.period == params["period"]
    assert rv.ecc == params["eccentricity"]
    assert rv.tau == params["tau"]
    assert rv.gamma == params["mean_val"]
    assert rv.omega == params["omega"]


def test_initalize_rv_class_from_file():
    paramfile = "tests/test_params.txt"
    rv = RV.from_file(paramfile)
    params = parse_paramfile(paramfile)

    assert rv.semi_amp == params["k1"]
    assert rv.period == params["period"]
    assert rv.ecc == params["eccentricity"]
    assert rv.tau == params["tau"]
    assert rv.gamma == params["mean_val"]
    assert rv.omega == params["omega"]


@pytest.mark.parametrize("semi_amp, period, tau, gamma, omega", [
    (1, 1, 1000, 10, 0),
    (10, 10, 3800, 100, 100),
])
def test_rv_class_max_amp_on_circle(semi_amp, period, tau, gamma, omega):
    ecc = 0
    rv = RV(semi_amp, period, ecc, tau, gamma, omega)
    assert rv.max_amp() == semi_amp


@pytest.mark.parametrize("semi_amp, period, ecc, tau, gamma, omega, expected_amp", [
    (1, 1, 0.25, 1000, 10, 0, 1.25),
    (10, 10, 0.5, 2800, 100, 300, 12.5),
    (20, 10, 0.75, 5800, 100, 180, 35.0),
])
def test_rv_class_max_amp_on_elipse(semi_amp, period, ecc, tau, gamma, omega, expected_amp):
    rv = RV(semi_amp, period, ecc, tau, gamma, omega)
    assert rv.max_amp() <= abs(semi_amp * (1 + ecc))  # omega = 0, 2pi etc
    assert rv.max_amp() == expected_amp


@given(st.floats(min_value=0, max_value=np.pi), st.floats(min_value=0.001, max_value=0.999))
@example(2, 0.5)  # example with an integer
def test_true_anomaly_with_scalar(ma, ecc):
    assume(abs(ma) > 0.001)
    ta = RV.true_anomaly(ma, ecc)
    E = np.arccos((ecc + np.cos(ta)) / (1 + ecc * np.cos(ta)))
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


@given(st.lists(st.floats(min_value=0, max_value=np.pi), min_size=1), st.floats(min_value=0.001, max_value=0.9999))
def test_true_anomaly(ma, ecc):
    ma = np.asarray(ma)
    assume(np.all(np.abs(ma) > 0.0001))
    ta = RV.true_anomaly(ma, ecc)
    E = np.arccos((ecc + np.cos(ta)) / (1 + ecc * np.cos(ta)))

    assert np.allclose(ma, E - ecc * np.sin(E), rtol=0.05)
    assert len(ta) == len(ma)


def test_create_companion():
    host = RV(semi_amp=1, m1=10, m2=5)
    companion = host.create_companion()

    assert companion.semi_amp == -host.semi_amp * host._params["m1"] / host._params["m2"]
    assert companion._params["k2"] == host._params["k1"]
    assert companion._params["k1"] == companion.semi_amp


@pytest.mark.parametrize("k2", [1, 10, 30])
def test_create_companion_with_k2(k2):
    host = RV(semi_amp=1, m1=10, m2=5, k2=k2)
    companion = host.create_companion()

    assert companion.semi_amp == k2  # does not depend on m1 and m2 if k2 given

    assert companion._params["k2"] == host._params["k1"]
    assert companion._params["k1"] == companion.semi_amp


# (Msun, Mjup)
@pytest.mark.parametrize("semi_amp, mass_ratio", [
    (1, 5.0), (0.8, 10), (10.2, 1024)])
def test_create_companion_with_mass_ratio(semi_amp, mass_ratio):
    host = RV(semi_amp=1, k1=1, m1=10, m2=5)
    print("host params", host._params)

    companion = host.create_companion(mass_ratio=mass_ratio)
    print("companion params", companion._params)
    assert companion.semi_amp == (- host.semi_amp * mass_ratio)
    assert host.period == companion.period
    assert host.gamma == companion.gamma
    assert host.ecc == companion.ecc
    assert host.tau == companion.tau

    assert companion._params["k2"] == host._params["k1"]
    assert companion._params["k1"] == companion.semi_amp


def test_double_create_companion_returns_host():
    host = RV(semi_amp=1, m1=10, m2=5)
    companion = host.create_companion()

    host_2 = companion.create_companion()

    assert host_2.semi_amp == host.semi_amp
    assert host == host_2
    assert host != companion


@pytest.mark.parametrize("mass_ratio", [
    5.0, 10, 1024])
def test_double_create_companion_with_ratio_returns_host(mass_ratio):
    host = RV(semi_amp=1, m1=10, m2=5)
    print("host parsm", host._params)
    companion = host.create_companion(mass_ratio=mass_ratio)
    print("comp params", companion._params)
    host_2 = companion.create_companion(1.0 / mass_ratio)
    print("host_2 params", host_2._params)
    assert host_2.semi_amp == host.semi_amp
    assert host == host_2
    assert host != companion


def test_from_dict_works_properly():
    params = {"k1": 10, "eccentricity": 0.5, "period": 5, "mean_val": 4, "tau": 1, "omega": 1, "m1": 4, "m2": 6}

    rv = RV.from_dict(params)
    assert rv._params.get("other_params") is None
    assert rv._params.get("m1") == 4
    assert rv._params.get("m2") == 6


def test_param_from_dict_and_to_dict_give_the_same_dict():
    params = {"k1": 10, "eccentricity": 0.5, "period": 5, "mean_val": 4, "tau": 1, "omega": 1,
              "m1": 4, "m2": 6, "k2": 100, "name": "test"}

    rv1 = RV.from_dict(params)
    rv1_params = rv1.to_dict()
    assert params == rv1_params
    assert params == rv1._params
    assert rv1_params == rv1._params

    rv2 = RV.from_dict(rv1_params)
    assert rv2 == rv1
