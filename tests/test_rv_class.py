import pytest

from utils.rv_utils import RV
from utils.parse import parse_paramfile


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
    # params = {"k1": 1, "period": 2, "tau":5000, "omega": 1, "eccentricity": 0.5, "mean_val": 5}
    rv = RV.from_file(paramfile)
    params = parse_paramfile(paramfile)
    print(rv)
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
    assert rv.max_amp() <= abs(semi_amp * (1 + ecc))   # omega = 0, 2pi etc
    assert rv.max_amp() == expected_amp
