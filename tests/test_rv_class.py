from utils.rv_utils import RV
from utils.parse import parse_paramfile


def test_initalize_RV_class():

    rv = RV()
    assert isinstance(rv, RV)


def test_initalize_RV_from_dict():
    params = {"k1": 1, "period": 2, "tau": 5000, "omega": 1, "eccentricity": 0.5, "mean_val": 5}
    rv = RV.from_dict(params)
    assert rv.semi_amp == params["k1"]
    assert rv.period == params["period"]
    assert rv.ecc == params["eccentricity"]
    assert rv.tau == params["tau"]
    assert rv.gamma == params["mean_val"]
    assert rv.omega == params["omega"]


def test_initalize_RV_from_file():
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

