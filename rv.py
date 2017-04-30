"""Radial Velocity calculations.

Goals
-----
To calcualte when the radial velocity is different by a certian value.
Plot radial velocity phase curves. Indicating obtained measurement locations.

"""
import os
import ephem
import argparse
import numpy as np
import logging
from logging import debug
from typing import Dict, List
from utils.rv_utils import RV_from_params
# try:
#     from ajplanet import pl_rv_array
#     use_ajplanet = False
# except:
#     use_ajplanet = False

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
try:
    from utils_debug import pv
except ImportError:
    from utils.utils_debug import pv


def _parser():
    parser = argparse.ArgumentParser(description='Radial velocity plotting')
    parser.add_argument('params', help='RV parameters filename')
    # parser.add_argument('-d', '--date', default='today',
    #                     help='Date in format YYYY-MM-DD. Default is today.')
    # #  prediciting separated lines
    # parser.add_argument('-r', '--rv_diff', help='RV difference threshold to find')
    parser.add_argument('-o', '--obs_times', help='Times of previous observations YYYY-MM-DD format',
                        nargs='+', default=None)
    parser.add_argument('-l', '--obs_list', help='File with list of obs times.', type=str, default=None)
    # parser.add_argument('-f', '--files', help='Params and obs-times are file'
    #                    ' names to open', action='store_true')
    parser.add_argument('-m', '--mode', help='Display mode '
                        ' e.g. phase or time plot. Default="phase"',
                        choices=['time', 'phase'], default='phase')
    parser.add_argument("--debug", help="Turning on debug output", action='store_true', default=False)
    return parser.parse_args()


def companion_amplitude(k_host: float, m_host: float, m_companion: float) -> float:
    """Calcualte the companion RV maximum amplitude.

    Parameters
    ----------
    k_host: float
        Amplitude of radial velocity variation of host.
    m_host: float
        Mass of host
    m_companion: float
        Mass of companion in consistent units.

    Returns
    -------
    k_companion: float
        RV amplitude of companion.

    """
    sun_jupiter_mass = 1047.56  # Solar mass in jupiter masses
    m_host *= sun_jupiter_mass  # Convert to jupiter mass
    return -k_host * m_host / m_companion


def parse_obslist(fname: str, path: str=None) -> List[str]:
    # Parse Obslist file containing list of dates/times
    """Parse Obslist file containing list of dates/times.

    Parameters
    ----------
    fname: str
        Filename of obs_list file.
    path: str [optional]
        Path to directory of filename.

    Returns
    --------
    times: list of strings
        Observation times in a list.
    """
    if path is not None:
        fname = os.path.join(path, fname)
    if not os.path.exists(fname):
        logging.warning("Obs_list file given does not exist. {}".format(fname))

    obstimes = list()
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                if "#" in line:   # Remove comment from end of line
                    line = line.split("#")[0]
                if "." in line:
                    line = line.split(".")[0]   # remove fractions of seconds.
                obstimes.append(line.strip())
        debug(pv("obstimes"))
    return obstimes


def parse_paramfile(param_file: str, path: str=None) -> Dict:
    """Extract orbit and stellar parameters from parameter file.

    Parameters
    ----------
    param_file: str
        Filename of parameter file.
    path: str [optional]
        Path to directory of filename.

    Returns
    --------
    parameters: dict
        Paramemters as a {param: value} dictionary.
    """
    if path is not None:
        param_file = os.path.join(path, param_file)
    parameters = dict()
    if not os.path.exists(param_file):
        logging.warning("Parameter file given does not exist. {}".format(param_file))

    with open(param_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                if '#' in line:   # Remove comment from end of line
                    line = line.split("#")[0]
                if line.endswith("="):
                    logging.warning(("Parameter missing value in {}.\nLine = {line}."
                                    " Value set to None.").format(param_file, line))
                    line = line + " None"   # Add None value when parameter is missing
                par, val = line.lower().split('=')
                par, val = par.strip(), val.strip()
                if (val.startswith("[") and val.endswith("]")) or ("," in val):  # Val is a list
                    parameters[par] = parse_list_string(val)
                else:
                    try:
                        parameters[par] = float(val)  # Turn parameters to floats if possible.
                    except ValueError:
                        parameters[par] = val

    return parameters


def parse_list_string(string):
    """Parse list of floats out of a string."""
    string = string.replace("[", "").replace("]", "").strip()
    list_str = string.split(",")
    list_str = [float(val) for val in list_str]
    return list_str


def obs_time_jd(obs_times=None, obs_list=None):
    """Combine observation dates and turn to jd.

    Parameters
    ----------
    obs_times: list of str
        List of dates entered at command line.
    obs_list: str or None
        Filename to observation list.

    Returns
    -------
    Dates: list of floats
        Combined dates converted to julian dates.
    """
    if obs_list is not None:
        obs_list_vals = parse_obslist(obs_list)
        debug(pv("obs_list_vals"))
        if obs_times is None:
            obs_times = obs_list_vals
        else:
            obs_times = obs_times + obs_list_vals

    debug(pv("obs_times"))
    jds = [ephem.julian_date(t) for t in obs_times]
    debug(pv("jds"))
    return jds


def main(params, mode="phase", obs_times=None, obs_list=None):  # obs_times=None, mode='phase', rv_diff=None
    r"""Radial velocity displays.

    Parameters
    ----------
    params: str
        Filename for text file containing the rv parameters. Format of 'param = value\n'.
    mode: str
        Mode for script to use. Phase, time, future.
    obs_times: list of str
        Dates of observations added manually at comand line of format YYYY-MM-DD.
    obs_list: str
        Filename for list which contains obs_times (YYY-MM-DD HH:MM:SS).
    """
    only_msini = True   # Use only the msini values not m_true.
    # Load in params and store as a dictionary
    parameters = parse_paramfile(params)

    # Test of core parameters
    for key in ["name", "k1", "eccentricity", "omega", "tau", "period"]:
        if key not in parameters.keys():
            raise ValueError("Core RV parameter not provided in param file, '{}'".format(key))

    if "mean_val" not in parameters.keys():
        logging.info("mean_val parameter was not provided so set to 0 km/s")
        parameters["mean_val"] = 0.0

    # combine obs_times and obs_list and turn into jd.
    obs_jd = obs_time_jd(obs_times, obs_list)

    # Calculate companion semi-major axis
    if mode in ("error", "indiv"):
        pass
    else:
        if "k2" in parameters.keys():
            pass
        else:
            if ('m_true' in parameters.keys()) and not only_msini:
                # Use true mass if given
                parameters['k2'] = companion_amplitude(parameters['k1'],
                                                       parameters['m_star'],
                                                       parameters['m_true'])
            else:
                parameters['k2'] = companion_amplitude(parameters['k1'],
                                                       parameters['m_star'],
                                                       parameters['msini'])

    if mode == "phase":
        RV_phase_curve(parameters)
    elif mode == "specdiff":
        differential_analysis(parameters, obs_jd)
    elif mode == "error":
        error_analysis(parameters)
    elif mode == "indiv":
        individual_errorbar(parameters)
    elif mode == "diff":
        # simple diff
        rv_diff(parameters, obs_jd)
    else:
        raise NotImplementedError("Other modes are not Implemented yet.")


def individual_errorbar(params):
    """Make plots just varying one parameter at a time."""
    for var_par in ["k1", "m_star", "period", "msini"]:
        p = [params["period"][0] if isinstance(params["period"], list) else params["period"]]
        today = ephem.julian_date(ephem.now())
        t_space = np.linspace(today, today + p[0], 100)
        plt.figure()
        rv_params = {}
        for key in params:
            if key != var_par:
                if isinstance(params[key], list):
                    rv_params[key] = params[key][0]
                else:
                    rv_params[key] = params[key]
        for value in min_mid_max(params[var_par]):
            rv_params[var_par] = value
            rv_params["k2"] = companion_amplitude(rv_params['k1'],
                                                  rv_params['m_star'],
                                                  rv_params['msini'])

            plt.plot(t_space, RV_from_params(t_space, rv_params, ignore_mean=True, companion=True))
            # plt.plot(t_space, RV_from_params(t_space, rv_params, ignore_mean=True, companion=False))
        plt.xlim()
        plt.ylabel("Companion RV")
        plt.xlabel("JD")
        plt.title("Bounds analysis for {0!s}".format(var_par))
        plt.legend()
        plt.savefig("Test_errors_{}.png".format(var_par))
    plt.show()


def error_analysis(params=None, obs_times=None):
    """Each parameter has error bars (or not).

    [val, lower, upper].
    """
    # simple test of altering one paramerter
    today = ephem.julian_date(ephem.now())
    t_p = params["period"]
    if isinstance(t_p, list):
        t_p = t_p[0]

    t_space = np.linspace(today, today + t_p, 100)
    if "k2" not in params.keys():
        params["k2"] = None
    # need to reduce number of parameters 4 at a time?.

    tau = params["tau"]
    tau = (tau[0] if isinstance(tau, list) else tau)
    mean_val = params["mean_val"]
    mean_val = (mean_val[0] if isinstance(mean_val, list) else mean_val)
    period = params["period"]
    period = (period[0] if isinstance(period, list) else period)

    p_list = []
    for key in ["k1", "eccentricity", "omega", "msini", "m_true", "m_star"]:
        p_list.append(min_mid_max(params[key]))

    error_iter = itertools.product(*p_list)
    for i, (kk, ee, ww, mm, mt, ms) in enumerate(error_iter):

        # Calculate k2
        kk2 = companion_amplitude(kk, ms, mm)

        these_params = {"msini": mm, "m_true": mt, "tau": tau, "period": period,
                        "k2": kk2, "k1": kk, "mean_val": mean_val, "omega": ww,
                        "eccentricity": ee}
        debug(pv("(kk, ee, ww, mm, mt, ms, kk2)"))
        comp_rvs = RV_from_params(t_space, these_params, ignore_mean=True, companion=True)

        plt.plot(t_space, comp_rvs, alpha=0.5, label="{}".format(i))
        debug(pv("i"))
    plt.legend()
    plt.title("test of many lines")
    plt.savefig("Test_errors_many_params.png")
    plt.show()

    return 0


def min_mid_max(param: List[float]) -> List[float]:
    """Param with error bars.

    Asolute values taken incease of error in parameter 2.
    """
    if param is None:
        return [None, None, None]

    if isinstance(param, (int, float, np.float32, np.float64, np.int32, np.int64)):
        # No errorbars
        return param
    elif len(param) == 3:
        # Unequal errorbars
        if (param[1] > 0):
            raise(ValueError("Lower error bar should be negative.  x - dx"))
        if param[2] < 0:
            raise(ValueError("Lower error bar value should be negative.  x +/- dx"))
        return [param[0] + param[1], param[0], param[0] + param[2]]
    elif len(param) == 2:
        # Equal errorbars
        if param[1] < 0:
            raise(ValueError("Equal error bar should be given positive.  x +/- dx"))
        return [param[0] - param[1], param[0], param[0] + param[1]]
    else:
        return param


def RV_phase_curve(params: Dict, cycle_fraction: float=1, ignore_mean: bool=False, t_past=False, t_future=False) -> int:
    """Plot RV phase curve centered on zero.

    Parameters
    ----------
    params: dict
        Parameters of system.
    cycle_fraction: float
        Fraction of phase space to plot. Default=1
    ignore_mean: bool
        Remove the contribution from the systems mean velocity.
    #t_vals: float, array-like
        Times of past observations.
    t_future: float, array-like
        Times of future observations.

    Returns
    -------
    exit_status: bool
        Returns 0 if successful.

        Displays matplotlib figure.
    """
    phase = np.linspace(-0.5, 0.5, 100) * cycle_fraction
    t = params["tau"] + phase * params["period"]
    # t = params[4] + phase * params[5]

    host_rvs = RV_from_params(t, params, ignore_mean=ignore_mean)
    companion_rvs = RV_from_params(t, params, ignore_mean=ignore_mean, companion=True)

    fig = plt.figure(figsize=(10, 7))
    fig.subplots_adjust()
    ax1 = host_subplot(111)

    ax1.plot(phase, host_rvs, label="Host", lw=2, color="k")
    ax1.set_xlabel("Orbital Phase")
    ax1.set_ylabel("Host RV (km/s)")

    ax2 = ax1.twinx()
    ax2.plot(phase, companion_rvs, '--', label="Companion", lw=2)
    ax2.set_ylabel("Companion RV (km/s)")

    if 'name' in params.keys():
        plt.title("RV Phase Curve for {}".format(params['name'].upper()))
    else:
        plt.title("RV Phase Curve")
    # if t_vals:
    #    for t_num, t_val in enumerate(t_vals):
    #         phi = ((t_val - params[4])/params[5] - 0.5) % 1 + 0.5
    #         rv_star = RV_from_params(t_val, params, ignore_mean=False)
    #         rv_planet = RV_from_params(t_val, params, ignore_mean=False, companion=True)
    #         ax1.plot(phi, rv_star, ".", markersize=12, markeredgewidth=3)
    #         ax2.plot(phi, rv_planet, ".", markersize=12, markeredgewidth=3)

    # if t_future:
    #     for t_num, t_val in enumerate(t_future):
    #         phi = ((t_val - params[4])/params[5] - 0.5) % 1 + 0.5
    #         rv_star = RV_from_params(t_val, params, ignore_mean=False)
    #         rv_planet = RV_from_params(t_val, params, ignore_mean=False, companion=True)
    #         ax1.plot(phi, rv_star, "+", markersize=12, markeredgewidth=3)
    #         ax2.plot(phi, rv_planet, "+", markersize=12, markeredgewidth=3)

    plt.legend(loc=0)
    plt.show()
    return 0


if __name__ == '__main__':
    args = vars(_parser())
    # star_name = args.pop('star_name')
    debug_on = args.pop('debug')
    opts = {k: args[k] for k in args}

    # debug = args.pop('debug')
    if debug_on:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(message)s')
    else:
        logging.basicConfig(level=logging.WARNING,
                            format='%(asctime)s %(levelname)s %(message)s')

    main(**opts)
