"""Radial Velocity calculations.

Goals
-----
To calculate when the radial velocity is different by a certain value.
Plot radial velocity phase curves. Indicating obtained measurement locations.

"""
import argparse
import logging
from datetime import datetime
from typing import Dict, List, Any, Union

import astropy.units as u
import ephem
import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import c
from mpl_toolkits.axes_grid1 import host_subplot

from utils.parse import parse_paramfile
from utils.rv_utils import companion_amplitude, RV_from_params, strtimes2jd, join_times
from utils.rv_utils import RV, JulianDate

# try:
#     from ajplanet import pl_rv_array
#     use_ajplanet = False
# except:
#     use_ajplanet = False

c_km_s = c.to(u.kilometer / u.second)  # Speed of light in km/s


def _parser():
    # type: () -> argparse.Namespace
    """RV Argparse parser."""
    parser = argparse.ArgumentParser(description='Radial velocity plotting')
    parser.add_argument('params', help='RV parameters filename', type=str)
    parser.add_argument('-d', '--date', default=None,
                        help='Reference date in format YYYY-MM-DD [HH:MM:SS]. Default=None uses time of now.')
    # pre-predicting parated lines
    # parser.add_argument('-r', '--rv_diff', help='RV difference threshold to find')
    parser.add_argument('-o', '--obs_times', help='Times of previous observations YYYY-MM-DD format',
                        nargs='+', default=None)
    parser.add_argument('-l', '--obs_list', help='File with list of obs times.', type=str, default=None)
    # parser.add_argument('-f', '--files', help='Params and obs-times are file'
    #                    ' names to open', action='store_true')
    parser.add_argument('-m', '--mode', help='Display mode '
                        ' e.g. phase or time plot. Default="phase"',
                        choices=['phase', 'time'], default='phase', type=str)
    parser.add_argument("--debug", help="Turning on debug output", action='store_true', default=False)
    return parser.parse_args()



def main(params, mode="phase", obs_times=None, obs_list=None, date=None):  # obs_times=None, mode='phase', rv_diff=None
    # type: (Dict[str, Union[str, float]], str, List[str], str, str) -> None
    r"""Radial velocity displays.

    Parameters
    ----------
    params: str
        Filename for text file containing the rv parameters. Format of 'param = value\n'.
    mode: str
        Mode for script to use. Phase, time, future.
    obs_times: list of str
        Dates of observations added manually at command line of format YYYY-MM-DD.
    obs_list: str
        Filename for list which contains obs_times (YYY-MM-DD HH:MM:SS).
    date: str
        Reference date for some modes. Defaults=None)
    """
    # only_msini = True   # Use only the msini values not m_true.
    only_msini = False
    # Load in params and store as a dictionary
    parameters = parse_paramfile(params)

    # Test of core parameters
    for key in ["name", "k1", "eccentricity", "omega", "tau", "period"]:
        if key not in parameters.keys():
            raise ValueError("A core parameter was not provided in param file, '{}'".format(key))

    if "mean_val" not in parameters.keys():
        logging.info("mean_val parameter was not provided so set to 0 km/s")
        parameters["mean_val"] = 0.0
    elif parameters["mean_val"] == "":
        logging.info("mean_val parameter was blank so set to 0 km/s")
        parameters["mean_val"] = 0.0

    # combine obs_times and obs_list and turn into jd.
    if obs_times:
        if (".txt" in obs_times) or (".dat" in obs_times):
            raise ValueError("Filename given instead of dates for obs_times.")

    obs_times = join_times(obs_times, obs_list)
    obs_jd = strtimes2jd(obs_times) # , reduced=True


    # Calculate companion semi-major axis
    if mode in ("error", "indiv"):
        pass
    else:
        if "k2" in parameters.keys():
            pass
        else:
            if ('m_true' in parameters.keys()) and not only_msini:
                # Use true mass if given
                if not parameters["m_true"] == "":
                    mass_used = parameters['m_true']
                    true_mass_flag = True   # parameter to indicate if the true mass was used or not
                else:
                    mass_used = parameters["msini"]
                    true_mass_flag = False
            else:
                mass_used = parameters["msini"]
                true_mass_flag = False

            parameters['k2'] = companion_amplitude(parameters['k1'],
                                                   parameters['m_star'],
                                                   mass_used)
            parameters["true_mass_flag"] = true_mass_flag   # True if true mass used

    if mode == "phase":
            RV_phase_curve(parameters, t_past=obs_jd)
    elif mode == "time":
        if date is not None:
            date = ephem.julian_date(date)
        RV_time_curve(parameters, t_past=obs_jd, start_day=date)

    else:
        raise NotImplementedError("Other modes are not Implemented yet.")


def RV_phase_curve(params, cycle_fraction=1, ignore_mean=False, t_past=False, t_future=False):
    # type: (Dict[str, Union[str, float]], float, bool, Any, Any) -> int
    """Plot RV phase curve centered on zero.

    Parameters
    ----------
    params: dict
        Parameters of system.
    cycle_fraction: float
        Fraction of phase space to plot. (Default=1)
    ignore_mean: bool
        Remove the contribution from the systems mean velocity. (Default=False)
    t_past: float, array-like
        Times of past observations.
    t_future: float, array-like
        Times of future observations.

    Returns
    -------
    exit_status: bool
        Returns 0 if successful.

        Displays matplotlib figure.
    """

    host = RV.from_dict(params)
    companion = RV.from_dict(params)
    companion.semi_amp = params["k2"]
    phase = np.linspace(-0.5, 0.5, 100) * cycle_fraction
    t = params["tau"] + phase * params["period"]

    host_rvs = RV_from_params(t, params, ignore_mean=ignore_mean)
    companion_rvs = RV_from_params(t, params, ignore_mean=ignore_mean, companion=True)

    fig = plt.figure(figsize=(10, 7))
    fig.subplots_adjust()
    ax1 = host_subplot(111)

    ax1.plot(phase, host_rvs, label="Host", lw=2, color="k")
    ax1.set_xlabel("Orbital Phase")
    ax1.set_ylabel("Host RV (km/s)")

    ax2 = ax1.twinx()
    if params["true_mass_flag"]:
        ax2.plot(phase, companion_rvs, '--', label="M_2 Companion", lw=2)
    else:
        ax2.plot(phase, companion_rvs, '--', label="M_2sini Companion", lw=2)
    ax2.set_ylabel("Companion RV (km/s)")

    if 'name' in params.keys():
        if "companion" in params.keys():
            plt.title("RV Phase Curve for {} {}".format(params['name'].upper(), params['companion']))
        else:
            plt.title("RV Phase Curve for {}".format(params['name'].upper()))
    else:
        plt.title("RV Phase Curve")
    if t_past:
        for t_num, t_val in enumerate(t_past):
            phi = ((t_val - params["tau"]) / params["period"] + 0.5) % 1 - 0.5
            rv_star = RV_from_params(t_val, params, ignore_mean=False)
            rv_planet = RV_from_params(t_val, params, ignore_mean=False, companion=True)
            ax1.plot(phi, rv_star, ".", markersize=10, markeredgewidth=2)
            ax2.plot(phi, rv_planet, "+", markersize=10, markeredgewidth=2)

    if t_future:
        raise NotImplementedError("Adding future observations times is not implemented yet.")
    # if t_future:
    #     for t_num, t_val in enumerate(t_future):
    #         phi = ((t_val - params[4])/params[5] - 0.5) % 1 + 0.5
    #         rv_star = RV_from_params(t_val, params, ignore_mean=False)
    #         rv_planet = RV_from_params(t_val, params, ignore_mean=False, companion=True)
    #         ax1.plot(phi, rv_star, "+", markersize=12, markeredgewidth=3)
    #         ax2.plot(phi, rv_planet, "+", markersize=12, markeredgewidth=3)

    # Determine rv max amplitudes.
    a_1 = params["k1"] * (1 + params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    b_1 = params["k1"] * (1 - params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    amp1 = max([abs(a_1), abs(b_1)])

    a_2 = params["k2"] * (1 + params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    b_2 = params["k2"] * (1 - params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    amp2 = max([abs(a_2), abs(b_2)])

    # Adjust axis limits
    ax1.set_ylim(params["mean_val"] - (amp1 * 1.1), params["mean_val"] + (amp1 * 1.1))
    ax2.set_ylim(params["mean_val"] - (amp2 * 1.1), params["mean_val"] + (amp2 * 1.1))

    ax1.axhline(params["mean_val"], color="black", linestyle="-.", alpha=0.5)
    ax2.axhline(params["mean_val"], color="black", linestyle="-.", alpha=0.5)

    plt.legend(loc=0)
    plt.show()
    return 0


# Lots of duplication - could be improved
def RV_time_curve(params, cycle_fraction=1, ignore_mean=False, t_past=False, t_future=False, start_day=None):
    # type: (Dict[str, Union[str, float]], float, bool, Any, Any, Any) -> int
    """Plot RV phase curve centered on zero.

    Parameters
    ----------
    params: dict
        Parameters of system.
    cycle_fraction: float
        Fraction of phase space to plot. (Default=1)
    ignore_mean: bool
        Remove the contribution from the systems mean velocity. (Default=False)
    t_past: float, array-like
        Times of past observations in julian days.
    t_future: float, array-like
        Times of future observations  in julian days.
    start_day: str
        Day to being RV curve from in julian days. The Default=None sets the start time to ephem.now().

    Returns
    -------
    exit_status: bool
        Returns 0 if successful.

        Displays matplotlib figure.
    """
    if start_day is None:
        ephem_now = ephem.now()
        t_start = ephem.julian_date(ephem_now)
    else:
        t_start = start_day

    # Make curve from start of t_past
    if isinstance(t_past, float):
        obs_start = t_past
    elif t_past is not None:
        obs_start = np.min(t_past)
    else:
        obs_start = t_start
    # Specify 100 points per period
    num_cycles = ((t_start + params["period"] * cycle_fraction) - np.min([t_start, obs_start])) / params["period"]
    num_points = np.ceil(500 * num_cycles)
    if num_points > 10000:
        logging.debug("num points" = {}.format(num_points))
        raise ValueError("num_points is to large")

    t_space = np.linspace(min([t_start, obs_start]), t_start + params["period"] * cycle_fraction, num_points)

    host_rvs = RV_from_params(t_space, params, ignore_mean=ignore_mean)
    companion_rvs = RV_from_params(t_space, params, ignore_mean=ignore_mean, companion=True)

    fig = plt.figure(figsize=(10, 7))
    fig.subplots_adjust()
    ax1 = host_subplot(111)
    ax1.plot(t_space - t_start, host_rvs, label="Host", lw=2, color="k")

    start_dt = JulianDate(t_start).to_datetime()
    if (start_dt.hour == 0) and (start_dt.minute == 0) and (start_dt.second == 0):
        start_string = datetime.strftime(start_dt, "%Y-%m-%d")
    else:
        start_string = datetime.strftime(start_dt, "%Y-%m-%d %H:%M:%S")   # Issue with 00:00:01 not appearing
    ax1.set_xlabel("Days from {!s}".format(start_string))

    ax1.set_ylabel("Host RV (km/s)")

    ax2 = ax1.twinx()
    if params["true_mass_flag"]:
        ax2.plot(t_space - t_start, companion_rvs, '--', label="M_2 Companion", lw=2)
    else:
        ax2.plot(t_space - t_start, companion_rvs, '--', label="M_2sini Companion", lw=2)
    ax2.set_ylabel("Companion RV (km/s)")

    if 'name' in params.keys():
        if "companion" in params.keys():
            plt.title("Radial Velocity Curve for {} {}".format(params['name'].upper(), params['companion']))
        else:
            plt.title("Radial Velocity Curve for {}".format(params['name'].upper()))
    else:
        plt.title("Radial Velocity Curve")
    if t_past is not None:
        t_past = np.asarray(t_past)
        # for t_num, t_val in enumerate(t_past):
        rv_star = RV_from_params(t_past, params, ignore_mean=False)
        rv_planet = RV_from_params(t_past, params, ignore_mean=False, companion=True)
        ax1.plot(t_past - t_start, rv_star, ".", markersize=10, markeredgewidth=2, label="Host Obs")
        ax2.plot(t_past - t_start, rv_planet, "+", markersize=10, markeredgewidth=2, label="Comp Obs")

    if t_future:
        raise NotImplementedError("Adding future observations times is not implemented yet.")
    # if t_future:
    #     for t_num, t_val in enumerate(t_future):
    #         phi = ((t_val - params[4])/params[5] - 0.5) % 1 + 0.5
    #         rv_star = RV_from_params(t_val, params, ignore_mean=False)
    #         rv_planet = RV_from_params(t_val, params, ignore_mean=False, companion=True)
    #         ax1.plot(phi, rv_star, "+", markersize=12, markeredgewidth=3)
    #         ax2.plot(phi, rv_planet, "+", markersize=12, markeredgewidth=3)

    # Determine rv max amplitudes.
    a_1 = params["k1"] * (1 + params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    b_1 = params["k1"] * (1 - params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    amp1 = max([abs(a_1), abs(b_1)])

    a_2 = params["k2"] * (1 + params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    b_2 = params["k2"] * (1 - params["eccentricity"] * np.cos(params["omega"] * np.pi / 180))
    amp2 = max([abs(a_2), abs(b_2)])

    # Adjust axis limits
    ax1.set_ylim(params["mean_val"] - (amp1 * 1.1), params["mean_val"] + (amp1 * 1.1))
    ax2.set_ylim(params["mean_val"] - (amp2 * 1.1), params["mean_val"] + (amp2 * 1.1))

    ax1.axhline(params["mean_val"], color="black", linestyle="-.", alpha=0.5)
    ax2.axhline(params["mean_val"], color="black", linestyle="-.", alpha=0.5)

    plt.legend(loc=0)
    plt.show()
    return 0


if __name__ == '__main__':
    args = vars(_parser())
    debug = args.pop('debug')
    opts = {k: args[k] for k in args}

    if debug:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(message)s')
    else:
        logging.basicConfig(level=logging.WARNING,
                            format='%(asctime)s %(levelname)s %(message)s')

    main(**opts)
