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
from utils.rv_utils import RV, JulianDate, prepare_mass_params
from utils.rv_utils import strtimes2jd, join_times, check_core_parameters

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
    parser.add_argument("--save_only", help="Only save the figure, do not show it.", action="store_true")
    parser.add_argument("--debug", help="Turning on debug output", action='store_true', default=False)
    return parser.parse_args()


def main(params, mode="phase", obs_times=None, obs_list=None, date=None,
         save_only=False):  # obs_times=None, mode='phase', rv_diff=None
    # type: (Dict[str, Union[str, float]], str, List[str], str, str, bool) -> None
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
    only_msini = False
    # Load in params and store as a dictionary
    parameters = parse_paramfile(params)

    # Convert msun to jupyter masses
    parameters["m1"] = parameters["m1"] * M_sun / M_jup

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
    obs_jd = strtimes2jd(obs_times)

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
                    true_mass_flag = True  # parameter to indicate if the true mass was used or not
                else:
                    mass_used = parameters["msini"]
                    true_mass_flag = False
            else:
                mass_used = parameters["msini"]
                true_mass_flag = False

            parameters['k2'] = companion_amplitude(parameters['k1'],
                                                   parameters['m_star'],
                                                   mass_used)
            parameters["true_mass_flag"] = true_mass_flag  # True if true mass used

    host_orbit = RV.from_dict(parameters)
    companion_orbit = host_orbit.create_companion()

    if mode == "phase":
        fig = binary_phase_curve(host_orbit, companion_orbit, t_past=obs_jd)
    elif mode == "time":
        if date is not None:
            print("Date doing into ehem.julian_date", date)
            date = ephem.julian_date(date)
        fig = binary_time_curve(host_orbit, companion_orbit, t_past=obs_jd, start_day=date)
    else:
        raise NotImplementedError("Other modes are not Implemented yet.")
    if not save_only:
        fig.show()

    return fig


def binary_phase_curve(host, companion, cycle_fraction=1, ignore_mean=False, t_past=False, t_future=False):
    # type: (RV, RV, float, bool, Any, Any) -> int
    """Plot RV phase curve centered on zero.

    Parameters
    ----------
    host: RV
        RV class for systems host.
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
    fig: matplotlib.Figure
        Returns figure object.

    """
    host.ignore_mean = ignore_mean
    companion.ignore_mean = ignore_mean

    phase = np.linspace(-0.5, 0.5, 100) * cycle_fraction
    t = host.tau + phase * host.period

    host_rvs = host.rv_at_phase(phase)
    companion_rvs = companion.rv_at_phase(phase)

    fig = plt.figure(figsize=(10, 7))
    fig.subplots_adjust()
    ax1 = host_subplot(111)
    ax1.plot(phase, host_rvs, label="Host", lw=2, color="k")
    ax1.set_xlabel("Orbital Phase")
    ax1.set_ylabel("Host RV (km/s)")

    ax2 = ax1.twinx()
    companion_label = generate_companion_label(companion)
    ax2.plot(phase, companion_rvs, '--', label=companion_label, lw=2)
    ax2.set_ylabel("Companion RV (km/s)")

    if 'name' in host._params.keys():
        if "companion" in host._params.keys():
            plt.title("RV Phase Curve for {} {}".format(host._params['name'].upper(), host._params['companion']))
        else:
            plt.title("RV Phase Curve for {}".format(host._params['name'].upper()))
    else:
        plt.title("RV Phase Curve")

    if t_past:
        for t_num, t_val in enumerate(t_past):
            phi = ((t_val - host.tau) / host.period + 0.5) % 1 - 0.5
            rv_star = host.rv_at_phase(phi)
            rv_planet = companion.rv_at_phase(phi)
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
    host_delta_y = host.max_amp() * 1.1
    comp_delta_y = companion.max_amp() * 1.1

    ax1.set_ylim(host.gamma - host_delta_y, host.gamma + host_delta_y)
    ax2.set_ylim(companion.gamma - comp_delta_y, companion.gamma + comp_delta_y)

    ax1.axhline(host.gamma, color="black", linestyle="-.", alpha=0.5)
    ax2.axhline(companion.gamma, color="black", linestyle="-.", alpha=0.5)

    plt.legend(loc=0)
    return fig


# Lots of duplication - could be improved


def binary_time_curve(host, companion, cycle_fraction=1, ignore_mean=False, t_past=False, t_future=False,
                      start_day=None):
    # type: (RV, RV, float, bool, Any, Any, Any) -> int
    """Plot RV phase curve centered on zero.

    Parameters
    ----------
    host: RV
        RV class for system.
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
    host.ignore_mean = ignore_mean
    companion.ignore_mean = ignore_mean

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
    num_cycles = ((t_start + host.period * cycle_fraction) - np.min([t_start, obs_start])) / host.period
    num_points = np.ceil(500 * num_cycles)
    if num_points > 10000:
        logging.debug("num points = {}".format(num_points))
        raise ValueError("num_points is to large")

    t_space = np.linspace(min([t_start, obs_start]), t_start + host.period * cycle_fraction, num_points)

    # host_rvs = RV_from_params(t_space, params, ignore_mean=ignore_mean)
    host_rvs = host.rv_at_times(t_space)
    companion_rvs = companion.rv_at_times(t_space)

    fig = plt.figure(figsize=(10, 7))
    fig.subplots_adjust()
    ax1 = host_subplot(111)
    ax1.plot(t_space - t_start, host_rvs, label="Host", lw=2, color="k")

    start_dt = JulianDate(t_start).to_datetime()
    if (start_dt.hour == 0) and (start_dt.minute == 0) and (start_dt.second == 0):
        start_string = datetime.strftime(start_dt, "%Y-%m-%d")
    else:
        start_string = datetime.strftime(start_dt, "%Y-%m-%d %H:%M:%S")  # Issue with 00:00:01 not appearing
    ax1.set_xlabel("Days from {!s}".format(start_string))

    ax1.set_ylabel("Host RV (km/s)")

    ax2 = ax1.twinx()
    companion_label = generate_companion_label(companion)
    ax2.plot(t_space - t_start, companion_rvs, '--', label=companion_label, lw=2)
    ax2.set_ylabel("Companion RV (km/s)")

    if 'name' in host._params.keys():
        if "companion" in host._params.keys():
            plt.title("Radial Velocity Curve for {} {}".format(host._params['name'].upper(), host._params['companion']))
        else:
            plt.title("Radial Velocity Curve for {}".format(host._params['name'].upper()))
    else:
        plt.title("Radial Velocity Curve")
    if t_past is not None:
        t_past = np.asarray(t_past)
        rv_star = host.rv_at_times(t_past)
        rv_planet = companion.rv_at_times(t_past)
        ax1.plot(t_past - t_start, rv_star, "b.", markersize=10, markeredgewidth=2, label="Host Obs")
        ax2.plot(t_past - t_start, rv_planet, "r+", markersize=10, markeredgewidth=2, label="Comp Obs")

    if t_future:
        t_past = np.asarray(t_future)
        rv_star = host.rv_at_times(t_future)
        rv_planet = companion.rv_at_times(t_future)
        ax1.plot(t_future - t_start, rv_star, "ko", markersize=10, markeredgewidth=2, label="Host Obs")
        ax2.plot(t_future - t_start, rv_planet, "g*", markersize=10, markeredgewidth=2, label="Comp Obs")
    # raise NotImplementedError("Adding future observations times is not implemented yet.")

    # Determine rv max amplitudes.
    amp1 = host.max_amp()
    amp2 = companion.max_amp()

    # Adjust axis limits
    ax1.set_ylim(host.gamma - (amp1 * 1.1), host.gamma + (amp1 * 1.1))
    ax2.set_ylim(companion.gamma - (amp2 * 1.1), companion.gamma + (amp2 * 1.1))

    ax1.axhline(host.gamma, color="black", linestyle="-.", alpha=0.5)
    ax2.axhline(companion.gamma, color="black", linestyle="-.", alpha=0.5)

    plt.legend(loc=0)

    return fig


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

    fig = main(**opts)
