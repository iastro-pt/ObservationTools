
# Radial Velocity calculations:

# Goals:
# To calcualte when the radial velocity is different by a certian value.
# Plot radial velocity phase curves. Indicating obtained measurement locations

import argparse
import numpy as np

try:
    from ajplanet import pl_rv_array
    use_ajplanet = True
except:
    use_ajplanet = False

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot


def _parser():
    parser = argparse.ArgumentParser(description='Radial velocity plotting')
    parser.add_argument('params', help='RV parameters filename')
    # parser.add_argument('-d', '--date', default='today',
    #                     help='Date in format YYYY-MM-DD. Default is today.')
    # #  prediciting separated lines
    # parser.add_argument('-r', '--rv_diff', help='RV difference threshold to find')
    # parser.add_argument('-o', '--obs_times',  help='Times of previous observations', nargs='+')
    # parser.add_argument('-f', '--files', help='Params and obs-times are file'
    #                    ' names to open', action='store_true')
    parser.add_argument('-m', '--mode', help='Display mode '
                        ' e.g. phase or time plot. Default="phase"',
                        choices=['time', 'phase'], default='phase')
    return parser.parse_args()


def companion_amplitude(k_host, m_host, m_companion):
    """ Calcualte the companion RV maximum amplitude.

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


def main(params, mode="phase"):  # obs_times=None, mode='phase', rv_diff=None
    """ Do main stuff.

    Parameters
    ----------
    params: str
        Filenamefor text file containing the rv parameters. Format of 'param = value\n'.
    mode: str
        Mode for script to use. Phase, time, future
    """

    # Load in params and store as a dictionary
    parameters = dict()
    with open(params, 'r') as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                par, val = line.lower().split('=')
                parameters[par.strip()] = val.strip()

    # Turn most parameters to floats.
    for key in parameters.keys():
        if key in ['mean_val', 'k1', 'omega', 'eccentricity', 'tau', 'period',
                   'm_star', 'msini ', 'm_true']:
            parameters[key] = float(parameters[key])

    # Calculate companion semi-major axis
    if 'k2' in parameters.keys():
        pass
    else:
        if 'm_true' in parameters.keys():
            # Use true mass if given
            parameters['k2'] = companion_amplitude(parameters['k1'],
                                                   parameters['m_star'],
                                                   parameters['m_true'])
        else:
            parameters['k2'] = companion_amplitude(parameters['k1'],
                                                   parameters['m_star'],
                                                   parameters['msini'])

    # print(parameters)

    if mode == "phase":
        RV_phase_curve(parameters)
    else:
        raise NotImplemented


def RV_phase_curve(params, cycle_fraction=1, ignore_mean=False, t_past=False, t_future=False):
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
    None:
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


# #######################################################
# Functions for RV calculations
# #######################################################
def true_anomaly(ma, ecc, niterationmax=10000):
    """ Compute the true anomaly using the Newton-Raphson method.

    Parameters
    ----------
    ma: array-like
        Mean anomaly.
    ecc: float
        Orbital eccentricity.
    niterationmax: int
        Maximum number of iterations for N-R method.

    Returns
    -------
    ta: array-like
        True anomaly

    Notes
    -----
    Adapted from Rodrigo Diaz.
    """

    if not isinstance(ma, float):
        ea = ma
    else:
        ea = np.array([ma, ])

    # Initialise at ea0 = ma
    niteration = 0
    ea0 = ma

    while np.linalg.norm(ea - ea0, ord=1) > 1e-5 or niteration == 0:
        ea0 = ea

        ff = ea - ecc*np.sin(ea) - ma   # Function
        dff = 1 - ecc*np.cos(ea)        # Derivative

        # Use Newton method
        ea = ea0 - ff / dff

        # Increase iteration number; if above limit, break with exception.
        niteration += 1
        if niteration >= niterationmax:
            raise RuntimeError('Eccentric anomaly computation'
                               'not converged.')

    # Compute true anomaly from eccentric anomaly
    return 2. * np.arctan2(np.sqrt(1. + ecc) * np.sin(ea/2.),
                           np.sqrt(1. - ecc) * np.cos(ea/2.))


def mean_anomaly(times, t0, period):
    """ Calculate mean anomaly using period, tau and a time value

    Parameters
    ----------
    times: array-like
        Times to compute mean anomaly.
    t0: float
        Time of periastron passage. (Julian days)
    period: float
        Period of orbit.

    Returns
    -------
    ma: array-like
        Mean anomaly.
    """
    if not isinstance(times, (int, float)):
        times = times
    else:
        times = np.array(times)
    return 2 * np.pi * (times - t0) / period


def radial_velocity(gamma, k, ta, omega, ecc):
    """ Radial velocity equation.
    Parameters
    ----------
    gamma: float
        Mean RV motion of system.
    k: float
        RV amplitude.
    ta: array-like, float
        True anomaly.
    omega: float
        Argument of periastron. (radians)
    ecc: float
        Eccentricity of orbit.

    Returns
    -------
    RV: array-like, float
        Radial veloctity values

    Notes
    -----
    RV = gamma + k *(np.cos(ta + omega) + ecc * np.cos(omega))
    """
    # Calculate radial velocity of star
    return gamma + k * (np.cos(ta + omega) + ecc * np.cos(omega))


# RV calculation done in python (for when ajplanet is not available)
def rv_curve_py(times, gamma, k, omega, ecc, t0, period):
    """Generate values for Radial Velocity curve.

    Parameters
    times: array-like, float
        Times to calcualte radial velocity values(Julian days))
    gamma: float
        Mean RV offset value.
    k: float
        Radial velocity amplitude.
    omega: float
        Argument of periastron (radians).
    ecc: float
        Eccentricity.
    t0: float
        Time of periastron passage.
    period: float
        Period of orbit.

    Returns
    -------
    rv: array-like
        Radial velocity values evaulated at the given times.

    """
    ma = mean_anomaly(times, t0, period)
    ta = true_anomaly(ma, ecc)
    rv = radial_velocity(gamma, k, ta, omega, ecc)
    return rv


def RV_from_params(t, params, ignore_mean=False, companion=False):
    """ Get radial velocity values at given times using the orbital parameters.

    Parameters
    ----------
    t: array-like
        The times at which to calculate the RV.
    params: dict
        Orbtial parameters required for RV.(mean_val, k1, omega, e, tau, period, optinal [k2]).
    ignore_mean: bool
        Ignore the average radial velocity motion of the system
    companion: bool
        Calculate RV for companion instead.

    Notes:
    1) Omega should be given in degrees. This function converts it to radians.
    2) The units of mean_val and k1,k2 should be the same e.g. both km/s

    Returns:
    rvs: array-like
        The radial velocity values evaluated at the provided times.

    """
    if not isinstance(t, np.ndarray):
        t = np.array(t)

    # Select needed entries from dict to calcualte rv
    if companion:
        param_list = [params["mean_val"], params["k2"], params["omega"],
                      params["eccentricity"], params["tau"], params["period"]]
    else:
        param_list = [params["mean_val"], params["k1"], params["omega"],
                      params["eccentricity"], params["tau"], params["period"]]

    param_list[2] = np.deg2rad(param_list[2])   # convert omega to radians

    if not ignore_mean:
        # Set the mean rv veloctiy to zero
        param_list[0] = 0

    if use_ajplanet:   # Use ajplanet if available
        rvs = pl_rv_array(t, *param_list[:])  # *unpacks parameters from list
    else:              # Use python version
        rvs = rv_curve_py(t, *param_list[:])  # *unpacks parameters from list

    return rvs


if __name__ == '__main__':

    args = vars(_parser())
    # star_name = args.pop('star_name')
    opts = {k: args[k] for k in args}

    main(**opts)
