import numpy as np
from astropy.constants import M_sun, M_jup

from utils.rv_utils import rv_curve_py


def RV_from_params(t, params, ignore_mean=False, companion=False):
    # type: (Any, Dict[str, Union[str, float]], bool, bool) -> Any
    """Get radial velocity values at given times using the orbital parameters.

    Parameters
    ----------
    t: array-like
        The times at which to calculate the RV.
    params: dict
        Orbital parameters required for RV.(mean_val, k1, omega, e, tau, period, optimal [k2]).
    ignore_mean: bool
        Ignore the average radial velocity motion of the system. Default=False.
    companion: bool
        Calculate RV for companion instead. Default=False. (k2 required)

    Returns
    -------
    rvs: array-like
        The radial velocity values evaluated at the provided times.

    Notes
    -----
    1) Omega should be given in degrees. This function converts it to radians.
    2) The units of mean_val and k1, k2 should be the same e.g. both km/s
    3) Tau should be the julian date, and the period given in days.

    """
    if not isinstance(t, np.ndarray):
        t = np.array(t)

    # Select needed entries from dict to calculate rv
    if companion:
        param_list = [params["mean_val"], params["k2"], params["omega"],
                      params["eccentricity"], params["tau"], params["period"]]
    else:
        param_list = [params["mean_val"], params["k1"], params["omega"],
                      params["eccentricity"], params["tau"], params["period"]]

    param_list[2] = np.deg2rad(param_list[2])  # convert omega to radians

    for param in param_list:
        if isinstance(param, str):
            raise TypeError("One of the params was not converted to float, {}. Check the parameter file.".format(param))
    if ignore_mean:
        # Set the mean rv to zero
        param_list[0] = 0

    # if use_ajplanet:   # Use ajplanet if available
    #     rvs = pl_rv_array(t, *param_list[:])  # *unpacks parameters from list
    # else:              # Use python version
    rvs = rv_curve_py(t, *param_list[:])  # *unpacks parameters from list

    return rvs


def companion_amplitude(k_host, m_host, m_companion):
    # type: (float, float, float) -> float
    """Calculate the companion RV maximum amplitude.

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
    sun_jupiter_mass = (M_sun / M_jup).value  # Solar mass in jupiter masses
    m_host *= sun_jupiter_mass  # Convert to jupiter mass
    return -k_host * m_host / m_companion
