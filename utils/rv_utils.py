import datetime
import numpy as np
from typing import Any, Dict
from utils.parse import parse_paramfile


from utils.parse import parse_paramfile

# TODO: Replace "Any" with numpy type hint when available


class RV(object):
    def __init__(self, semi_amp=0.0, period=0.0, ecc =0.0,
                 tau=0.0, gamma=0.0, omega=0.0, **other_params):
        self.semi_amp = semi_amp
        self.period = period
        self.ecc = ecc
        self.tau = tau
        self.gamma = gamma
        self.omega = omega
        self.params = self.param_dict()

        if other_params is not None:
            self.params.update(other_params)

    def __repr__(self):
        return "RV(semi_amp={1}, period={2}, ecc={3}, tau={4}, omega={5}, gamma={6}, params={7})".format(
            self.__class__, self.semi_amp, self.period, self.ecc, self.tau, self.omega, self.gamma, self.params)

    def param_dict(self):
        return {"k1": self.semi_amp, "period": self.period, "eccentricity": self.ecc,
                "tau": self.tau, "mean_val": self.gamma, "omega": self.omega}

    @classmethod
    def from_dict(cls, params):
        return cls(semi_amp=params["k1"], period=params["period"], ecc=params["eccentricity"],
                  tau=params["tau"], gamma=params["mean_val"], omega=params["omega"], other_params=params)

    @classmethod
    def from_file(cls, filename):
        """Parameters in key = val\n text file."""
        param_dict = parse_paramfile(filename)
        return cls.from_dict(param_dict)

    def rv_at_phase(self, phase):
        t = phase * self.period + self.tau
        return self.rv_at_times(t)

    def rv_at_times(self, t):
        """Evaluate RV at the provided times."""
        true_anomaly = self.true_anomaly(self.mean_anomaly(t, self.tau, self.period))
        return self.radial_velocity(self.gamma, self.semi_amp,
                                    true_anomaly, self.omega, self.ecc )

    def rv_full_phase(self, points=100):
        """Return RV curve evaluated one full phase."""
        phase = np.linspace(0, 2*np.pi, points)
        times = phase * self.params["period"] + self.params["t0"]
        return self.rv_at_times(times)

    # #######################################################
    # Functions for RV calculations
    # #######################################################
    @staticmethod
    def true_anomaly(ma, ecc, niterationmax=10000):
        # type: (Any, float, int) -> Any
        """Compute the true anomaly using the Newton-Raphson method.

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
        if not isinstance(ma, (int, float)):
            ea = ma
        else:
            ea = np.array([ma, ])

        if isinstance(ea, list):
            raise TypeError("Unsupported type 'list', input a numpy array or an int/float.")
        if len(ea) == 0:
            raise ValueError("A empty array was given.")

        # Initialise at ea0 = ma
        niteration = 0
        ea0 = ma

        while np.linalg.norm(ea - ea0, ord=1) > 1e-5 or niteration == 0:
            ea0 = ea

            ff = ea - ecc * np.sin(ea) - ma   # Function
            dff = 1 - ecc * np.cos(ea)        # Derivative

            # Use Newton method
            ea = ea0 - ff / dff

            # Increase iteration number; if above limit, break with exception.
            niteration += 1
            if niteration >= niterationmax:
                raise RuntimeError('Eccentric anomaly computation'
                                   'not converged.')

        # Compute true anomaly from eccentric anomaly
        return 2. * np.arctan2(np.sqrt(1. + ecc) * np.sin(ea / 2.),
                               np.sqrt(1. - ecc) * np.cos(ea / 2.))

    @staticmethod
    def mean_anomaly(times, t0, period):
        # type: (Any, float, float) -> Any
        """Calculate mean anomaly using period, tau and a time value.

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

    @staticmethod
    def radial_velocity(gamma, k, ta, omega, ecc):
        # type: (float, float, Any, float, float, float) -> Any
        """Radial velocity equation.

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
            Radial velocity values

        Notes
        -----
        RV = gamma + k *(np.cos(ta + omega) + ecc * np.cos(omega)).

        """
        # Calculate radial velocity of star
        return gamma + k * (np.cos(ta + omega) + ecc * np.cos(omega))


# RV calculation done in python (for when ajplanet is not available)
def rv_curve_py(times, gamma, k, omega, ecc, t0, period):
    # type: (Any, float, float, float, float, float, float) -> Any
    """Generate values for Radial Velocity curve.

    Parameters
    times: array-like, float
        Times to calculate radial velocity values(Julian days))
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
        Radial velocity values evaluated at the given times.

    """
    ma = RV.mean_anomaly(times, t0, period)
    ta = RV.true_anomaly(ma, ecc)
    rv = RV.radial_velocity(gamma, k, ta, omega, ecc)
    return rv


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

    param_list[2] = np.deg2rad(param_list[2])   # convert omega to radians

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
    sun_jupiter_mass = 1047.56  # Solar mass in jupiter masses
    m_host *= sun_jupiter_mass  # Convert to jupiter mass
    return -k_host * m_host / m_companion


def jd2datetime(jd, reduced=False):
    # type: (float) -> Any
    """Convert from a julian-date to a datetime object.

    Parameters
    ----------
    jd: float
        Julian date to calculate datetime for.
    reduced: bool
        Is input jd in Reduced JD format, (JD-2400000)

    Returns
    -------
    dt: datetime object
        Datetime of julian date.
    Inspiration from https://stackoverflow.com/questions/13943062/
    """
    if reduced:
        jd = jd + 2400000
    julian_epoch = datetime.datetime(2000, 1, 1, 12)  # noon (the epoch name is unrelated)
    j2000_jd = datetime.timedelta(2451545)            # julian epoch in julian dates

    dt = datetime.timedelta(jd) + julian_epoch - j2000_jd
    return dt


def datetime2jd(dt, reduced=False):
    # type: (Any) -> float
    """Convert from a datetime to a jd object.

    Test against pyehem.julian_date()

    Parameters
    ----------
    dt: datetime object
        Datetime for date to calculate jd.
    reduced: bool
        Return reduced JD, (JD-2400000)

    Returns
    -------
    jd: float
        Julian date time
    Inspiration from https://stackoverflow.com/questions/13943062/
    """
    julian_epoch = datetime.datetime(2000, 1, 1, 12)  # noon (the epoch name is unrelated)
    j2000_jd = datetime.timedelta(2451545)            # julian epoch in julian dates

    jd = dt - julian_epoch + j2000_jd

    jd = jd.total_seconds() / (24 * 60 * 60)  # Turn timedelta into a float
    if reduced:
        jd = jd - 2400000
    return jd
