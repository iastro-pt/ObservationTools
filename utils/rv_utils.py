import copy
import datetime
import logging
from typing import Any

import numpy as np
from astropy.constants import M_jup, M_sun

from utils.parse import parse_obslist, parse_paramfile


# TODO: Replace "Any" with numpy type hint when available


class RV(object):
    """

    Notes
    -----
    1) Omega should be given in degrees. This function converts it to radians.
    2) The units of mean_val and k1, k2 should be the same e.g.both km / s
    3) Tau should be the julian date, and the period given in days.
    """
    def __init__(self, semi_amp=0.0, period=0.0, ecc=0.0,
                 tau=0.0, gamma=0.0, omega=0.0, **other_params):
        self.semi_amp = semi_amp
        self.period = period
        self.ecc = ecc
        self.tau = tau
        self.gamma = gamma
        self.omega = omega
        self._params = self.orbit_dict()

        self.ignore_mean = other_params.get("ignore_mean", False)
        if other_params is not None:
            self._params.update(other_params)

    def __repr__(self):
        other_params = ""
        for key, value in self._params.items():
            if key not in ["k1", "eccentricity", "period", "mean_val", "tau", "omega"]:
                other_params += ", {}={}".format(key, value)

        return "RV(semi_amp={1}, period={2}, ecc={3}, tau={4}, omega={5}, gamma={6}{7})".format(
            self.__class__, self.semi_amp, self.period, self.ecc, self.tau, self.omega, self.gamma, other_params)
    
    def orbit_dict(self):
        return {"k1": self.semi_amp, "period": self.period, "eccentricity": self.ecc,
                "tau": self.tau, "mean_val": self.gamma, "omega": self.omega}

    def to_dict(self):
        self._params.update(self.orbit_dict())
        self._params.update({"ignore_mean": self.ignore_mean})
        return self._params

    @classmethod
    def from_dict(cls, params):
        other_params = params.copy()
        for par in ["k1", "eccentricity", "period", "mean_val", "tau", "omega"]:
            other_params.pop(par)
        return cls(semi_amp=params["k1"], period=params["period"], ecc=params["eccentricity"],
                   tau=params["tau"], gamma=params["mean_val"], omega=params["omega"], **other_params)

    @classmethod
    def from_file(cls, filename):
        """Parameters in key = val\n text file."""
        param_dict = parse_paramfile(filename)
        return cls.from_dict(param_dict)

    def create_companion(self, mass_ratio=None):
        """Create companion RV object.

        It has 3 ways to determine the amplitude of the companion in order of priority:
            1: using a mass_ratio m1/m2 passed into the method.
            2: If "k2" is already a parameter use that.
            3: Calculate the mass ratio from m1 and m2. (In same units)

        Switches k1 and k2 and m1 an m2 parameters. (m1 refers to self, while m2 the other body in orbit.)

        Inputs
        ------
        mass_ratio: float
            m_1 / m_2

        Returns
        ------
        companion: RV
            Rv object for the companion.
        """
        params = copy.copy(self._params)
        if mass_ratio is not None:
            k2 = -params["k1"] * mass_ratio
            params["ratio_flag"] = True
        elif params.get("k2") is not None:
            k2 = params.get("k2")
            params["k2_flag"] = True
        else:
            # Make from masses in parameters
            M1 = params.get("m1")
            M2 = params.get("m2")
            if (M1 is None) or (M2 is None):
                print("params =", params)
                raise ValueError("A mass parameter (m1 or m2) was not provided")
            else:
                # M1_jup = M1 * (M_sun / M_jup).value
                mass_ratio = M1 / M2  # assuming m1 and m2 have same units
            k2 = -params["k1"] * mass_ratio
            params["m1"], params["m2"] = params.get("m2"), params.get("m1")

        params["k2"], params["k2_old"] = k2, params.get("k2")

        params["k1"], params["k2"] = params["k2"], params["k1"]  # Switch to Companion
        return RV.from_dict(params)

    @property
    def ignore_mean(self):
        try:
            return self._ignore_mean
        except AttributeError:
            return False
    @ignore_mean.setter
    def ignore_mean(self, value=None):
        if value is None:
            val = self._params.get("ignore_mean", False)
        self._ignore_mean = value

    def rv_at_phase(self, phase):
        t = phase * self.period + self.tau
        return self.rv_at_times(t)

    def rv_at_times(self, t):
        """Evaluate RV at the provided times."""
        true_anomaly = self.true_anomaly(self.mean_anomaly(t, self.tau, self.period), self.ecc)
        return self.radial_velocity(self.gamma, self.semi_amp,
                                    true_anomaly, self.omega, self.ecc)

    def rv_full_phase(self, center=0, points=100):
        """Return RV curve evaluated one full phase."""
        phase = np.linspace(0, 1, points) + center
        return self.rv_at_phase(phase)

    def max_amp(self):
        amp_1 = self.semi_amp * (1 + self.ecc * np.cos(self.omega * np.pi / 180))
        amp_2 = self.semi_amp * (1 - self.ecc * np.cos(self.omega * np.pi / 180))
        return np.max([np.abs(amp_1), np.abs(amp_2)])

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
            ea = np.array([ma, ], dtype=np.float)

        if isinstance(ea, list):
            raise TypeError("Unsupported type 'list', input a numpy array or an int/float.")
        if len(ea) == 0:
            raise ValueError("A empty array was given.")

        # Initialise at ea0 = ma
        niteration = 0
        ea0 = ma

        while np.linalg.norm(ea - ea0, ord=1) > 1e-5 or niteration == 0:
            ea0 = ea

            ff = ea - ecc * np.sin(ea) - ma  # Function
            dff = 1 - ecc * np.cos(ea)  # Derivative

            # Use Newton method
            ea = ea0 - ff / dff

            # Increase iteration number; if above limit, break with exception.
            niteration += 1
            if niteration >= niterationmax:
                raise RuntimeError('Eccentric anomaly computation '
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

    def __eq__(self, other):
        # all properties the same
        checks = [self.semi_amp == other.semi_amp,
                  self.period == other.period,
                  self.omega == other.omega,
                  self.tau == other.tau,
                  self.ecc == other.ecc,
                  self.gamma == other.gamma]

        return all(checks)

    def __ne__(self, other):
        return not self.__eq__(other)


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


class JulianDate(object):
    """Handle julian dates."""
    julian_epoch_dt = datetime.datetime(2000, 1, 1, 12)  # noon
    julian_epoch_jd = datetime.timedelta(2451545)  # julian epoch in julian dates
    reduce_jd = 2400000
    strformat = "%Y-%m-%d %H:%M:%S"

    def __init__(self, jd, reduced=False):
        self.jd = jd
        self.reduced = reduced

    @classmethod
    def from_datetime(cls, dt, reduced=False):
        # type: (Any, bool) -> JulianDate
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
        jd: JulianDate
            JulianDate object.
        Inspiration from https://stackoverflow.com/questions/13943062/
        """
        jd = dt - cls.julian_epoch_dt + cls.julian_epoch_jd
        jd = float(jd.total_seconds()) / (24 * 60 * 60)  # Turn timedelta into a float
        if reduced:
            jd -= cls.reduce_jd
        return cls(jd, reduced=reduced)

    def to_datetime(self):
        # type: None -> datetime.datetime
        """ Return JulianDate as a datetime.datetime object.

        Returns
        -------
        dt: datetime object
            Datetime of julian date.
        Inspiration from https://stackoverflow.com/questions/13943062/
        """
        print("self.jd", self.jd)
        if self.reduced:
            _jd = self.jd + self.reduce_jd
        else:
            _jd = self.jd
        print("_jd", _jd)
        dt = datetime.timedelta(_jd) + self.julian_epoch_dt - self.julian_epoch_jd
        return dt

    @classmethod
    def from_str(cls, time_str, format=strformat):
        """Return JulianDate from a time string.

        Inputs
        ------
        time_str: str
        format: str
            Format of time string.

        Returns
        -------
        dt: datetime object
            Datetime of julian date.
        """
        if format is None:
            format = cls.strformat
        dt = datetime.datetime.strptime(time_str, format)
        return cls.from_datetime(dt)

    def to_str(self, format=None):
        """Return date string from a JulianDate.

        Input
        -----
        format: str
             String datetime format.

        Returns
        -------
        datesting: str
            String with date.
        """
        if format is None:
            format = self.strformat
        dt = self.to_datetime()
        return dt.strftime(format)

    def reduce(self):
        if not self.reduced:
            self.jd -= self.reduce_jd
            self.reduced = True


def strtimes2jd(obs_times, reduced=False, format=None):
    # type: (Union[str, List[str]], bool, Union[None, str]) -> List[float]
    """Convenience function for convert str times to reduced JD.
    If reduced=True returns JD-2400000
    """
    reduce_value = 2400000 if reduced else 0

    if obs_times is not None:
        print("obs times", obs_times)

        if isinstance(obs_times, str):
            if reduced:
                jds = JulianDate.from_str(obs_times, format)
                jds.reduce()
                jds = jds.jd
            else:
                jds = JulianDate.from_str(obs_times, format).jd
        elif isinstance(obs_times, (list, tuple)):
            if reduced:
                jds = []
                for obs in obs_times:
                    jd = JulianDate.from_str(obs, format)
                    jd.reduce()
                    jds.append(jd.jd)
                    # jds = [JulianDate.from_str(obs).reduce().jd for obs in obs_times]
            else:
                jds = [JulianDate.from_str(obs, format).jd for obs in obs_times]
        print("obs jd times", jds)
        return jds
    else:
        return None


def join_times(obs_times=None, obs_list=None):
    # type: (List[str], str) -> List[str]
    """Combine observation dates and turn to jd.

    Parameters
    ----------
    obs_times: list of str or None
        List of dates entered at command line.
    obs_list: str or None
        Filename to observation list.

    Returns
    -------
    obs_times: list of str
        Combined list of date strings.

    """
    if obs_times is None:
        obs_times = []

    if obs_list is None:
        obs_list = []
    else:
        obs_list = parse_obslist(obs_list)

    logging.debug("obs list = {}", format(obs_list))
    obs_times = obs_times + obs_list

    logging.debug("obs_times = {}".format(obs_times))
    if not obs_times:  # An empty list is "Falsely"
        return None
    else:
        return obs_times


def prepare_mass_params(params, only_msini=True):
    """Update parameter dictionary to set m1 and m2 if not given."""
    if params.get("m1") is None:
        params["m1"] = params["m_sun"]  # solar mass

    # Convert m_sun to jupyter masses
    params["m1"] = params["m1"] * M_sun / M_jup  # jupyter mass

    if params.get("m2") is None:
        params["m2"] = params["m_sini"] if only_msini else params["m_true"]
        # should give a key error if the correct mass not given

    params["msini_flag"] = only_msini

    params["k2_flag"] = False if params.get("k2") is None else True

    return params


def check_core_parameters(params):
    """Test core rv parameters."""
    for key in ["name", "k1", "eccentricity", "omega", "tau", "period"]:
        if key not in params.keys():
            raise ValueError("A core parameter was not provided in the param file, '{}'".format(key))

    if "mean_val" not in params.keys():
        logging.info("mean_val parameter was not provided so set to 0 km/s")
        params["mean_val"] = 0.0
    elif params["mean_val"] == "":
        logging.info("mean_val parameter was blank so set to 0 km/s")
        params["mean_val"] = 0.0

    return params


def generate_companion_label(companion):
    msini_flag = companion._params.get("msini_flag", False)
    k2_flag = companion._params.get("k2_flag", False)
    ratio_flag = companion._params.get("ratio_flag", False)

    if not msini_flag and not k2_flag and not ratio_flag:
        label = "M2 Companion"
    elif msini_flag and not k2_flag and not ratio_flag:
        label = "M2sini Companion"
    elif k2_flag and not ratio_flag:
        label = "Given k2 Companion"
    else:
        label = "Mass ratio Companion"

    return label
