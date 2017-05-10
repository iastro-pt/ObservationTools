"""Functions for parsing the parameter files."""
import os
from utils.utils_debug import pv
from logging import debug
from typing import List, Dict


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
    try:
        return [float(val) for val in list_str]
    except ValueError as e:
        # Can't turn into floats.
        return [val.strip() for val in list_str]
