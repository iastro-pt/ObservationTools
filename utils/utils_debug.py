"""Function for Debuging."""
import inspect


def pv(name):
    # type: (str) -> str
    """Analysis an expresion 'expresion : evaulation'.

    Used to help debuging values.
    """
    frame = inspect.currentframe().f_back
    val = eval(name, frame.f_globals, frame.f_locals)
    return '{0}: {1}'.format(name, val)
