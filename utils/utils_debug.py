"""Function for Debuging."""
import inspect


def pv(name):
    # type: (str) -> str
    """Analysis an expresion 'expresion : evaulation'.

    Used to help debuging values.
    """
    if "__" in name:
        raise ValueError("Double underscores not allowed for saftey reasons.")
    frame = inspect.currentframe().f_back
    val = eval(name, frame.f_globals, frame.f_locals)
    return '{0}: {1}'.format(name, val)
