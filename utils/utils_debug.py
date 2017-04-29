import inspect

def pv(name):
    # record=inspect.getouterframes(inspect.currentframe())[1]
    # frame=record[0]
    frame = inspect.currentframe().f_back
    val = eval(name, frame.f_globals, frame.f_locals)
    # print('{0}: {1}'.format(name, val))
    return '{0}: {1}'.format(name, val)
