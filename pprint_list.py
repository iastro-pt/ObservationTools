"""
Pretty-print a Python list
https://stackoverflow.com/a/36085705
"""

import math

def pprint(obj, cols=4, columnwise=True, gap=4):
    """
    Print the given list in evenly-spaced columns.

    Parameters
    ----------
    obj : list
        The list to be printed.
    cols : int
        The number of columns in which the list should be printed.
    columnwise : bool, default=True
        If True, the items in the list will be printed column-wise.
        If False the items in the list will be printed row-wise.
    gap : int
        The number of spaces that should separate the longest column
        item/s from the next column. This is the effective spacing
        between columns based on the maximum len() of the list items.
    """

    sobj = [str(item) for item in obj]
    if cols > len(sobj): 
        cols = len(sobj)

    max_len = max([len(item) for item in sobj])
    if columnwise: 
        cols = int(math.ceil(float(len(sobj)) / float(cols)))
    
    plist = [sobj[i: i+cols] for i in range(0, len(sobj), cols)]
    
    if columnwise:
        if not len(plist[-1]) == cols:
            plist[-1].extend(['']*(len(sobj) - len(plist[-1])))
        plist = zip(*plist)
    
    printer = '\n'.join([
        ''.join([c.ljust(max_len + gap) for c in p])
        for p in plist])

    print (printer)