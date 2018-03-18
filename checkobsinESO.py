# from __future__ import print_function
import numpy as np
from astroquery.eso import Eso
import argparse
import warnings
from pprint_list import pprint

def _parser():
    parser = argparse.ArgumentParser(
        description='Check if targets have observations in the ESO archive.')
    parser.add_argument('targets', nargs='+',
                        help='which targets to search for, '
                             'e.g., HD20010 or HD20010 HD41248')
    return parser.parse_args()


if __name__ == '__main__':
    args = _parser()

    insts = Eso.list_instruments()
    print ('Checking observations in ESO Archive for instruments')
    pprint (insts, cols=6)


    for star in args.targets:
        print ('\n'+star+':')
        found = False

        # Search in all the instruments in ESO Archive  
        for inst in Eso.list_instruments():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                x = Eso.query_instrument(inst, 
                        column_filters={'target': star,
                                        'dp_cat': 'SCIENCE',
                                        'box':'00 05 00'})
            try:
                if len(x)>0:
                    print (' - %s (%d)' % (inst.upper(), len(x)))
                    found = True
            except:
                continue
        
        if not found:
            print ('No observations found')
                
