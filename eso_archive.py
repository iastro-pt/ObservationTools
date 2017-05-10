#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from astroquery.eso import Eso
import argparse


def _parser():
    eso = Eso()
    instruments = eso.list_instruments()
    parser = argparse.ArgumentParser(description='Search and download from the ESO archive.')
    parser.add_argument('username', help='ESO username')
    parser.add_argument('--store_password', help='Store ESO password', default=False, action='store_true')
    parser.add_argument('--object', help='Object, e.g. "Aldebaran"', default=False)
    parser.add_argument('--instrument', help='ESO instruments', choices=instruments, type=str.lower)
    return parser.parse_args()


if __name__ == '__main__':
    args = _parser()
    # Login
    eso = Eso()
    eso.login(args.username, store_password=args.store_password)

    target = args.object
    instrument = args.instrument

    if instrument is None:
        raise ValueError("Add an instrument to query.")
    if target is False:
        raise ValueError("Include a target to query.")

    table = eso.query_instrument(instrument, column_filters={'target': target})

    if table is None:
        print("The query 'instrument={}, target={}'. Returned no results.".format(instrument, target))
    else:
        table.pprint()


    # data_files = eso.retrieve_data(table['DP.ID'])
    # print data_files
