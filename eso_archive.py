#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import os
import argparse
from shutil import copyfile
from astroquery.eso import Eso


def _parser():
    eso = Eso()
    instruments = eso.list_instruments()
    parser = argparse.ArgumentParser(description='Search and download from the ESO archive.')
    parser.add_argument('-u', '--username', help='ESO username', default=False)
    parser.add_argument('-d', '--download', help='Download data', default=False, action='store_true')
    parser.add_argument('-n', '--number', help='How many data files will be downloaded (-1=all)', default=1, type=int)
    parser.add_argument('--store_password', help='Store ESO password', default=False, action='store_true')
    parser.add_argument('-o', '--object', help='Object, e.g. "Aldebaran"', default=False)
    parser.add_argument('-i', '--instrument', help='ESO instruments', choices=instruments, type=str.lower)
    return parser.parse_args()


def move_file(pwd, data_file):
    fname = data_file.rpartition('/')[-1]
    dst = '{}/{}'.format(pwd, fname)
    copyfile(data_file, dst)
    print('Downloaded: {}'.format(fname))


if __name__ == '__main__':
    args = _parser()
    # Login
    eso = Eso()
    username = args.username
    if args.download:
        if not username:
            username = raw_input('Download requested, please provide ESO username: ')
        try:
            eso.login(username, store_password=args.store_password)
        except:
            raise SystemExit('Invalid username or password for ESO')

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

    if args.download:
        pwd = os.getcwd()
        if args.number == -1:
            data_files = eso.retrieve_data(table['DP.ID'])
        else:
            data_files = eso.retrieve_data(table['DP.ID'][:args.number])
        if isinstance(data_files, list):
            for data_file in data_files:
                move_file(pwd, data_file)
        else:
            move_file(pwd, data_files)
