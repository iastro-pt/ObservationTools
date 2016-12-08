
# Radial Velocity calculations:

# Goals:
# To calcualte when the radial velocity is different by a certian value.
# Plot radial velocity phase curves. Indicating obtained measurement locations

import argparse


def _parser():
    parser = argparse.ArgumentParser(description='Radial velocity plotting')
    parser.add_argument('params', help='RV parameters filename')
    # parser.add_argument('-d', '--date', default='today',
    #                     help='Date in format YYYY-MM-DD. Default is today.')
    # #  prediciting separated lines
    # parser.add_argument('-r', '--rv_diff', help='RV difference threshold to find')
    # parser.add_argument('-o', '--obs_times',  help='Times of previous observations', nargs='+')
    # parser.add_argument('-f', '--files', help='Params and obs-times are file'
    #                    ' names to open', action='store_true')
    # parser.add_argument('-m', '--mode', help='Display mode '
    #                     ' e.g. phase or time plot. Default="phase"',
    #                     choices=['time', 'phase'], default='phase')
    return parser.parse_args()


def main(params):  # obs_times=None, mode='phase', rv_diff=None
    """ Do main stuff """

    # Load in params and store as a dictionary
    parameters = dict()
    with open(params, 'r') as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                par, val = line.lower().split('=')
                parameters[par.strip()] = val.strip()

    # Turn most parameters to floats.
    for key in parameters.keys():
        if key in ['mean_val', 'k1', 'omega', 'eccentricity', 'tau', 'Period',
                     'm_star', 'msini ', 'm_true']:
            parameters[key] = float(parameters[key])


    print(parameters)


if __name__ == '__main__':

    args = vars(_parser())
    # star_name = args.pop('star_name')
    opts = {k: args[k] for k in args}

    main(**opts)
