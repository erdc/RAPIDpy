#! /usr/bin/env python

import sys

from RAPIDpy.postprocess import find_goodness_of_fit

if __name__=='__main__':
    obs_file = sys.argv[1]
    sim_file = sys.argv[2]
    out_file = sys.argv[3]

    find_goodness_of_fit_csv_2(obs_file, sim_file, out_file)
