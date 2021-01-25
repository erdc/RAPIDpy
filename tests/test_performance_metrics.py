# -*- coding: utf-8 -*-
#
#  test_performance_metrics.py
#
#  RAPIDpy
#
#  Created by Alan D. Snow.
#  Copyright © 2016 Alan D Snow. All rights reserved.
#

import numpy as np
import tempfile
import os

from RAPIDpy.postprocess import find_goodness_of_fit_csv_2

def test_find_goodness_of_fit_csv_2():
    testdir = os.path.dirname(os.path.abspath(__file__))
    
    tmpfile = tempfile.TemporaryFile()
    obsfile = os.path.join(testdir, 'data/obs_multiple_rivid.csv')
    simfile = os.path.join(testdir, 'data/sim_multiple_rivid.csv')
    compare = os.path.join(testdir, 'compare/rapid_performance_metrics.csv')

    find_goodness_of_fit_csv_2(obsfile, simfile, tmpfile)

    # Simulate closing and reopening tmpfile.
    tmpfile.seek(0)
  
    out = np.genfromtxt(tmpfile, delimiter=',', skip_header=1, 
                        usecols=(1,2), dtype=float)
    benchmark = np.genfromtxt(compare, delimiter=',', skip_header=1, 
                              usecols=(1,2), dtype=float)

    np.testing.assert_array_equal(out, benchmark) 

    tmpfile.close()

if __name__=='__main__':
    st = time.time()
    test_find_goodness_of_fit_csv_2()
    e = time.time() - st
    print('time:', e)
