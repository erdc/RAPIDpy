# -*- coding: utf-8 -*-
##
##  test_postprocess.py
##  RAPIDpy
##
##  Created by Alan D. Snow.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##

from datetime import datetime
from glob import glob
from netCDF4 import Dataset
from nose.tools import ok_
from numpy.testing import assert_almost_equal
import os
from shutil import copy
import unittest

#local import
from RAPIDpy.rapid import RAPID
from RAPIDpy.postprocess.generate_return_periods import generate_return_periods
from RAPIDpy.postprocess.generate_seasonal_averages import generate_seasonal_averages

from RAPIDpy.helper_functions import (compare_csv_decimal_files,
                                      remove_files)

class TestRAPIDInflow(unittest.TestCase):
    def setUp(self):
        #define global variables
        MAIN_TESTS_FOLDER = os.path.dirname(os.path.abspath(__file__))
        self.INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'data')
        self.COMPARE_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'compare')
        self.OUTPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'output')

    def test_generate_return_periods(self):
        """
        Checks generating return period data from RAPID Qout
        """
        print("TEST 1: TEST GENERATE RETURN PERIOD FROM RAPID QOUT")
        return_periods_file_name = 'return_periods_erai_t511_24hr_19800101to19861231.nc'
        generated_return_periods_file = os.path.join(self.OUTPUT_DATA_PATH, return_periods_file_name)
        generate_return_periods(qout_file=os.path.join(self.INPUT_DATA_PATH, 'Qout_erai_t511_24hr_19800101to19861231.nc'), 
                                return_period_file=generated_return_periods_file
                                )
                                
        compare_return_periods_file = os.path.join(self.COMPARE_DATA_PATH, return_periods_file_name)
        
        #check other info in netcdf file
        d1 = Dataset(generated_return_periods_file)
        d2 = Dataset(compare_return_periods_file)
        assert_almost_equal(d1.variables['max_flow'][:], d2.variables['max_flow'][:], decimal=5)
        assert_almost_equal(d1.variables['return_period_20'][:], d2.variables['return_period_20'][:], decimal=5)
        assert_almost_equal(d1.variables['return_period_10'][:], d2.variables['return_period_10'][:], decimal=5)
        assert_almost_equal(d1.variables['return_period_2'][:], d2.variables['return_period_2'][:], decimal=5)
        if 'rivid' in d2.variables.keys():
            ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
        if 'lat' in d2.variables.keys():
            ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
        if 'lon' in d2.variables.keys():
            ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
        d1.close()
        d2.close()
        
    def test_generate_seasonal_averages(self):
        """
        Checks generating seasonal average data from RAPID Qout
        """
        print("TEST 2: TEST GENERATE SEASONAL AVERAGE FROM RAPID QOUT")
        seasonal_averages_file_name = 'seasonal_averages_erai_t511_24hr_19800101to19861231.nc'
        generated_seasonal_averages_file = os.path.join(self.OUTPUT_DATA_PATH, seasonal_averages_file_name)
        generate_seasonal_averages(qout_file=os.path.join(self.INPUT_DATA_PATH, 'Qout_erai_t511_24hr_19800101to19861231.nc'), 
                                   seasonal_average_file=generated_seasonal_averages_file
                                   )
                                
        compare_seasonal_averages_file = os.path.join(self.COMPARE_DATA_PATH, seasonal_averages_file_name)
        
        #check other info in netcdf file
        d1 = Dataset(generated_seasonal_averages_file)
        d2 = Dataset(compare_seasonal_averages_file)
        assert_almost_equal(d1.variables['average_flow'][:], d2.variables['average_flow'][:], decimal=5)
        assert_almost_equal(d1.variables['std_dev_flow'][:], d2.variables['std_dev_flow'][:], decimal=5)
        if 'rivid' in d2.variables.keys():
            ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
        if 'lat' in d2.variables.keys():
            ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
        if 'lon' in d2.variables.keys():
            ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
        d1.close()
        d2.close()
        
        
    def test_generate_seasonal_qinit(self):
        """
        Checks generating seasonal qinit from rapid Qout
        """
        print("TEST 3: TEST GENERATE SEASONAL QINIT FROM RAPID QOUT")
    
        rapid_manager = RAPID(Qout_file=os.path.join(self.INPUT_DATA_PATH, 'Qout_erai_t511_24hr_19800101to19861231.nc'),
                              rapid_connect_file=os.path.join(self.COMPARE_DATA_PATH, 'gis', 'x-x', 'rapid_connect.csv')
                              )
        
        seasonal_init_file_name = 'Qinit_seasonal_avg_jan_1.csv'
        generated_seasonal_init_file = os.path.join(self.OUTPUT_DATA_PATH, seasonal_init_file_name)

        rapid_manager.generate_seasonal_intitialization(qinit_file=generated_seasonal_init_file,
                                                        datetime_start_initialization=datetime(1984,1,1))
        
        compare_seasonal_init_file = os.path.join(self.COMPARE_DATA_PATH, generated_seasonal_init_file)
        
        compare_csv_decimal_files(generated_seasonal_init_file,compare_seasonal_init_file)
        
        
    def tearDown(self):
        #remove unused data
        remove_files(*[f for f in glob(os.path.join(self.OUTPUT_DATA_PATH,"*")) if not f.endswith(".gitignore")])


if __name__ == '__main__':
    unittest.main()
