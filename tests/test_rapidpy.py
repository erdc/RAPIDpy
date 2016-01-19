# -*- coding: utf-8 -*-
##
##  test_rapidpy.py
##  RAPIDpy
##
##  Created by Alan D. Snow 2016.
##  Copyright © 2015 Alan D Snow. All rights reserved.
##

from filecmp import cmp as fcmp
from nose.tools import raises, ok_
import os
from shutil import copy


#local import
import sys
from RAPIDpy.rapid import RAPID

def test_generate_rapid_input_file():
    """
    Checks RAPID input file generation with valid input
    """
    main_tests_folder = os.path.dirname(os.path.abspath(__file__))
    
    compare_data_path = os.path.join(main_tests_folder, 'compare_data')
    output_data_path = os.path.join(main_tests_folder, 'tmp_output_files')

    print "TEST 2: GENERATE NAMELIST FILE"
    rapid_manager = RAPID(rapid_executable_location="",
                          use_all_processors=True,                          
                          ZS_TauR = 24*3600, #duration of routing procedure (time step of runoff data)
                          ZS_dtR = 15*60, #internal routing time step
                          ZS_TauM = 12*24*3600, #total simulation time 
                          ZS_dtM = 24*3600 #input time step 
                         )
    rapid_manager.update_parameters(rapid_connect_file='rapid_connect.csv',
                                    Vlat_file='m3_riv.nc',
                                    riv_bas_id_file='riv_bas_id.csv',
                                    k_file='k.csv',
                                    x_file='x.csv',
                                    Qout_file='Qout.nc'
                                    )
    generated_input_file = os.path.join(output_data_path, 
                                        "rapid_namelist-GENERATE.txt")
    rapid_manager.generate_namelist_file(generated_input_file)
    generated_input_file_solution = os.path.join(compare_data_path, 
                                                 "rapid_namelist-GENERATE-SOLUTION.txt")


    ok_(fcmp(generated_input_file, generated_input_file_solution))
    
    try:
        os.remove(generated_input_file)
    except OSError:
        pass

if __name__ == '__main__':
    import nose
    nose.main()

