# -*- coding: utf-8 -*-
##
##  test_inflow.py
##  RAPIDpy
##
##  Created by Alan D. Snow.
##  Copyright © 2016 Alan D Snow. All rights reserved.
##

from datetime import datetime
from glob import glob
import multiprocessing
from netCDF4 import Dataset
from nose.tools import ok_
from numpy.testing import assert_almost_equal
import os
from past.builtins import xrange
from shutil import copytree, rmtree
import unittest

#local import
from RAPIDpy.inflow import run_lsm_rapid_process
from RAPIDpy.inflow.CreateInflowFileFromERAInterimRunoff import CreateInflowFileFromERAInterimRunoff
from RAPIDpy.inflow.CreateInflowFileFromLDASRunoff import CreateInflowFileFromLDASRunoff
from RAPIDpy.inflow.CreateInflowFileFromWRFHydroRunoff import CreateInflowFileFromWRFHydroRunoff

from RAPIDpy.helper_functions import (compare_csv_decimal_files,
                                      remove_files)
#GLOBAL VARIABLES


class TestRAPIDInflow(unittest.TestCase):
    def setUp(self):
        #define global variables
        MAIN_TESTS_FOLDER = os.path.dirname(os.path.abspath(__file__))
        self.COMPARE_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'compare')
        self.INFLOW_COMPARE_DATA_PATH = os.path.join(self.COMPARE_DATA_PATH, 'inflow')
        self.LSM_INPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'data','lsm_grids')
        self.OUTPUT_DATA_PATH = os.path.join(MAIN_TESTS_FOLDER, 'output')
        self.RAPID_DATA_PATH = os.path.join(self.OUTPUT_DATA_PATH, 'input')
        self.RAPID_EXE_PATH = os.path.join(MAIN_TESTS_FOLDER,
                                           "..", "..",
                                           "rapid", "src", "rapid")
        self.CYGWIN_BIN_PATH = 'C:\\cygwin64\\bin'
        
    def _setup_automated(self, directory_name):
        """
        setup for automated method
        """
        rapid_input_path = os.path.join(self.RAPID_DATA_PATH, directory_name)
        rapid_output_path = os.path.join(self.OUTPUT_DATA_PATH, "output", directory_name)
        try:
            os.mkdir(self.RAPID_DATA_PATH)
        except OSError:
            pass
        
        try:
            os.mkdir(os.path.join(self.OUTPUT_DATA_PATH, "output"))
        except OSError:
            pass
        
        try:
            copytree(os.path.join(self.COMPARE_DATA_PATH, "gis", directory_name),
                     rapid_input_path)    
        except OSError:
            pass

        return(rapid_input_path, rapid_output_path) 
        
    def _setup_manual(self, directory_name):
        """
        setup for manual method
        """
        rapid_input_path, rapid_output_path = self._setup_automated(directory_name)

        try:
            os.mkdir(rapid_output_path)
        except OSError:
            pass

        return(rapid_input_path, rapid_output_path) 

    def _run_automatic(self, lsm_folder_name):
        """
        run for automatic method
        """
        #run main process    
        run_lsm_rapid_process(
            rapid_executable_location=self.RAPID_EXE_PATH,
            cygwin_bin_location=self.CYGWIN_BIN_PATH,
            rapid_io_files_location=self.OUTPUT_DATA_PATH,
            lsm_data_location=os.path.join(self.LSM_INPUT_DATA_PATH, lsm_folder_name), 
            simulation_start_datetime=datetime(1980, 1, 1),
            simulation_end_datetime=datetime(2014, 12, 31),
            generate_rapid_namelist_file=False,
            run_rapid_simulation=False,
            use_all_processors=True,
        )
        

    def test_run_era_interim_inflow(self):
        """
        Checks generating inflow file from ERA Interim LSM
        """
        print("TEST 1: TEST GENERATE INFLOW FILE FROM ERA INTERIM DATA")
        
        rapid_input_path, rapid_output_path = self._setup_automated("x-x")

        #run main process    
        run_lsm_rapid_process(
            rapid_executable_location=self.RAPID_EXE_PATH,
            cygwin_bin_location=self.CYGWIN_BIN_PATH,
            rapid_io_files_location=self.OUTPUT_DATA_PATH,
            lsm_data_location=os.path.join(self.LSM_INPUT_DATA_PATH, 'erai3'), 
            simulation_start_datetime=datetime(1980, 1, 1),
            simulation_end_datetime=datetime(2014, 1, 31),
            generate_rapid_namelist_file=False,
            run_rapid_simulation=True,
            generate_return_periods_file=False,
            generate_seasonal_initialization_file=False,
            generate_initialization_file=True,
            use_all_processors=True,
        )
        
        #CHECK OUTPUT    
        #m3_riv
        m3_file_name = "m3_riv_bas_erai_t511_3hr_20030121to20030122.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        self._compare_m3(generated_m3_file,generated_m3_file_solution)
        
        #qout file
        qout_file_name = "Qout_erai_t511_3hr_20030121to20030122.nc"
        generated_qout_file = os.path.join(rapid_output_path, qout_file_name)
        generated_qout_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, qout_file_name)
        d1 = Dataset(generated_qout_file)
        d2 = Dataset(generated_qout_file_solution)
        assert_almost_equal(d1.variables['Qout'][:], d2.variables['Qout'][:], decimal=5)
        ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
        ok_((d1.variables['time'][:] == d2.variables['time'][:]).all())
        if 'lat' in d2.variables.keys():
            ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
        if 'lon' in d2.variables.keys():
            ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
        d1.close()
        d2.close()
                                                     
        #initialization file
        qinit_file_name = "qinit_erai_t511_3hr_20030121to20030122.csv"
        generated_qinit_file = os.path.join(rapid_input_path, qinit_file_name)
        generated_qinit_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, qinit_file_name)
    
        ok_(compare_csv_decimal_files(generated_qinit_file, generated_qinit_file_solution))
        
        #additional cleanup
        remove_files(generated_qinit_file)
    
    def test_generate_erai_t511_inflow_manual(self):
        """
        Checks generating inflow file from ERA Interim t511 LSM manually
        """
        print("TEST 1.1: TEST GENERATE INFLOW FILE FROM ERA Interim t511 DATA MANUALLY")
        
        rapid_input_path, rapid_output_path = self._setup_manual("x-x")

        lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'erai3', '*.nc')))                 
        mp_lock = multiprocessing.Manager().Lock()
    
        inf_tool = CreateInflowFileFromERAInterimRunoff()
    
        m3_file_name = "m3_riv_bas_erai_t511_3hr_20030121to20030122.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    
        inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                          start_datetime_utc=datetime(2003,1,21),
                                          number_of_timesteps=len(lsm_file_list)*8,
                                          simulation_time_step_seconds=3*3600,
                                          in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                          in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                          land_surface_model_description="RAPID Inflow from ERA Interim (T511 Grid) 3 Hourly Runoff",
                                          modeling_institution="US Army Engineer Research and Development Center"
                                          )
                                          
        inf_tool.execute(nc_file_list=lsm_file_list, 
                         index_list=list(xrange(len(lsm_file_list))), 
                         in_weight_table=os.path.join(rapid_input_path, 'weight_era_t511.csv'), 
                         out_nc=generated_m3_file, 
                         grid_type='t511', 
                         mp_lock=mp_lock)
                          
        #CHECK OUTPUT
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        self._compare_m3(generated_m3_file,generated_m3_file_solution)

    def test_generate_nldas2_inflow(self):
        """
        Checks generating inflow file from NLDAS V2 LSM
        """
        print("TEST 2: TEST GENERATE INFLOW FILE FROM NLDAS V2 DATA")
        rapid_input_path, rapid_output_path = self._setup_automated("x-x")
        
        self._run_automatic('nldas2')
        
        #CHECK OUTPUT    
        #m3_riv
        m3_file_name = "m3_riv_bas_nasa_nldas_3hr_20030121to20030121.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        self._compare_m3(generated_m3_file,generated_m3_file_solution)
    
    def test_generate_nldas2_inflow2(self):
        """
        Checks generating inflow file from NLDAS V2 LSM manually
        """
        print("TEST 3.1: TEST GENERATE INFLOW FILE FROM ERA NLDAS V2 MANUALLY")
        rapid_input_path, rapid_output_path = self._setup_manual("x-x")
        
        lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'nldas2', '*.nc')))                 
        lsm_file_list = [lsm_file_list[nldas_index:nldas_index+3] for nldas_index in range(0, len(lsm_file_list), 3)\
                         if len(lsm_file_list[nldas_index:nldas_index+3])==3]
        mp_lock = multiprocessing.Manager().Lock()
    
        inf_tool = CreateInflowFileFromLDASRunoff(lat_dim="lat_110", 
                                                  lon_dim="lon_110", 
                                                  lat_var="lat_110", 
                                                  lon_var="lon_110", 
                                                  surface_runoff_var="SSRUNsfc_110_SFC_ave2h",
                                                  subsurface_runoff_var="BGRUNsfc_110_SFC_ave2h",
                                                  time_step_seconds=3*3600)
    
        m3_file_name = "m3_riv_bas_nasa_nldas_3hr_20030121to20030121.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    
        inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                          start_datetime_utc=datetime(2003,1,21),
                                          number_of_timesteps=len(lsm_file_list),
                                          simulation_time_step_seconds=3*3600,
                                          in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                          in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                          land_surface_model_description="RAPID Inflow from NLDAS Hourly Runoff",
                                          modeling_institution="US Army Engineer Research and Development Center"
                                          )
                                          
        inf_tool.execute(nc_file_list=lsm_file_list, 
                         index_list=list(xrange(len(lsm_file_list))), 
                         in_weight_table=os.path.join(rapid_input_path, 'weight_nldas.csv'), 
                         out_nc=generated_m3_file, 
                         grid_type='nldas', 
                         mp_lock=mp_lock)
                          
        #CHECK OUTPUT
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)

        self._compare_m3(generated_m3_file,generated_m3_file_solution)
        
    def test_generate_era20cm_inflow(self):
        """
        Checks generating inflow file from ERA 20CM LSM
        """
        print("TEST 3: TEST GENERATE INFLOW FILE FROM ERA 20CM DATA")
    
        rapid_input_path, rapid_output_path = self._setup_automated("x-x")
    
        run_lsm_rapid_process(
            rapid_executable_location=self.RAPID_EXE_PATH,
            cygwin_bin_location=self.CYGWIN_BIN_PATH,
            rapid_io_files_location=self.OUTPUT_DATA_PATH,
            lsm_data_location=os.path.join(self.LSM_INPUT_DATA_PATH, 'era20cm'), 
            simulation_start_datetime=datetime(1980, 1, 1),
            simulation_end_datetime=datetime(2014, 1, 31),
            ensemble_list=range(10),
            generate_rapid_namelist_file=False,
            run_rapid_simulation=False,
            use_all_processors=True,
        )
        
        for i in range(10):
            #CHECK OUTPUT    
            #m3_riv
            m3_file_name = "m3_riv_bas_era_20cm_t159_3hr_20000129to20000130_{0}.nc".format(i)
            generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
            generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)

            self._compare_m3(generated_m3_file,generated_m3_file_solution)
    
    def test_generate_era20cm_inflow2(self):
        """
        Checks generating inflow file from ERA 20CM LSM manually
        """
        print("TEST 3.1: TEST GENERATE INFLOW FILE FROM ERA 20CM DATA MANUALLY")

        rapid_input_path, rapid_output_path = self._setup_manual("x-x")

        lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'era20cm', '*_0.nc')))                 
        mp_lock = multiprocessing.Manager().Lock()
    
        inf_tool = CreateInflowFileFromERAInterimRunoff()
    
        m3_file_name = "m3_riv_bas_era_20cm_t159_3hr_20000129to20000130_0.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    
        inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                          start_datetime_utc=datetime(2014,8,20),
                                          number_of_timesteps=len(lsm_file_list)*8,
                                          simulation_time_step_seconds=3*3600,
                                          in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                          in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                          land_surface_model_description="RAPID Inflow from ERA 20CM (T159 Grid) 3 Hourly Runoff",
                                          modeling_institution="US Army Engineer Research and Development Center"
                                          )
                                          
        inf_tool.execute(nc_file_list=lsm_file_list, 
                         index_list=list(xrange(len(lsm_file_list))), 
                         in_weight_table=os.path.join(rapid_input_path, 'weight_era_t159.csv'), 
                         out_nc=generated_m3_file, 
                         grid_type='t159', 
                         mp_lock=mp_lock)
                          
        #CHECK OUTPUT
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        self._compare_m3(generated_m3_file,generated_m3_file_solution)
        
    def test_generate_erai_t255_inflow(self):
        """
        Checks generating inflow file from ERA Interim t255 LSM
        """
        print("TEST 4: TEST GENERATE INFLOW FILE FROM ERA Interim t255 DATA")
        rapid_input_path, rapid_output_path = self._setup_automated("x-x")

        #run main process
        self._run_automatic('erai3t255')
        
        #CHECK OUTPUT    
        #m3_riv
        m3_file_name = "m3_riv_bas_erai_t255_3hr_20140820to20140821.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        
        self._compare_m3(generated_m3_file,generated_m3_file_solution)

        
    def test_generate_erai_t255_inflow2(self):
        """
        Checks generating inflow file from ERA Interim t255 LSM manually
        """
        print("TEST 4.1: TEST GENERATE INFLOW FILE FROM ERA Interim t255 DATA MANUALLY")
        
        rapid_input_path, rapid_output_path = self._setup_manual("x-x")

        lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'erai3t255', '*.nc')))                 
        mp_lock = multiprocessing.Manager().Lock()
    
        inf_tool = CreateInflowFileFromERAInterimRunoff()
    
        m3_file_name = "m3_riv_bas_erai_t255_3hr_20140820to20140821.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    
        inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                          start_datetime_utc=datetime(2014,8,20),
                                          number_of_timesteps=len(lsm_file_list)*8,
                                          simulation_time_step_seconds=3*3600,
                                          in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                          in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                          land_surface_model_description="RAPID Inflow from ERA Interim (T255 Grid) 3 Hourly Runoff",
                                          modeling_institution="US Army Engineer Research and Development Center"
                                          )
                                          
        inf_tool.execute(nc_file_list=lsm_file_list, 
                         index_list=list(xrange(len(lsm_file_list))), 
                         in_weight_table=os.path.join(rapid_input_path, 'weight_era_t255.csv'), 
                         out_nc=generated_m3_file, 
                         grid_type='t255', 
                         mp_lock=mp_lock)
                          
        #CHECK OUTPUT
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        self._compare_m3(generated_m3_file,generated_m3_file_solution)
        
    def test_generate_gldas2_inflow(self):
        """
        Checks generating inflow file from GLDAS V2 LSM
        """
        print("TEST 5: TEST GENERATE INFLOW FILE FROM GLDAS V2 DATA")
        rapid_input_path, rapid_output_path = self._setup_automated("x-x")
        
        #run main process
        self._run_automatic('gldas2')

        #CHECK OUTPUT    
        #m3_riv
        m3_file_name = "m3_riv_bas_nasa_gldas2_3hr_20101231to20101231.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)

        self._compare_m3(generated_m3_file,generated_m3_file_solution)
        
    def test_generate_gldas2_inflow2(self):
        """
        Checks generating inflow file from GLDAS V2 LSM manually
        """
        print("TEST 5.1: TEST GENERATE INFLOW FILE FROM GLDAS V2 DATA MANUALLY")
        rapid_input_path, rapid_output_path = self._setup_manual("x-x")

        lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'gldas2', '*.nc4')))                 
        mp_lock = multiprocessing.Manager().Lock()
    
        inf_tool = CreateInflowFileFromLDASRunoff(lat_dim="lat", 
                                                  lon_dim="lon", 
                                                  lat_var="lat", 
                                                  lon_var="lon", 
                                                  surface_runoff_var="Qs_acc",
                                                  subsurface_runoff_var="Qsb_acc",
                                                  time_step_seconds=3*3600)
    
        m3_file_name = "m3_riv_bas_nasa_gldas2_3hr_20101231to20101231.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    
        inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                          start_datetime_utc=datetime(2010,12,31),
                                          number_of_timesteps=len(lsm_file_list),
                                          simulation_time_step_seconds=3*3600,
                                          in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                          in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                          land_surface_model_description="RAPID Inflow from GLDAS2.0 LIS land surface model 3 hourly runoff",
                                          modeling_institution="US Army Engineer Research and Development Center"
                                          )
                                          
        inf_tool.execute(nc_file_list=lsm_file_list, 
                         index_list=list(xrange(len(lsm_file_list))), 
                         in_weight_table=os.path.join(rapid_input_path, 'weight_gldas2.csv'), 
                         out_nc=generated_m3_file, 
                         grid_type='gldas2', 
                         mp_lock=mp_lock)
                          
        #CHECK OUTPUT
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)

        self._compare_m3(generated_m3_file,generated_m3_file_solution)
        
    def test_generate_lis_inflow(self):
        """
        Checks generating inflow file from LIS LSM
        """
        print("TEST 6: TEST GENERATE INFLOW FILE FROM LIS DATA")
        rapid_input_path, rapid_output_path = self._setup_automated("u-k")
        
        #run main process
        self._run_automatic('lis')

        #CHECK OUTPUT    
        #m3_riv
        m3_file_name = "m3_riv_bas_nasa_lis_3hr_20110121to20110121.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)

        self._compare_m3(generated_m3_file,generated_m3_file_solution)

    def test_generate_lis_inflow2(self):
        """
        Checks generating inflow file from LIS LSM manually
        """
        print("TEST 6.1: TEST GENERATE INFLOW FILE FROM LIS DATA MANUALLY")
        rapid_input_path, rapid_output_path = self._setup_manual("u-k")
        
        lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'lis', '*.nc')))                 
        lsm_file_list = [lsm_file_list[nldas_index:nldas_index+3] for nldas_index in range(0, len(lsm_file_list), 3)\
                         if len(lsm_file_list[nldas_index:nldas_index+3])==3]
        mp_lock = multiprocessing.Manager().Lock()
    
        inf_tool = CreateInflowFileFromLDASRunoff(lat_dim="north_south", 
                                                  lon_dim="east_west", 
                                                  lat_var="lat", 
                                                  lon_var="lon", 
                                                  surface_runoff_var="Qs_inst",
                                                  subsurface_runoff_var="Qsb_inst",
                                                  time_step_seconds=3*3600)
    
        m3_file_name = "m3_riv_bas_nasa_lis_3hr_20110121to20110121.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    
        inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                          start_datetime_utc=datetime(2008,3,3),
                                          number_of_timesteps=len(lsm_file_list),
                                          simulation_time_step_seconds=3*3600,
                                          in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                          in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                          land_surface_model_description="RAPID Inflow from NASA GSFC LIS hourly runoff",
                                          modeling_institution="US Army Engineer Research and Development Center"
                                          )
                                          
        inf_tool.execute(nc_file_list=lsm_file_list, 
                         index_list=list(xrange(len(lsm_file_list))), 
                         in_weight_table=os.path.join(rapid_input_path, 'weight_lis.csv'), 
                         out_nc=generated_m3_file, 
                         grid_type='lis', 
                         mp_lock=mp_lock)
                          
        #CHECK OUTPUT
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)

        self._compare_m3(generated_m3_file,generated_m3_file_solution)
    
    def test_generate_joules_inflow(self):
        """
        Checks generating inflow file from Joules LSM
        """
        print("TEST 7: TEST GENERATE INFLOW FILE FROM Joules DATA")
        rapid_input_path, rapid_output_path = self._setup_automated("u-k")
        
        #run main process    
        run_lsm_rapid_process(
            rapid_executable_location=self.RAPID_EXE_PATH,
            cygwin_bin_location=self.CYGWIN_BIN_PATH,
            rapid_io_files_location=self.OUTPUT_DATA_PATH,
            lsm_data_location=os.path.join(self.LSM_INPUT_DATA_PATH, 'joules'), 
            simulation_start_datetime=datetime(1980, 1, 1),
            simulation_end_datetime=datetime(2014, 1, 31),
            file_datetime_re_pattern = r'\d{8}_\d{2}',
            file_datetime_pattern = "%Y%m%d_%H",      
            generate_rapid_namelist_file=False,
            run_rapid_simulation=False,
            use_all_processors=True,
        )
        
        #CHECK OUTPUT
        m3_file_name = "m3_riv_bas_met_office_joules_3hr_20080803to20080803.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        
        #check other info in netcdf file
        self._compare_m3(generated_m3_file, generated_m3_file_solution)

    
    def test_generate_joules_inflow2(self):
        """
        Checks generating inflow file from Joules LSM manually
        """
        print("TEST 7.1: TEST GENERATE INFLOW FILE FROM Joules DATA MANUALLY")
        rapid_input_path, rapid_output_path = self._setup_manual("u-k")
        
        lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'joules', '*.nc')))                 
        lsm_file_list = [lsm_file_list[nldas_index:nldas_index+3] for nldas_index in range(0, len(lsm_file_list), 3)\
                         if len(lsm_file_list[nldas_index:nldas_index+3])==3]
        mp_lock = multiprocessing.Manager().Lock()
    
        inf_tool = CreateInflowFileFromLDASRunoff(lat_dim="north_south", 
                                                  lon_dim="east_west", 
                                                  lat_var="north_south", 
                                                  lon_var="east_west", 
                                                  surface_runoff_var="Qs_inst",
                                                  subsurface_runoff_var="Qsb_inst",
                                                  time_step_seconds=3*3600)
    
        m3_file_name = "m3_riv_bas_met_office_joules_3hr_20080803to20080803.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    
        inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                          start_datetime_utc=datetime(2008,3,3),
                                          number_of_timesteps=len(lsm_file_list),
                                          simulation_time_step_seconds=3*3600,
                                          in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                          in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                          land_surface_model_description="RAPID Inflow from Met Office Joules Hourly Runoff",
                                          modeling_institution="US Army Engineer Research and Development Center"
                                          )
                                          
        inf_tool.execute(nc_file_list=lsm_file_list, 
                         index_list=list(xrange(len(lsm_file_list))), 
                         in_weight_table=os.path.join(rapid_input_path, 'weight_joules.csv'), 
                         out_nc=generated_m3_file, 
                         grid_type='joules', 
                         mp_lock=mp_lock)
                          
        #CHECK OUTPUT
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        
        self._compare_m3(generated_m3_file, generated_m3_file_solution)


    def test_generate_erai_t511_24_inflow(self):
        """
        Checks generating inflow file from ERA Interim t511 24hr LSM
        """
        print("TEST 8: TEST GENERATE INFLOW FILE FROM ERA Interim t511 24hr DATA")
        rapid_input_path, rapid_output_path = self._setup_automated("x-x")

        #run main process
        self._run_automatic('erai24')
        
        #CHECK OUTPUT    
        #m3_riv
        m3_file_name = "m3_riv_bas_erai_t511_24hr_19990109to19990110.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        
        self._compare_m3(generated_m3_file,generated_m3_file_solution)

        
    def test_generate_erai_t511_24_inflow2(self):
        """
        Checks generating inflow file from ERA Interim t511 24hr LSM manually
        """
        print("TEST 8.1: TEST GENERATE INFLOW FILE FROM ERA Interim t511 24hr DATA MANUALLY")
        
        rapid_input_path, rapid_output_path = self._setup_manual("x-x")

        lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'erai24', '*.nc')))                 
        mp_lock = multiprocessing.Manager().Lock()
    
        inf_tool = CreateInflowFileFromERAInterimRunoff()
    
        m3_file_name = "m3_riv_bas_erai_t511_24hr_19990109to19990110.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
    
        inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                          start_datetime_utc=datetime(1999,1,9),
                                          number_of_timesteps=len(lsm_file_list),
                                          simulation_time_step_seconds=24*3600,
                                          in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                          in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                          land_surface_model_description="RAPID Inflow from ERA Interim (T511 Grid) 24 Hourly Runoff",
                                          modeling_institution="US Army Engineer Research and Development Center"
                                          )
                                          
        inf_tool.execute(nc_file_list=lsm_file_list, 
                         index_list=list(xrange(len(lsm_file_list))), 
                         in_weight_table=os.path.join(rapid_input_path, 'weight_era_t511.csv'), 
                         out_nc=generated_m3_file, 
                         grid_type='t511', 
                         mp_lock=mp_lock)
                          
        #CHECK OUTPUT
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        self._compare_m3(generated_m3_file,generated_m3_file_solution)

    def test_generate_wrf_inflow(self):
        """
        Checks generating inflow file from WRF LSM
        """
        print("TEST 9: TEST GENERATE INFLOW FILE FROM WRF DATA")
        rapid_input_path, rapid_output_path = self._setup_automated("m-s")

        #run main process
        self._run_automatic('wrf')
        
        #CHECK OUTPUT    
        #m3_riv
        m3_file_name = "m3_riv_bas_wrf_wrf_1hr_20080601to20080601.nc"
        generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
        generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
        
        self._compare_m3(generated_m3_file,generated_m3_file_solution)

    def test_generate_wrf_inflow2(self):
         """
         Checks generating inflow file from WRF LSM manually
         """
         print("TEST 9.1: TEST GENERATE INFLOW FILE FROM WRF DATA MANUALLY")
         
         rapid_input_path, rapid_output_path = self._setup_manual("m-s")
 
         lsm_file_list =  sorted(glob(os.path.join(self.LSM_INPUT_DATA_PATH, 'wrf', '*.nc')))                 
         mp_lock = multiprocessing.Manager().Lock()
     
         inf_tool = CreateInflowFileFromWRFHydroRunoff(lat_dim="south_north", 
                                                       lon_dim="west_east", 
                                                       lat_var="XLAT", 
                                                       lon_var="XLONG", 
                                                       surface_runoff_var="SFROFF",
                                                       subsurface_runoff_var="UDROFF",
                                                       time_step_seconds=3*3600)
     
         m3_file_name = "m3_riv_bas_wrf_wrf_1hr_20080601to20080601.nc"
         generated_m3_file = os.path.join(rapid_output_path, m3_file_name)
     
         inf_tool.generateOutputInflowFile(out_nc=generated_m3_file,
                                           start_datetime_utc=datetime(2008,6,1),
                                           number_of_timesteps=len(lsm_file_list),
                                           simulation_time_step_seconds=3*3600,
                                           in_rapid_connect_file=os.path.join(rapid_input_path, 'rapid_connect.csv'),
                                           in_rivid_lat_lon_z_file=os.path.join(rapid_input_path, 'comid_lat_lon_z.csv'),
                                           land_surface_model_description="RAPID Inflow from WRF Hourly Runoff",
                                           modeling_institution="US Army Engineer Research and Development Center"
                                           )
                                           
         inf_tool.execute(nc_file_list=lsm_file_list, 
                          index_list=list(xrange(len(lsm_file_list))), 
                          in_weight_table=os.path.join(rapid_input_path, 'weight_wrf.csv'), 
                          out_nc=generated_m3_file, 
                          grid_type='wrf', 
                          mp_lock=mp_lock)
                           
         #CHECK OUTPUT
         generated_m3_file_solution = os.path.join(self.INFLOW_COMPARE_DATA_PATH, m3_file_name)
         self._compare_m3(generated_m3_file,generated_m3_file_solution)
         
    def _compare_m3(self, generated_m3_file, generated_m3_file_solution):
        
        #check other info in netcdf file
        d1 = Dataset(generated_m3_file)
        d2 = Dataset(generated_m3_file_solution)
        assert_almost_equal(d1.variables['m3_riv'][:], d2.variables['m3_riv'][:], decimal=5)
        if 'rivid' in d2.variables.keys():
            ok_((d1.variables['rivid'][:] == d2.variables['rivid'][:]).all())
        if 'lat' in d2.variables.keys():
            ok_((d1.variables['lat'][:] == d2.variables['lat'][:]).all())
        if 'lon' in d2.variables.keys():
            ok_((d1.variables['lon'][:] == d2.variables['lon'][:]).all())
        d1.close()
        d2.close()
        
    def tearDown(self):
        rmtree(os.path.join(self.OUTPUT_DATA_PATH, "input"))
        rmtree(os.path.join(self.OUTPUT_DATA_PATH, "output"))

if __name__ == '__main__':
    unittest.main()
