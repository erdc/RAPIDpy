# -*- coding: utf-8 -*-
##
##  lsm_rapid_process.py
##  RAPIDpy
##
##  Created by Alan D. Snow.
##  Copyright © 2015-2016 Alan D Snow. All rights reserved.
##  License: BSD 3-Clause

from datetime import datetime
import multiprocessing
from netCDF4 import Dataset
import os
import re
import traceback
#local imports
from ..rapid import RAPID
from .CreateInflowFileFromERAInterimRunoff import CreateInflowFileFromERAInterimRunoff
from .CreateInflowFileFromLDASRunoff import CreateInflowFileFromLDASRunoff
from .CreateInflowFileFromWRFHydroRunoff import CreateInflowFileFromWRFHydroRunoff
from .generate_return_periods import generate_return_periods
from .utilities import (case_insensitive_file_search,
                        get_valid_watershed_list,
                        get_watershed_subbasin_from_folder,
                        partition)


#------------------------------------------------------------------------------
#MULTIPROCESSING FUNCTION
#------------------------------------------------------------------------------
def generate_inflows_from_runoff(args):
    """
    prepare runoff inflow file for rapid
    """
    watershed = args[0]
    subbasin = args[1]
    runoff_file_list = args[2]
    file_index_list = args[3]
    weight_table_file = args[4]
    grid_type = args[5]
    rapid_inflow_file = args[6]
    RAPID_Inflow_Tool = args[7]
    mp_lock = args[8]

    time_start_all = datetime.utcnow()

    if not isinstance(runoff_file_list, list): 
        runoff_file_list = [runoff_file_list]
    else:
        runoff_file_list = runoff_file_list

    if not isinstance(file_index_list, list): 
        file_index_list = [file_index_list]
    else:
        file_index_list = file_index_list
    if runoff_file_list and file_index_list:
        #prepare ECMWF file for RAPID
        print("Runoff downscaling for: {0} {1}".format(watershed, subbasin))
        index_string = "Index: {0}".format(file_index_list[0])
        if len(file_index_list) > 1:
            index_string += " to {0}".format(file_index_list[-1])
        print(index_string)
        runoff_string = "File(s): {0}".format(runoff_file_list[0])
        if len(runoff_file_list) > 1:
            runoff_string += " to {0}".format(runoff_file_list[-1])
        print(runoff_string)
          
           
        print("Converting inflow")
        try:
            RAPID_Inflow_Tool.execute(nc_file_list=runoff_file_list,
                                      index_list=file_index_list,
                                      in_weight_table=weight_table_file,
                                      out_nc=rapid_inflow_file,
                                      grid_type=grid_type,
                                      mp_lock=mp_lock,
                                      )
        except Exception:
            # This prints the type, value, and stack trace of the
            # current exception being handled.
            traceback.print_exc()
            raise
            
        time_finish_ecmwf = datetime.utcnow()
        print("Time to convert inflows: {0}".format(time_finish_ecmwf-time_start_all))

#------------------------------------------------------------------------------
#MAIN PROCESS
#------------------------------------------------------------------------------
def run_lsm_rapid_process(rapid_executable_location,
                          rapid_io_files_location,
                          lsm_data_location,
                          simulation_start_datetime,
                          simulation_end_datetime=datetime.utcnow(),
                          ensemble_list=[None],
                          generate_rapid_namelist_file=True,
                          run_rapid_simulation=True,
                          generate_return_periods_file=False,
                          generate_seasonal_initialization_file=False,
                          generate_initialization_file=False,
                          use_all_processors=True,
                          num_processors=1,
                          cygwin_bin_location="",
                          modeling_institution="US Army Engineer Research and Development Center"
                          ):
    """
    This is the main process to generate inflow for RAPID and to run RAPID
    """
    time_begin_all = datetime.utcnow()

    #use all processors makes precedent over num_processors arg
    if use_all_processors == True:
        NUM_CPUS = multiprocessing.cpu_count()
    elif num_processors > multiprocessing.cpu_count():
        print("WARNING: Num processors requested exceeded max. Set to max ...")
        NUM_CPUS = multiprocessing.cpu_count()
    else:
        NUM_CPUS = num_processors

    #get list of correclty formatted rapid input directories in rapid directory
    rapid_input_directories = get_valid_watershed_list(os.path.join(rapid_io_files_location, 'input'))

    for ensemble in ensemble_list:
        ensemble_file_ending = ".nc"
        ensemble_file_ending4 = ".nc4"
        if ensemble != None:
            ensemble_file_ending = "_{0}.nc".format(ensemble)
            ensemble_file_ending4 = "_{0}.nc4".format(ensemble)
        #get list of files
        lsm_file_list = []
        for subdir, dirs, files in os.walk(lsm_data_location, followlinks=True):
            for lsm_file in files:
                if lsm_file.endswith(ensemble_file_ending) or lsm_file.endswith(ensemble_file_ending4):
                    lsm_file_list.append(os.path.join(subdir, lsm_file))
        
        lsm_file_list_subset = []
        
        for lsm_file in sorted(lsm_file_list):
            match = re.search(r'\d{8}', lsm_file)
            file_date = datetime.strptime(match.group(0), "%Y%m%d")
            if file_date > simulation_end_datetime:
                break
            if file_date >= simulation_start_datetime:
                lsm_file_list_subset.append(os.path.join(subdir, lsm_file))
        print(lsm_file_list_subset[0])
        actual_simulation_start_datetime = datetime.strptime(re.search(r'\d{8}', lsm_file_list_subset[0]).group(0), "%Y%m%d")
        print(lsm_file_list_subset[-1])
        actual_simulation_end_datetime = datetime.strptime(re.search(r'\d{8}', lsm_file_list_subset[-1]).group(0), "%Y%m%d")
        
        lsm_file_list = sorted(lsm_file_list_subset)
        
        #check to see what kind of file we are dealing with
        lsm_example_file = Dataset(lsm_file_list[0])

        
        #INDENTIFY LAT/LON DIMENSIONS
        dim_list = lsm_example_file.dimensions.keys()

        latitude_dim = "lat"
        if 'latitude' in dim_list:
            latitude_dim = 'latitude'
        elif 'g0_lat_0' in dim_list:
            #GLDAS/NLDAS MOSAIC
            latitude_dim = 'g0_lat_0'
        elif 'lat_110' in dim_list:
            #NLDAS NOAH/VIC
            latitude_dim = 'lat_110'
        elif 'north_south' in dim_list:
            #LIS/Joules
            latitude_dim = 'north_south'
        elif 'south_north' in dim_list:
            #WRF Hydro
            latitude_dim = 'south_north'
        
        longitude_dim = "lon"
        if 'longitude' in dim_list:
            longitude_dim = 'longitude'
        elif 'g0_lon_1' in dim_list:
            #GLDAS/NLDAS MOSAIC
            longitude_dim = 'g0_lon_1'
        elif 'lon_110' in dim_list:
            #NLDAS NOAH/VIC
            longitude_dim = 'lon_110'
        elif 'east_west' in dim_list:
            #LIS/Joules
            longitude_dim = 'east_west'
        elif 'west_east' in dim_list:
            #WRF Hydro
            longitude_dim = 'west_east'

        lat_dim_size = len(lsm_example_file.dimensions[latitude_dim])
        lon_dim_size = len(lsm_example_file.dimensions[longitude_dim])

        #IDENTIFY VARIABLES
        var_list = lsm_example_file.variables.keys()

        latitude_var="lat"
        if 'latitude' in var_list:
            latitude_var = 'latitude'
        elif 'g0_lat_0' in var_list:
            latitude_var = 'g0_lat_0'
        elif 'lat_110' in var_list:
            latitude_var = 'lat_110'
        elif 'north_south' in var_list:
            latitude_var = 'north_south'
        elif 'XLAT' in var_list:
            latitude_var = 'XLAT'
    
        longitude_var="lon"
        if 'longitude' in var_list:
            longitude_var = 'longitude'
        elif 'g0_lon_1' in var_list:
            longitude_var = 'g0_lon_1'
        elif 'lon_110' in var_list:
            longitude_var = 'lon_110'
        elif 'east_west' in var_list:
            longitude_var = 'east_west'
        elif 'XLONG' in var_list:
            longitude_var = 'XLONG'

        surface_runoff_var=""
        subsurface_runoff_var=""
        snowmelt_runoff_var=""
        for var in var_list:
            if var.startswith("SSRUN"):
                #NLDAS/GLDAS
                surface_runoff_var = var
            elif var.startswith("BGRUN"):
                #NLDAS/GLDAS
                subsurface_runoff_var = var
            elif var == "Qs_acc":
                #GLDAS v2
                surface_runoff_var = var
            elif var == "Qsb_acc":
                #GLDAS v2
                subsurface_runoff_var = var
# SNOWMELT IS INCLUDED IN RUNOFF,
# HOWEVER YOU CAN ADD RUNOFF OUTPUT
# ANYWAY BY UNCOMMENTING:
#            elif var == "Qsm_acc":
#                #GLDAS v2
#                snowmelt_runoff_var = var
            elif var == "Qs_inst":
                #LIS
                surface_runoff_var = var
            elif var == "Qsb_inst":
                #LIS
                subsurface_runoff_var = var
            elif var == "SFROFF":
                #WRF Hydro
                surface_runoff_var = var
            elif var == "UDROFF":
                #WRF Hydro
                subsurface_runoff_var = var
            elif var.lower() == "ro":
                #ERA Interim
                surface_runoff_var = var
            


        #IDENTIFY GRID TYPE & TIME STEP

        try:
            time_dim = "time"
            if "Time" in lsm_example_file.dimensions:
                time_dim = "Time"
            file_size_time = len(lsm_example_file.dimensions[time_dim])
        except Exception as ex:
            print("ERROR: {0}".format(ex))
            print("Assuming time dimension is 1")
            file_size_time = 1

        out_file_ending = "{0}to{1}{2}".format(actual_simulation_start_datetime.strftime("%Y%m%d"), 
                                               actual_simulation_end_datetime.strftime("%Y%m%d"), 
                                               ensemble_file_ending)
            
        weight_file_name = ''
        grid_type = ''
        model_name = ''
        time_step = 0
        description = ""
        RAPID_Inflow_Tool = None
        total_num_time_steps = 0
        institution = ""
        title = ""
        try:
            institution = lsm_example_file.getncattr("institution")
        except AttributeError:
            pass
        try:
            title = lsm_example_file.getncattr("title")
        except AttributeError:
            pass
        
        if institution == "European Centre for Medium-Range Weather Forecasts" \
            or surface_runoff_var.lower() == "ro":
            #these are the ECMWF models
            if lat_dim_size == 361 and lon_dim_size == 720:
                print("Runoff file identified as ERA Interim Low Res (T255) GRID")
                #A) ERA Interim Low Res (T255)
                #Downloaded as 0.5 degree grid
                # dimensions:
                #	 longitude = 720 ;
                #	 latitude = 361 ;
                description = "ERA Interim (T255 Grid)"
                model_name = "erai"
                weight_file_name = r'weight_era_t255\.csv'
                grid_type = 't255'

            elif lat_dim_size == 512 and lon_dim_size == 1024:
                print("Runoff file identified as ERA Interim High Res (T511) GRID")
                #B) ERA Interim High Res (T511)
                # dimensions:
                #  lon = 1024 ;
                #  lat = 512 ;
                description = "ERA Interim (T511 Grid)"
                weight_file_name = r'weight_era_t511\.csv'
                model_name = "erai"
                grid_type = 't511'
            elif lat_dim_size == 161 and lon_dim_size == 320:
                print("Runoff file identified as ERA 20CM (T159) GRID")
                #C) ERA 20CM (T159) - 3hr - 10 ensembles
                #Downloaded as 1.125 degree grid
                # dimensions:
                #  longitude = 320 ;
                #  latitude = 161 ;    
                description = "ERA 20CM (T159 Grid)"
                weight_file_name = r'weight_era_t159\.csv'
                model_name = "era_20cm"
                grid_type = 't159'
            else:
                lsm_example_file.close()
                raise Exception("Unsupported grid size.")

            #time units are in hours
            if file_size_time == 1:
                time_step = 24*3600 #daily
                description += " Daily Runoff"
            elif file_size_time == 8:
                time_step = 3*3600 #3 hourly
                description += " 3 Hourly Runoff"
            else:
                lsm_example_file.close()
                raise Exception("Unsupported ECMWF time step.")

            total_num_time_steps=file_size_time*len(lsm_file_list)
            RAPID_Inflow_Tool = CreateInflowFileFromERAInterimRunoff()
                 
        elif institution == "NASA GSFC":
            if title == "GLDAS2.0 LIS land surface model output":
                print("Runoff file identified as GLDAS v2 LIS GRID")
                #TODO: SNOWMELT
                #this is the LIS model
                weight_file_name = r'weight_gldas2\.csv'
                grid_type = 'gldas2'
                description = "GLDAS2.0 LIS land surface model 3 hourly runoff"
                model_name = "nasa"
                
                #time units are in minutes
                if file_size_time == 1:
                    time_step = 3*3600 #3-hourly
                else:
                    lsm_example_file.close()
                    raise Exception("Unsupported GLDAS 2.0 time step.")
                
            else:
                print("Runoff file identified as LIS GRID")
                #this is the LIS model
                weight_file_name = r'weight_lis\.csv'
                grid_type = 'lis'
                description = "NASA GFC LIS hourly runoff"
                model_name = "nasa"
                #time units are in minutes
                if file_size_time == 1:
                    #time_step = 1*3600 #hourly
                    time_step = 3*3600 #3-hourly
                else:
                    lsm_example_file.close()
                    raise Exception("Unsupported LIS time step.")
            
            total_num_time_steps=file_size_time*len(lsm_file_list)
            RAPID_Inflow_Tool = CreateInflowFileFromLDASRunoff(latitude_dim,
                                                               longitude_dim,
                                                               latitude_var,
                                                               longitude_var,
                                                               surface_runoff_var,
                                                               subsurface_runoff_var,
                                                               time_step,
                                                               snowmelt_runoff_var)

        elif institution == "Met Office, UK":
            print("Runoff file identified as Joules GRID")
            #this is the LIS model
            weight_file_name = r'weight_joules\.csv'
            grid_type = 'joules'
            description = "Met Office Joules Hourly Runoff"
            model_name = "met_office"
            #time units are in minutes
            if file_size_time == 1:
                #time_step = 1*3600 #hourly
                time_step = 3*3600 #3-hourly
            else:
                lsm_example_file.close()
                raise Exception("Unsupported LIS time step.")
        
            total_num_time_steps=file_size_time*len(lsm_file_list)

            RAPID_Inflow_Tool = CreateInflowFileFromLDASRunoff(latitude_dim,
                                                               longitude_dim,
                                                               latitude_var,
                                                               longitude_var,
                                                               surface_runoff_var,
                                                               subsurface_runoff_var,
                                                               time_step)
        elif surface_runoff_var.startswith("SSRUN") \
            and subsurface_runoff_var.startswith("BGRUN"):

            model_name = "nasa"
            if lat_dim_size == 600 and lon_dim_size == 1440:
                print("Runoff file identified as GLDAS GRID")
                #GLDAS NC FILE
                #dimensions:
                #    g0_lat_0 = 600 ;
                #    g0_lon_1 = 1440 ;
                #variables
                #SSRUN_GDS0_SFC_ave1h (surface), BGRUN_GDS0_SFC_ave1h (subsurface)
                # or
                #SSRUNsfc_GDS0_SFC_ave1h (surface), BGRUNsfc_GDS0_SFC_ave1h (subsurface)
                description = "GLDAS 3 Hourly Runoff"
                weight_file_name = r'weight_gldas\.csv'
                grid_type = 'gldas'
 
                if file_size_time == 1:
                    time_step = 3*3600 #3 hourly
                else:
                    lsm_example_file.close()
                    raise Exception("Unsupported GLDAS time step.")
                
                total_num_time_steps=file_size_time*len(lsm_file_list)

            elif lat_dim_size <= 224 and lon_dim_size <= 464:
                print("Runoff file identified as NLDAS GRID")
                #NLDAS MOSAIC FILE
                #dimensions:
                #    g0_lat_0 = 224 ;
                #    g0_lon_1 = 464 ;
                #NLDAS NOAH/VIC FILE
                #dimensions:
                #    lat_110 = 224 ;
                #    lon_110 = 464 ;

                description = "NLDAS Hourly Runoff"
                weight_file_name = r'weight_nldas\.csv'
                grid_type = 'nldas'

                if file_size_time == 1:
                    #time_step = 1*3600 #hourly
                    time_step = 3*3600 #3 hourly
                else:
                    lsm_example_file.close()
                    raise Exception("Unsupported NLDAS time step.")
                
            else:
                lsm_example_file.close()
                raise Exception("Unsupported runoff grid.")

            RAPID_Inflow_Tool = CreateInflowFileFromLDASRunoff(latitude_dim,
                                                               longitude_dim,
                                                               latitude_var,
                                                               longitude_var,
                                                               surface_runoff_var,
                                                               subsurface_runoff_var,
                                                               time_step)
        else:
            title = ""
            try:
                title = lsm_example_file.getncattr("TITLE")
            except AttributeError:
                pass

            if "WRF" in title:
                description = "WRF-Hydro Hourly Runoff"
                weight_file_name = r'weight_wrf\.csv'
                grid_type = 'wrf_hydro'
                time_step = 1*3600 #1 hourly
                total_num_time_steps=file_size_time*len(lsm_file_list)

                RAPID_Inflow_Tool = CreateInflowFileFromWRFHydroRunoff(latitude_dim,
                                                                       longitude_dim,
                                                                       latitude_var,
                                                                       longitude_var,
                                                                       surface_runoff_var,
                                                                       subsurface_runoff_var,
                                                                       time_step)
            else:
                lsm_example_file.close()
                raise Exception("Unsupported runoff grid.")
        
        lsm_example_file.close()

    	#VALIDATING INPUT IF DIVIDING BY 3
        if grid_type == 'nldas' or grid_type == 'lis' or grid_type == 'joules':
            num_extra_files = file_size_time*len(lsm_file_list) % 3
            if num_extra_files != 0:
                print("WARNING: Number of files needs to be divisible by 3. Remainder is {0}".format(num_extra_files))
                print("This means your simulation will be truncated")
            total_num_time_steps=int(file_size_time*len(lsm_file_list)/3)

        out_file_ending = "{0}_{1}_{2}hr_{3}".format(model_name, grid_type, int(time_step/3600), out_file_ending)
        #set up RAPID manager
        rapid_manager = RAPID(rapid_executable_location=rapid_executable_location,
                              cygwin_bin_location=cygwin_bin_location,
                              num_processors=NUM_CPUS,
                              ZS_TauR=time_step, #duration of routing procedure (time step of runoff data)
                              ZS_dtR=15*60, #internal routing time step
                              ZS_TauM=total_num_time_steps*time_step, #total simulation time
                              ZS_dtM=time_step #RAPID recommended internal time step (1 day)
                             )
    
        #run ERA Interim processes
        for rapid_input_directory in rapid_input_directories:
            watershed, subbasin = get_watershed_subbasin_from_folder(rapid_input_directory)

            master_watershed_input_directory = os.path.join(rapid_io_files_location,
                                                            'input',
                                                            rapid_input_directory)
            master_watershed_output_directory = os.path.join(rapid_io_files_location,
                                                             'output',
                                                             rapid_input_directory)
            try:
                os.makedirs(master_watershed_output_directory)
            except OSError:
                pass
    
            #create inflow to dump data into
            master_rapid_runoff_file = os.path.join(master_watershed_output_directory, 
                                                    'm3_riv_bas_{0}'.format(out_file_ending))
            
            weight_table_file = case_insensitive_file_search(master_watershed_input_directory,
                                                             weight_file_name)

            try:
                in_rivid_lat_lon_z_file = case_insensitive_file_search(master_watershed_input_directory,
                                                                    r'comid_lat_lon_z\.csv')
            except Exception:
                in_rivid_lat_lon_z_file = ""
                print("WARNING: comid_lat_lon_z file not found. These will not be added ...")
                pass
            
            RAPID_Inflow_Tool.generateOutputInflowFile(out_nc=master_rapid_runoff_file,
                                                       start_datetime_utc=actual_simulation_start_datetime,
                                                       number_of_timesteps=total_num_time_steps,
                                                       simulation_time_step_seconds=time_step,
                                                       in_rapid_connect_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                                                          r'rapid_connect\.csv'),
                                                       in_rivid_lat_lon_z_file=in_rivid_lat_lon_z_file,
                                                       land_surface_model_description=description,
                                                       modeling_institution=modeling_institution
                                                       )
            job_combinations = []
            if grid_type == 'nldas' or grid_type == 'lis' or grid_type == 'joules':
                print("Grouping {0} in threes".format(grid_type))
                lsm_file_list = [lsm_file_list[nldas_index:nldas_index+3] for nldas_index in range(0, len(lsm_file_list), 3)\
                                 if len(lsm_file_list[nldas_index:nldas_index+3])==3]
            if len(lsm_file_list) < NUM_CPUS:
                NUM_CPUS = len(lsm_file_list)
            mp_lock = multiprocessing.Manager().Lock()
            partition_list, partition_index_list = partition(lsm_file_list, NUM_CPUS)
            for loop_index, cpu_grouped_file_list in enumerate(partition_list):
                if cpu_grouped_file_list and partition_index_list[loop_index]:
                    job_combinations.append((watershed.lower(),
                                             subbasin.lower(),
                                             cpu_grouped_file_list,
                                             partition_index_list[loop_index],
                                             weight_table_file,
                                             grid_type,
                                             master_rapid_runoff_file,
                                             RAPID_Inflow_Tool,
                                             mp_lock))
                    #COMMENTED CODE IS FOR DEBUGGING
##                    generate_inflows_from_runoff((watershed.lower(),
##                                                  subbasin.lower(),
##                                                  cpu_grouped_file_list,
##                                                  partition_index_list[loop_index],
##                                                  weight_table_file,
##                                                  grid_type,
##                                                  master_rapid_runoff_file,
##                                                  RAPID_Inflow_Tool,
##                                                  mp_lock))
            pool = multiprocessing.Pool(NUM_CPUS)
            pool.map(generate_inflows_from_runoff,
                     job_combinations)
            pool.close()
            pool.join()

            #run RAPID for the watershed
            lsm_rapid_output_file = os.path.join(master_watershed_output_directory,
                                                 'Qout_{0}'.format(out_file_ending))
            rapid_manager.update_parameters(rapid_connect_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                                            r'rapid_connect\.csv'),
                                            Vlat_file=master_rapid_runoff_file,
                                            riv_bas_id_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                                         r'riv_bas_id\.csv'),
                                            k_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                                r'k\.csv'),
                                            x_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                                r'x\.csv'),
                                            Qout_file=lsm_rapid_output_file
                                            )
                                            
            rapid_manager.update_reach_number_data()

            if generate_rapid_namelist_file:
                rapid_manager.generate_namelist_file(os.path.join(master_watershed_input_directory,
                                                                  "rapid_namelist_{}".format(out_file_ending[:-3])))
            if run_rapid_simulation:
                rapid_manager.run()

                rapid_manager.make_output_CF_compliant(simulation_start_datetime=actual_simulation_start_datetime,
                                                       comid_lat_lon_z_file=in_rivid_lat_lon_z_file,
                                                       project_name="{0} Based Historical flows by {1}".format(description,
                                                                                                               modeling_institution))

                #generate return periods
                if generate_return_periods_file and os.path.exists(lsm_rapid_output_file) and lsm_rapid_output_file:
                    return_periods_file = os.path.join(master_watershed_output_directory,
                                                       'return_periods_{0}'.format(out_file_ending))
                    #assume storm has 3 day length
                    storm_length_days = 3
                    generate_return_periods(lsm_rapid_output_file,
                                            return_periods_file,
                                            NUM_CPUS,
                                            storm_length_days)
                        
                if generate_seasonal_initialization_file and os.path.exists(lsm_rapid_output_file) and lsm_rapid_output_file:
                    seasonal_qinit_file = os.path.join(master_watershed_input_directory,
                                                       'seasonal_qinit_{0}.csv'.format(out_file_ending[:-3]))
                    rapid_manager.generate_seasonal_intitialization(seasonal_qinit_file)
    
                if generate_initialization_file and os.path.exists(lsm_rapid_output_file) and lsm_rapid_output_file:
                    qinit_file = os.path.join(master_watershed_input_directory,
                                              'qinit_{0}.csv'.format(out_file_ending[:-3]))
                    rapid_manager.generate_qinit_from_past_qout(qinit_file)


    #print info to user
    time_end = datetime.utcnow()
    print("Time Begin All: {0}".format(time_begin_all))
    print("Time Finish All: {0}".format(time_end))
    print("TOTAL TIME: {0}".format(time_end-time_begin_all))
