# -*- coding: utf-8 -*-
##
##  rapid.py
##  RAPIDpy
##
##  Created by Alan D Snow, 2015.
##  Copyright © 2015 Alan D Snow. All rights reserved.
##

from csv import writer as csvwriter
import datetime
from dateutil.parser import parse
from dateutil.tz import tzoffset
from multiprocessing import cpu_count
import numpy as np
import os
#USGS not returning tzinfo, so this is no longer needed
#from pytz import utc
from requests import get
from subprocess import Popen, PIPE

#local
from .dataset import RAPIDDataset
from .helper_functions import csv_to_list, log, open_csv
from .postprocess import ConvertRAPIDOutputToCF

#------------------------------------------------------------------------------
#Main Dataset Manager Class
#------------------------------------------------------------------------------
class RAPID(object):
    """
    This class is designed to prepare the rapid_namelist file and run 
    the RAPID program.
    """
    def __init__(self, rapid_executable_location="", num_processors=1, 
                 use_all_processors=False, cygwin_bin_location="",
                 mpiexec_command="mpiexec", ksp_type="richardson",
                 **kwargs):
        """
        Initialize the class with variables given by the user
        """
        self._rapid_executable_location = rapid_executable_location
        if os.name == "nt" and (not cygwin_bin_location or not os.path.exists(cygwin_bin_location)):
            raise Exception("Required to have cygwin_bin_location set if using windows!")
        self._cygwin_bin_location = cygwin_bin_location
        self._cygwin_bash_exe_location = os.path.join(cygwin_bin_location, "bash.exe")
        self._mpiexec_command = mpiexec_command
        self._ksp_type = ksp_type
        
        #use all processors makes precedent over num_processors arg
        if use_all_processors == True:
            self._num_processors = cpu_count()
        elif num_processors > cpu_count():
            log("Num processors requested exceeded max. Set to max ...",
                "WARNING")
            self._num_processors = cpu_count()
        else:
            self._num_processors = num_processors
    
        #*******************************************************************************
        #Runtime options 
        #*******************************************************************************
        self.BS_opt_Qinit = False
        #!.false. --> no read initial flow    .true. --> read initial flow
        self.BS_opt_Qfinal = False
        #!.false. --> no write final flow     .true. --> write final flow 
        self.BS_opt_dam = False
        #!.false. --> no dam model used       .true. --> dam model used
        self.BS_opt_for = False
        #!.false. --> no forcing              .true. --> forcing
        self.BS_opt_influence = False
        #!.false. --> no output influence     .true. --> output influence
        self.IS_opt_routing = 1
        #!1       --> matrix-based Muskingum  2      --> traditional Muskingum
        #!3       --> Transbnd. matrix-based
        self.IS_opt_run = 1
        #!1       --> regular run             2      --> parameter optimization
        self.IS_opt_phi = 1
        #!1       --> phi1                    2      --> phi2
        
        #*******************************************************************************
        #Temporal information
        #*******************************************************************************
        #NOTE: ALL TIME IN SECONDS!
        #ALWAYS USED
        self.ZS_TauR = 0 #duration of routing procedure (time step of runoff data)
        self.ZS_dtR = 0 #internal routing time step
        #ONLY FOR REGULAR RUN
        self.ZS_TauM = 0 #total simulation time 
        self.ZS_dtM = 0 #input time step 
        #ONLY FOR OPTIMIZATION RUN
        self.ZS_TauO = 0 #total optimization time  
        self.ZS_dtO = 0 #observation time step
        #FORCING MODE (replace some values with observations) 
        self.ZS_dtF = 0 #time step of forcing data
        
        #*******************************************************************************
        #Domain in which input data is available
        #*******************************************************************************
        self.IS_riv_tot = 0 #number of river reaches in rapid connect file
        self.rapid_connect_file = '' #path to rapid_connect file
        self.IS_max_up = 2 #maximum number of ustream segments
        self.Vlat_file = '' #path to runoff file
        
        #*******************************************************************************
        #Domain in which model runs
        #*******************************************************************************
        self.IS_riv_bas = 0 #number of river reaches in subbasin
        self.riv_bas_id_file = '' #subbasin reach id file
        
        #*******************************************************************************
        #Initial instantaneous flow file
        #*******************************************************************************
        self.Qinit_file = '' #initial flow file (same order as rapid_connect)
        
        #*******************************************************************************
        #Final instantaneous flow file
        #*******************************************************************************
        self.Qfinal_file = '' #path to output final flow file
        
        #*******************************************************************************
        #Available dam data
        #*******************************************************************************
        self.IS_dam_tot = 0 #number of dams
        self.dam_tot_id_file = '' #ids of dam location
        
        #*******************************************************************************
        #Dam data used
        #*******************************************************************************
        self.IS_dam_use = 0 #number in subset of dam data to use
        self.dam_use_id_file = '' #ids of subset of dams
        
        #*******************************************************************************
        #Available forcing data
        #*******************************************************************************
        self.IS_for_tot = 0
        self.for_tot_id_file = ''
        self.Qfor_file = ''
        
        #*******************************************************************************
        #Forcing data used as model runs
        #*******************************************************************************
        self.IS_for_use = 0
        self.for_use_id_file = ''
        
        #*******************************************************************************
        #File where max (min) of absolute values of b (QoutR) are stored
        #*******************************************************************************
        self.babsmax_file = ''
        self.QoutRabsmin_file = ''
        self.QoutRabsmax_file = ''
        
        #*******************************************************************************
        #Regular model run
        #*******************************************************************************
        self.k_file = ''
        self.x_file = ''
        self.Qout_file = ''
        
        #*******************************************************************************
        #Optimization
        #*******************************************************************************
        self.ZS_phifac = 0
        #------------------------------------------------------------------------------
        #Routing parameters
        #------------------------------------------------------------------------------
        self.kfac_file = ''
        self.xfac_file = '' 
        self.ZS_knorm_init = 0
        self.ZS_xnorm_init = 0
        #------------------------------------------------------------------------------
        #Gage observations
        #------------------------------------------------------------------------------
        self.IS_obs_tot = 0
        self.obs_tot_id_file = ''
        self.Qobs_file = ''
        self.Qobsbarrec_file = ''
        self.IS_obs_use = 0
        self.obs_use_id_file = ''
        self.IS_strt_opt = 0
        
        
        
        self.update_parameters(**kwargs)
        

    def _get_cygwin_path(self, windows_path):
        """
        Convert windows path to cygpath
        """
        conv_cmd = [os.path.join(self._cygwin_bin_location, "cygpath.exe"),
                    "-u", windows_path]
        process = Popen(conv_cmd, 
                        stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            print(err)
            raise Exception(err)
        
        return out.strip()

    def _create_symlink_cygwin(self, initial_path, final_path):
        """
        Use cygqin to generate symbolic link
        """
        symlink_cmd = [os.path.join(self._cygwin_bin_location, "ln.exe"),
                       "-s", self._get_cygwin_path(initial_path), 
                       self._get_cygwin_path(final_path)]
        process = Popen(symlink_cmd, 
                        stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            print(err)
            raise Exception(err)
        
        return out.strip()

    def _dos2unix_cygwin(self, file_path):
        """
        Use cygwin to convert file to unix format
        """
        dos2unix_cmd = [os.path.join(self._cygwin_bin_location, "dos2unix.exe"),
                       self._get_cygwin_path(file_path)]
        process = Popen(dos2unix_cmd, 
                        stdout=PIPE, stderr=PIPE, shell=False)
        process.communicate()
    
    def update_parameters(self, **kwargs):
        """
        Update RAPID parameters
        """
        #set arguments based off of user input
        for key, value in list(kwargs.items()):
            if key in dir(self) and not key.startswith('_'):
                setattr(self, key, value)
            else:
                log("Invalid RAPID parameter %s." % key,
                    "ERROR")
    
    def update_reach_number_data(self):
        """
        Updates the reach number data based on input files
        """
        
        if not self.rapid_connect_file or not self.rapid_connect_file:
            log("Missing rapid_connect_file. Please set before running this function ...",
                "ERROR")

        if not self.riv_bas_id_file or not self.riv_bas_id_file:
            log("Missing riv_bas_id_file. Please set before running this function ...",
                "ERROR")

        #get rapid connect info
        rapid_connect_table = csv_to_list(self.rapid_connect_file)
        self.IS_riv_tot = len(rapid_connect_table)
        self.IS_max_up = max([int(float(row[2])) for row in rapid_connect_table])
    
        #get riv_bas_id info
        riv_bas_id_table = csv_to_list(self.riv_bas_id_file)
        self.IS_riv_bas = len(riv_bas_id_table)
    
    def update_simulation_runtime(self):
        """
        Updates the total simulation runtime from
        the m3 file and the time step
        """
        if not self.Vlat_file or not os.path.exists(self.Vlat_file):
            log("Need Vlat_file to proceed ...",
                "ERROR")

        if self.ZS_TauR <= 0:
            log("Missing routing time step ...",
                "ERROR")
    
        try:
            self.ZS_TauR = int(self.ZS_TauR)
        except Exception:
            log("Invalid routing time step: {0} ...".format(self.ZS_TauR),
                "ERROR")

        with RAPIDDataset(self.Vlat_file) as m3_nc:
            self.ZS_TauM = m3_nc.size_time*self.ZS_TauR
            self.ZS_TauO = m3_nc.size_time*self.ZS_TauR

    def generate_namelist_file(self, file_path):
        """
        Generate rapid_namelist file
        """
        log("Generating RAPID namelist file ...",
            "INFO")
        try:
            os.remove(file_path)
        except OSError:
            pass
        
        with open(file_path,'w') as new_file:
            new_file.write('&NL_namelist\n')
            for attr, value in sorted(list(self.__dict__.items())):
                if not attr.startswith('_'):
                    if attr.startswith('BS'):
                        new_file.write("%s = .%s.\n" % (attr, str(value).lower()))
                    elif isinstance(value, int):
                        new_file.write("%s = %s\n" % (attr, value))
                    else:
                        if value:
                            #file path
                            if os.name == "nt":
                                #if windows generate file with cygpath
                                value = self._get_cygwin_path(value)
                            new_file.write("%s = \'%s\'\n" % (attr, value))
            new_file.write("/\n")
        
    def update_namelist_file(self, file_path):
        """
        Update existing namelist file with new parameters
        """
        if os.path.exists(file_path) and file_path:
            log("Adding missing inputs from RAPID input file ...",
                "INFO")
            old_file = open(file_path, 'r')
            for line in old_file:
                line = line.strip()
                if not line[:1].isalpha() or not line:
                    continue
                line_split = line.split("=")
                attr = line_split[0].strip()
                value = None
                if len(line_split)>1:
                    value = line_split[1].strip().replace("'", "").replace('"', "")
                    #convert integers to integers
                    try:
                        value = int(value)
                    except Exception:
                        pass
                    #remove dots from beginning & end of value
                    if attr.startswith('BS'):
                        value = value.replace(".", "") 
                elif attr in self._no_value_attr_list:
                    value = True
                #add attribute if exists
                if attr in dir(self) \
                    and not attr.startswith('_'):
                    #set attribute if not set already
                    if not getattr(self, attr):
                        setattr(self, attr, value)
                else:
                    log("Invalid argument {0}. Skipping ...".format(attr),
                        "INFO")
            old_file.close()
            
            self.generate_namelist_file(file_path)
        else:
            log("RAPID namelist file to update not found.",
                "ERROR")
            
    def make_output_CF_compliant(self, 
                                 simulation_start_datetime,
                                 comid_lat_lon_z_file="",
                                 project_name="Normal RAPID project"):
        """
        Converts RAPID output to be CF compliant
        """
        need_to_convert = True
        with RAPIDDataset(self.Qout_file) as qout_nc:
            need_to_convert = not qout_nc.is_time_variable_valid()
        if not need_to_convert:
            log("RAPID Qout file already CF compliant ...",
                "INFO")
        else:
            cv = ConvertRAPIDOutputToCF(rapid_output_file=self.Qout_file, #location of timeseries output file
                                        start_datetime=simulation_start_datetime, #time of the start of the simulation time
                                        time_step=self.ZS_TauR, #time step of simulation in seconds
                                        qinit_file=self.Qinit_file, #RAPID qinit file
                                        comid_lat_lon_z_file=comid_lat_lon_z_file, #path to comid_lat_lon_z file
                                        rapid_connect_file=self.rapid_connect_file, #path to RAPID connect file
                                        project_name=project_name, #name of your project
                                        output_id_dim_name='rivid', #name of ID dimension in output file, typically COMID or FEATUREID
                                        output_flow_var_name='Qout', #name of streamflow variable in output file, typically Qout or m3_riv
                                        print_debug=False)
            cv.convert()
    
    def run(self, rapid_namelist_file=""):
        """
        Run RAPID program and generate file based on inputs
        """
    
        if not self._rapid_executable_location or not self._rapid_executable_location:
            log("Missing rapid_executable_location. Please set before running this function ...",
                "ERROR")

        time_start = datetime.datetime.utcnow()
    
        if not rapid_namelist_file or not os.path.exists(rapid_namelist_file):
            #generate input file if it does not exist
            if not rapid_namelist_file:
                rapid_namelist_file = os.path.join(os.getcwd(), "rapid_namelist")
            self.generate_namelist_file(rapid_namelist_file)
        else:
            #update existing file
            self.update_namelist_file(rapid_namelist_file)

        local_rapid_executable_location = os.path.join(os.path.dirname(rapid_namelist_file), "rapid_exe_symlink")

        def rapid_cleanup(*args):
            """
            Cleans up the rapid files generated by the process
            """
            for arg in args:
                #remove files
                try:
                    os.remove(arg)
                except OSError:
                    pass

        #create link to RAPID if needed
        temp_link_to_rapid = ""
        if not self._rapid_executable_location == local_rapid_executable_location:
            rapid_cleanup(local_rapid_executable_location)
            if os.name == "nt":
                self._create_symlink_cygwin(self._rapid_executable_location, 
                                            local_rapid_executable_location)
            else:
                os.symlink(self._rapid_executable_location, local_rapid_executable_location)
            temp_link_to_rapid = local_rapid_executable_location

        
        #run RAPID
        log("Running RAPID ...",
            "INFO")
        run_rapid_script = ""
        if os.name == "nt":
            run_rapid_script = os.path.join(os.getcwd(), "run_rapid.sh")
            with open(run_rapid_script, "w") as run_rapid:
                run_rapid.write("#!/bin/sh\n")
                run_rapid.write("cd {}\n".format(self._get_cygwin_path(os.getcwd())))
                if self._num_processors > 1:
                    run_rapid.write("{0} -np {1} {2} -ksp_type {3}\n".format(self._mpiexec_command,
                                                                             self._num_processors,
                                                                             self._get_cygwin_path(local_rapid_executable_location),
                                                                             self._ksp_type))
                else:
                    #htcondor will not allow mpiexec for single processor jobs
                    #this was added for that purpose
                    run_rapid.write("{0} -ksp_type {1}\n".format(self._get_cygwin_path(local_rapid_executable_location),
                                                                 self._ksp_type))
                
            
            self._dos2unix_cygwin(run_rapid_script)
            run_rapid_command = [self._cygwin_bash_exe_location, "-l", "-c", 
                                 self._get_cygwin_path(run_rapid_script)]

        else:
            #htcondor will not allow mpiexec for single processor jobs
            #this was added for that purpose
            run_rapid_command = [local_rapid_executable_location, 
                                 "-ksp_type", self._ksp_type]
                                 
            if self._num_processors > 1:
                run_rapid_command = [self._mpiexec_command, "-n", str(self._num_processors),
                                     local_rapid_executable_location, 
                                     "-ksp_type", self._ksp_type]

        process = Popen(run_rapid_command, 
                        stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            rapid_cleanup(temp_link_to_rapid, rapid_namelist_file, run_rapid_script)
            raise Exception(err)
        else:
            log('RAPID output:',
                "INFO")
            for line in out.split(b'\n'):
                print(line)
        rapid_cleanup(temp_link_to_rapid, rapid_namelist_file, run_rapid_script)
        log("Time to run RAPID: %s" % (datetime.datetime.utcnow()-time_start),
            "INFO")

    def generate_qinit_from_past_qout(self, qinit_file, time_index=-1):
        """
        Generate qinit from qout file
        """
        if not self.Qout_file or not os.path.exists(self.Qout_file):
            log('Missing Qout_file. Please set before running this function ...',
                "ERROR")

        if not self.rapid_connect_file or not self.rapid_connect_file:
            log('Missing rapid_connect file. Please set before running this function ...',
                "ERROR")
 
        log("Generating qinit file from qout file ...",
            "INFO")
        #get information from dataset
        with RAPIDDataset(self.Qout_file) as qout_nc:
            log("Extracting data ...",
                "INFO")
            streamflow_values = qout_nc.get_qout(time_index=time_index)
    
            log("Reordering data ...",
                "INFO")
            rapid_connect_array = csv_to_list(self.rapid_connect_file)
            stream_id_array = np.array([row[0] for row in rapid_connect_array], np.int32)
            init_flows_array = np.zeros(len(rapid_connect_array))
            for riv_bas_index, riv_bas_id in enumerate(qout_nc.get_river_id_array()):
                try:
                    data_index = np.where(stream_id_array==riv_bas_id)[0][0]
                    init_flows_array[data_index] = streamflow_values[riv_bas_index]
                except Exception:
                    log('riv bas id {0} not found in connectivity list.'.format(riv_bas_id),
                        "WARNING")
        
        log("Writing to file ...",
            "INFO")
        with open_csv(qinit_file, 'w') as qinit_out:
            for init_flow in init_flows_array:
                qinit_out.write('{}\n'.format(init_flow))

        self.Qinit_file = qinit_file
        self.BS_opt_Qinit = True
        log("Initialization Complete!",
            "INFO")

    def generate_seasonal_intitialization(self, qinit_file,
                                          datetime_start_initialization=datetime.datetime.utcnow()):
        """
        This function loops through a CF compliant rapid streamflow
        file to produce estimates for current streamflow based on
        the seasonal average over the data within the historical streamflow
        file.
        """
        #get information from datasets
        if not self.Qout_file or not os.path.exists(self.Qout_file):
            log("Missing Qout_file. Please set before running this function ...",
                "ERROR")

        if not self.rapid_connect_file or not self.rapid_connect_file:
            log("Missing rapid_connect file. Please set before running this function ...",
                "ERROR")
        
        with RAPIDDataset(self.Qout_file) as qout_hist_nc:
            if not qout_hist_nc.is_time_variable_valid():
                log("File must be CF 1.6 compliant with valid time variable ...",
                    "ERROR")

            log("Generating seasonal average qinit file from qout file ...",
                "INFO")
            
            log("Determining dates with streamflows of interest ...",
                "INFO")

            datetime_min = datetime_start_initialization - datetime.timedelta(3)
            datetime_max = datetime_start_initialization + datetime.timedelta(3)
            
            time_indices = []
            for idx, t in enumerate(qout_hist_nc.get_time_array()):
                var_time = datetime.datetime.utcfromtimestamp(t)
                #check if date within range of season
                if var_time.month >= datetime_min.month and var_time.month <= datetime_max.month:
                    if var_time.month > datetime_min.month:
                        if var_time.day < datetime_max.day:
                            time_indices.append(idx)
                    elif var_time.day >= datetime_min.day and var_time.day < datetime_max.day:
                        time_indices.append(idx)

            if not time_indices:
                log("No time steps found within range ...",
                    "ERROR")
            
            log("Extracting data ...",
                "INFO")
            streamflow_array = qout_hist_nc.get_qout(time_index_array=time_indices)

            log("Reordering data...",
                "INFO")
            rapid_connect_array = csv_to_list(self.rapid_connect_file)
            stream_id_array = np.array([row[0] for row in rapid_connect_array], dtype=np.int32)
            init_flows_array = np.zeros(len(rapid_connect_array))
            for riv_bas_index, riv_bas_id in enumerate(qout_hist_nc.get_river_id_array()):
                try:
                    data_index = np.where(stream_id_array==riv_bas_id)[0][0]
                    init_flows_array[data_index] = np.mean(streamflow_array[riv_bas_index])
                except Exception:
                    log('riv_bas_id {0} not found in connectivity list.'.format(riv_bas_id),
                        "WARNING")

            log("Writing to file ...",
                "INFO")
            with open_csv(qinit_file, 'w') as qinit_out:
                for init_flow in init_flows_array:
                    qinit_out.write('{}\n'.format(init_flow))

            log("Initialization Complete!",
                "INFO")

    def generate_usgs_avg_daily_flows_opt(self, reach_id_gage_id_file,
                                          start_datetime, end_datetime,
                                          out_streamflow_file, out_stream_id_file):
        """
        Generate streamflow file and stream id file required for calibration 
        based on usgs gage ids associated with stream ids
        """
        log("Generating avg streamflow file and stream id file required for calibration ...",
            "INFO")
        reach_id_gage_id_list = csv_to_list(reach_id_gage_id_file)
# USGS not returning tzinfo anymore, so removed tzinfo operations
#       if start_datetime.tzinfo is None or start_datetime.tzinfo.utcoffset(start_datetime) is None:
#            start_datetime = start_datetime.replace(tzinfo=utc)
#        if end_datetime.tzinfo is None or end_datetime.tzinfo.utcoffset(end_datetime) is None:
#            end_datetime = end_datetime.replace(tzinfo=utc)
        gage_data_matrix = []
        valid_comid_list = []
        
        #add extra day as it includes the start date (e.g. 7-5 is 2 days, but have data for 5,6,7, so +1)
        num_days_needed = (end_datetime-start_datetime).days + 1

        gage_id_list = []
        for row in reach_id_gage_id_list[1:]:
            station_id = row[1]
            if len(row[1]) == 7:
                station_id = '0' + row[1]
            gage_id_list.append(station_id)
        
        num_gage_id_list = np.array(gage_id_list, dtype=np.int32)
        log("Querying Server for Data ..." ,
            "INFO")
    
        query_params = {
                        'format': 'json',
                        'sites': ",".join(gage_id_list),
# USGS not returning tzinfo anymore, so removed tzinfo operations 
#                        'startDT': start_datetime.astimezone(tzoffset(None, -18000)).strftime("%Y-%m-%d"),
#                        'endDT': end_datetime.astimezone(tzoffset(None, -18000)).strftime("%Y-%m-%d"),
                        'startDT': start_datetime.strftime("%Y-%m-%d"),
                        'endDT': end_datetime.strftime("%Y-%m-%d"),
                        'parameterCd': '00060', #streamflow
                        'statCd': '00003' #average
                       }
        response = get("http://waterservices.usgs.gov/nwis/dv", params=query_params)
        if response.ok:
            data_valid = True
            try:
                requested_data = response.json()['value']['timeSeries']
            except IndexError:
                data_valid = False
                pass
            
            if data_valid:
                for time_series in enumerate(requested_data):
                    usgs_station_full_name = time_series[1]['name']
                    usgs_station_id = usgs_station_full_name.split(":")[1]
                    gage_data = []
                    for time_step in time_series[1]['values'][0]['value']:
                        local_datetime = parse(time_step['dateTime'])
                        if local_datetime > end_datetime:
                            break
                        
                        if local_datetime >= start_datetime:
                            if not time_step['value']:
                                log("MISSING DATA for USGS Station {0} {1} {2}".format(station_id,
                                                                                       local_datetime,
                                                                                       time_step['value']),
                                    "WARNING")
                            gage_data.append(float(time_step['value'])/35.3146667)
    
                    try:
                        #get where streamids assocated with USGS sation id is
                        streamid_index = np.where(num_gage_id_list==int(float(usgs_station_id)))[0][0]+1
                    except Exception:
                        log("USGS Station {0} not found in list ...".format(usgs_station_id),
                            "WARNING")
                        raise
                        
                    if len(gage_data) == num_days_needed:
                        gage_data_matrix.append(gage_data)
                        valid_comid_list.append(reach_id_gage_id_list[streamid_index][0])
                    else:
                        log("StreamID {0} USGS Station {1} MISSING {2} DATA VALUES".format(reach_id_gage_id_list[streamid_index][0],
                                                                                           usgs_station_id,
                                                                                           num_days_needed-len(gage_data)),
                            "WARNING")

            if gage_data_matrix and valid_comid_list:
                log("Writing Output ...",
                    "INFO")
                np_array = np.array(gage_data_matrix).transpose()  
                with open_csv(out_streamflow_file, 'w') as gage_data:
                    wf = csvwriter(gage_data)
                    for row in np_array:
                        wf.writerow(row)
                        
                with open_csv(out_stream_id_file, 'w') as comid_data:
                    cf = csvwriter(comid_data)
                    for row in valid_comid_list:
                        cf.writerow([int(float(row))])
                        
                #set parameters for RAPID run
                self.IS_obs_tot = len(valid_comid_list)
                self.obs_tot_id_file = out_stream_id_file
                self.Qobs_file = out_streamflow_file
                self.IS_obs_use = len(valid_comid_list)
                self.obs_use_id_file = out_stream_id_file
            else:
                log("No valid data returned ...",
                    "WARNING")
        else:
                log("USGS query error ...",
                    "WARNING")
