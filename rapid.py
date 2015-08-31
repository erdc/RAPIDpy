# -*- coding: utf-8 -*-

import datetime
from multiprocessing import cpu_count
import os
from subprocess import Popen, PIPE

#local
from helper_functions import csv_to_list
from make_CF_RAPID_output import ConvertRAPIDOutputToCF

#------------------------------------------------------------------------------
#Main Dataset Manager Class
#------------------------------------------------------------------------------
class RAPID(object):
    """
    This class is designed to prepare the rapid_namelist file and run 
    the RAPID program.
    """
    def __init__(self, rapid_executable_location, num_processors=1, 
                 use_all_processors=False, **kwargs):
        """
        Initialize the class with variables given by the user
        """
        self._rapid_executable_location = rapid_executable_location
        self._num_processors = num_processors
        #use all processors akes precedent over num_processors arg
        if use_all_processors == True:
            self._num_processors = cpu_count()
            
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
        

    def update_parameters(self, **kwargs):
        """
        Update AutoRoute parameters
        """
        #set arguments based off of user input
        for key, value in kwargs.iteritems():
            if key in dir(self) and not key.startswith('_'):
                setattr(self, key, value)
            else:
                raise Exception("Invalid RAPID parameter %s." % key)
    
    def update_reach_number_data(self):
        """
        Updates the reach number data based on input files
        """
        #get rapid connect info
        rapid_connect_table = csv_to_list(self.rapid_connect_file)
        self.IS_riv_tot = len(rapid_connect_table)
        self.IS_max_up = max([int(float(row[2])) for row in rapid_connect_table])
    
        #get riv_bas_id info
        riv_bas_id_table = csv_to_list(self.riv_bas_id_file)
        self.IS_riv_bas = len(riv_bas_id_table)


    def generate_namelist_file(self, file_path):
        """
        Generate rapid_namelist file
        """
        print "Generating RAPID namelist file ..."
        try:
            os.remove(file_path)
        except OSError:
            pass
        
        with open(file_path,'w') as new_file:
            new_file.write('&NL_namelist\n')
            for attr, value in self.__dict__.iteritems():
                if not attr.startswith('_') and value:
                    if attr.startswith('BS'):
                        new_file.write("%s = .%s.\n" % (attr, str(value).lower()))
                    elif isinstance(value, int):
                        new_file.write("%s = %s\n" % (attr, value))
                    else:
                        new_file.write("%s = \'%s\'\n" % (attr, value))
            new_file.write("/\n")
        
    def update_namelist_file(self, file_path):
        """
        Update existing namelist file with new parameters
        """
        if os.path.exists(file_path) and file_path:
            print "Adding missing inputs from RAPID input file ..."
            old_file = open(file_path, 'r')
            for line in old_file:
                line = line.strip()
                if not line[:1].isalpha() or not line:
                    continue
                line_split = line.split()
                attr = line_split[0]
                value = None
                if len(line_split)>1:
                    value = line_split[1]
                elif attr in self._no_value_attr_list:
                    value = True
                #add attribute if exists
                if attr in dir(self) \
                    and not attr.startswith('_'):
                    #set attribute if not set already
                    if not getattr(self, attr):
                        setattr(self, attr, value)
                else:
                    print "Invalid argument" , attr, ". Skipping ..."
            old_file.close()
            
            self.generate_input_file(file_path)
        else:
            raise Exception("RAPID namelist file to update not found.")
            
    def make_output_CF_compliant(self, 
                                 simulation_start_datetime,
                                 comid_lat_lon_z_file="",
                                 project_name="Normal RAPID project"):
        """
        Converts RAPID output to be CF compliant
        """
        print self.Qout_file
        cv = ConvertRAPIDOutputToCF(rapid_output_file=self.Qout_file, #location of timeseries output file
                                    start_datetime=simulation_start_datetime, #time of the start of the simulation time
                                    time_step=self.ZS_TauR, #time step of simulation in seconds
                                    qinit_file=self.Qinit_file, #RAPID qinit file
                                    comid_lat_lon_z_file=comid_lat_lon_z_file, #path to comid_lat_lon_z file
                                    rapid_connect_file=self.rapid_connect_file, #path to RAPID connect file
                                    project_name=project_name, #name of your project
                                    output_id_dim_name='COMID', #name of ID dimension in output file, typically COMID or FEATUREID
                                    output_flow_var_name='Qout', #name of streamflow variable in output file, typically Qout or m3_riv
                                    print_debug=False)
        cv.convert()
        
        
    def run(self, rapid_namelist_file=""):
        """
        Run RAPID program and generate file based on inputs
        """
    
        time_start = datetime.datetime.utcnow()
    
        if not rapid_namelist_file or not os.path.exists(rapid_namelist_file):
            #generate input file if it does not exist
            if not rapid_namelist_file:
                rapid_namelist_file = os.path.join(os.getcwd(), "rapid_namelist")
            self.generate_namelist_file(rapid_namelist_file)
        else:
            #update existing file
            self.update_namelist_file(rapid_namelist_file)

        local_rapid_executable_location = os.path.join(os.path.dirname(rapid_namelist_file), "rapid")

        #create link to RAPID if needed
        temp_link_to_rapid = ""
        if not os.path.exists(local_rapid_executable_location) \
            and not self._rapid_executable_location == local_rapid_executable_location:
            os.symlink(self._rapid_executable_location, local_rapid_executable_location)
            temp_link_to_rapid = local_rapid_executable_location

        def rapid_cleanup(local_rapid_executable, rapid_namelist_file):
            """
            Cleans up the rapid files generated by the process
            """
            #remove rapid link
            try:
                os.unlink(local_rapid_executable)
                os.remove(local_rapid_executable)
            except OSError:
                pass
    
            #remove namelist file
            try:
                os.remove(rapid_namelist_file)
            except OSError:
                pass
        
        #run RAPID
        print "Running RAPID ..."
        run_rapid_command = []
        if self._num_processors > 1:
            run_rapid_command = ["mpiexec", "-n", str(self._num_processors)]
        run_rapid_command.append(local_rapid_executable_location)
        run_rapid_command.append("-ksp_type")
        run_rapid_command.append("richardson")
        
        process = Popen(run_rapid_command, 
                        stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            print err
            rapid_cleanup(temp_link_to_rapid, rapid_namelist_file)
            raise
        else:
            print 'RAPID output:'
            for line in out.split('\n'):
                print line

        rapid_cleanup(temp_link_to_rapid, rapid_namelist_file)
        print "Time to run RAPID: %s" % (datetime.datetime.utcnow()-time_start)


"""
if __name__ == "__main__":
    rapid_manager = RAPID(rapid_executable_location=rapid_executable_location,
                          use_all_processors=True,                          
                          ZS_TauR = 24*3600, #duration of routing procedure (time step of runoff data)
                          ZS_dtR = 15*60, #internal routing time step
                          ZS_TauM = len(era_interim_file_list)*24*3600, #total simulation time 
                          ZS_dtM = 24*3600 #input time step 
                         )
    era_rapid_output_file = os.path.join(master_watershed_output_directory,
                                                           'Qout_erai.nc')
    rapid_manager.update_parameters(rapid_connect_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                                 r'rapid_connect\.csv'),
                                    Vlat_file=master_rapid_runoff_file,
                                    riv_bas_id_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                                 r'riv_bas_id\.csv'),
                                    k_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                        r'k\.csv'),
                                    x_file=case_insensitive_file_search(master_watershed_input_directory,
                                                                        r'x\.csv'),
                                    Qout_file=era_rapid_output_file
                                    )

    comid_lat_lon_z_file = case_insensitive_file_search(master_watershed_input_directory,
                                                        r'comid_lat_lon_z\.csv')

    rapid_manager.update_reach_number_data()
    rapid_manager.run()
    rapid_manager.make_output_CF_compliant(simulation_start_datetime=datetime.datetime(1980, 1, 1),
                                           comid_lat_lon_z_file=comid_lat_lon_z_file,
                                           project_name="ERA Interim Historical flows by US Army ERDC")     
"""
            
            