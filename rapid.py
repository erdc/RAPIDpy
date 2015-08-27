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
    #*******************************************************************************
    #Runtime options 
    #*******************************************************************************
    BS_opt_Qinit = False
    #!.false. --> no read initial flow    .true. --> read initial flow
    BS_opt_Qfinal = False
    #!.false. --> no write final flow     .true. --> write final flow 
    BS_opt_dam = False
    #!.false. --> no dam model used       .true. --> dam model used
    BS_opt_for = False
    #!.false. --> no forcing              .true. --> forcing
    BS_opt_influence = False
    #!.false. --> no output influence     .true. --> output influence
    IS_opt_routing = 1
    #!1       --> matrix-based Muskingum  2      --> traditional Muskingum
    #!3       --> Transbnd. matrix-based
    IS_opt_run = 1
    #!1       --> regular run             2      --> parameter optimization
    IS_opt_phi = 1
    #!1       --> phi1                    2      --> phi2
    
    #*******************************************************************************
    #Temporal information
    #*******************************************************************************
    #NOTE: ALL TIME IN SECONDS!
    #ALWAYS USED
    ZS_TauR = 0 #duration of routing procedure (time step of runoff data)
    ZS_dtR = 0 #internal routing time step
    #ONLY FOR REGULAR RUN
    ZS_TauM = 0 #total simulation time 
    ZS_dtM = 0 #input time step 
    #ONLY FOR OPTIMIZATION RUN
    ZS_TauO = 0 #total optimization time  
    ZS_dtO = 0 #observation time step
    #FORCING MODE (replace some values with observations) 
    ZS_dtF = 0 #time step of forcing data
    
    #*******************************************************************************
    #Domain in which input data is available
    #*******************************************************************************
    IS_riv_tot = 0 #number of river reaches in rapid connect file
    rapid_connect_file = '' #path to rapid_connect file
    IS_max_up = 2 #maximum number of ustream segments
    Vlat_file = '' #path to runoff file
    
    #*******************************************************************************
    #Domain in which model runs
    #*******************************************************************************
    IS_riv_bas = 0 #number of river reaches in subbasin
    riv_bas_id_file = '' #subbasin reach id file
    
    #*******************************************************************************
    #Initial instantaneous flow file
    #*******************************************************************************
    Qinit_file = '' #initial flow file (same order as rapid_connect)
    
    #*******************************************************************************
    #Final instantaneous flow file
    #*******************************************************************************
    Qfinal_file = '' #path to output final flow file
    
    #*******************************************************************************
    #Available dam data
    #*******************************************************************************
    IS_dam_tot = 0 #number of dams
    dam_tot_id_file = '' #ids of dam location
    
    #*******************************************************************************
    #Dam data used
    #*******************************************************************************
    IS_dam_use = 0 #number in subset of dam data to use
    dam_use_id_file = '' #ids of subset of dams
    
    #*******************************************************************************
    #Available forcing data
    #*******************************************************************************
    IS_for_tot = 0
    for_tot_id_file = ''
    Qfor_file = ''
    
    #*******************************************************************************
    #Forcing data used as model runs
    #*******************************************************************************
    IS_for_use = 0
    for_use_id_file = ''
    
    #*******************************************************************************
    #File where max (min) of absolute values of b (QoutR) are stored
    #*******************************************************************************
    babsmax_file = ''
    QoutRabsmin_file = ''
    QoutRabsmax_file = ''
    
    #*******************************************************************************
    #Regular model run
    #*******************************************************************************
    k_file = ''
    x_file = ''
    Qout_file = ''
    
    #*******************************************************************************
    #Optimization
    #*******************************************************************************
    ZS_phifac = 0
    #------------------------------------------------------------------------------
    #Routing parameters
    #------------------------------------------------------------------------------
    kfac_file = ''
    xfac_file = '' 
    ZS_knorm_init = 0
    ZS_xnorm_init = 0
    #------------------------------------------------------------------------------
    #Gage observations
    #------------------------------------------------------------------------------
    IS_obs_tot = 0
    obs_tot_id_file = ''
    Qobs_file = ''
    Qobsbarrec_file = ''
    IS_obs_use = 0
    obs_use_id_file = ''
    IS_strt_opt = 0

    
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
        print "Generating AutoRoute input file ..."
        try:
            os.remove(file_path)
        except OSError:
            pass
        
        with open(file_path,'w') as new_file:
            new_file.write('&NL_namelist\n')
            for attr, value in self.__dict__.iteritems():
                if not attr.startswith('_') and value:
                    new_file.write("%s = %s\n" % (attr, value))

        
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
        cv = ConvertRAPIDOutputToCF(self.Qfinal_file, #location of timeseries output file
                                   simulation_start_datetime, #time of the start of the simulation time
                                   self.ZS_TauR, #time step of simulation in seconds
                                   qinit_file=self.Qinit_file, #RAPID qinit file
                                   comid_lat_lon_z_file=comid_lat_lon_z_file, #path to comid_lat_lon_z file
                                   project_name=project_name, #name of your project
                                   output_id_dim_name='COMID', #name of ID dimension in output file, typically COMID or FEATUREID
                                   output_flow_var_name='Qout', #name of streamflow variable in output file, typically Qout or m3_riv
                                   print_debug=False)
        cv.convert()
        
        
    def run(self, rapid_namelist_file=""):
        """
        Run AutoRoute program and generate file based on inputs
        """
    
        time_start = datetime.datetime.utcnow()
    
        if not rapid_namelist_file or not os.path.exists(rapid_namelist_file):
            #generate input file if it does not exist
            if not rapid_namelist_file:
                rapid_namelist_file = "rapid_namelist"
            self.generate_namelist_file(rapid_namelist_file)
        else:
            #update existing file
            self.update_namelist_file(rapid_namelist_file)

        local_rapid_executable_location = os.path.join(os.path.basename(rapid_namelist_file), "rapid")
        #create link to RAPID
        os.symlink(self._rapid_executable_location, local_rapid_executable_location)


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
            run_rapid_command = ["mpiexec", "-n", self._num_processors]
        run_rapid_command.append(local_rapid_executable_location)
        run_rapid_command.append("-ksp_type")
        run_rapid_command.append("richardson")
                                 
        process = Popen(run_rapid_command, 
                        stdout=PIPE, stderr=PIPE, shell=False)
        out, err = process.communicate()
        if err:
            print err
            rapid_cleanup(local_rapid_executable_location, rapid_namelist_file)
            raise
        else:
            print 'RAPID output:'
            for line in out.split('\n'):
                print line

        rapid_cleanup(local_rapid_executable_location, rapid_namelist_file)
        print "Time to run RAPID: %s" % (datetime.datetime.utcnow()-time_start)


if __name__ == "__main__":
    input_folder = "/Users/rdchlads/autorapid/rapid-io/input/nfie_texas_gulf_region-huc_2_12"
    rapid_mng = RAPID('/Users/rdchlads/autorapid/rapid/run/rapid',
                      stream_file=os.path.join(input_folder, "streamflow_raster.tif"),
                      dem_file=os.path.join(input_folder, "Korea_DEMs", "merged_dems.tif"),
                      spatial_units="deg",
                      SHP_Out_File=os.path.join(input_folder,"tmp", "flood.tif"),
                      SHP_Out_Shapefile=os.path.join(input_folder,"tmp", "flood_Shp.shp"),
                     )
                         
    rapid_mng.run(autoroute_input_file=os.path.join(input_folder, "Flood", "AUTOROUTE_INPUT_FILE.txt"))
            
            
            