# RAPIDpy

Python scripting interface for the RAPID program.
More information about installation and the input parameters can be found at http://rapid-hub.org.
The source code for RAPID is located at https://github.com/c-h-david/rapid/.

[![DOI](https://zenodo.org/badge/19918/erdc-cm/RAPIDpy.svg)](https://zenodo.org/badge/latestdoi/19918/erdc-cm/RAPIDpy)

[![Build Status](https://travis-ci.org/erdc-cm/RAPIDpy.svg?branch=master)](https://travis-ci.org/erdc-cm/RAPIDpy)

[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-yellow.svg)](https://github.com/erdc-cm/RAPIDpy/blob/master/LICENSE)

#Installation

##Step 1: Install RAPID

### Before Installation Steps:
#### Ubuntu:
```
$ apt-get install gfortran g++ openmpi-bin
```
#### Windows with Cygwin:
Downloaded Cygwin (64-bit) (https://www.cygwin.com/) with these dependencies:
- dos2unix
- gcc-core
- gcc-fortran
- gcc-g++
- gdb
- git
- make
- openmpi
- time

### Installation Steps:
- Here is a script to download and install prereqs: http://rapid-hub.org/data/rapid_install_prereqs.tar.gz
- Follow the instructions on page 10-14: http://rapid-hub.org/docs/RAPID_Azure.pdf

##Step 2: Install netCDF4
###On Ubuntu:
```
$ apt-get install python-dev zlib1g-dev libhdf5-serial-dev libnetcdf-dev
```
###On Redhat/CentOS:
```
$ yum install epel-release
$ yum install netcdf4-python hdf5-devel netcdf-devel
```
###On OSX:
```
$ brew install homebrew/science/netcdf
```
##Step 3: Install RAPIDpy
```
$ sudo su
$ pip install netCDF4 RAPIDpy
$ exit
```
#How to use
##Running RAPID
### Step 1: Initialize the RAPID manager class. 
- First, add the path to you rapid executable location. 
- Next, you need to either tell it how many processors to use using the *num_processors* input variable or to use all available processors set *use_all_processors* to true.
- After that, add any other parameters you would like to use that would normally be in the rapid namelist file (this is case sensitive).


Example:
```python
from RAPIDpy.rapid import RAPID
rapid_manager = RAPID(rapid_executable_location='~/work/rapid/run/rapid'
                      use_all_processors=True,                          
                      ZS_TauR=24*3600, #duration of routing procedure (time step of runoff data)
                      ZS_dtR=15*60, #internal routing time step
                      ZS_TauM=365*24*3600, #total simulation time 
                      ZS_dtM=24*3600 #input time step 
                     )
```
If you are using Cygwin on Windows:
```python
from RAPIDpy.rapid import RAPID
rapid_manager = RAPID(rapid_executable_location='C:\\cygwin64\\home\\username\\work\\rapid\\run\\rapid',
                      cygwin_bin_location='C:\\cygwin64\\bin',
                      use_all_processors=True,                          
                      ZS_TauR=24*3600, #duration of routing procedure (time step of runoff data)
                      ZS_dtR=15*60, #internal routing time step
                      ZS_TauM=365*24*3600, #total simulation time 
                      ZS_dtM=24*3600 #input time step 
                     )
```

### Step 2 (optional): Add/update additional parameters later
You can add or update parameters using the *update_parameters* function by using the name of the variable in the rapid namelist file (this is case sensitive).


Example:
```python
rapid_manager.update_parameters(rapid_connect_file='../rapid_input_directory/rapid_connect.csv',
                                Vlat_file='../rapid_input_directory/m3_riv.nc',
                                riv_bas_id_file='../rapid_input_directory/riv_bas_id.csv,
                                k_file='../rapid_input_directory/k.csv',
                                x_file='../rapid_input_directory/x.csv',
                                Qout_file='../OUTPUT/Qout.nc'
                                )
```
### Step 3 (optional): Update reach number data
If you don't want to manually count the numbers for the rapid_connect or riv_bas_id files, use the *update_reach_number_data* function.


Example:
```python
rapid_manager.update_reach_number_data()
```

### Step 4: Run RAPID
This will generate your rapid_namelist file and run RAPID from wherever you call this script (your working directory).

Example:
```python
rapid_manager.run()
```

###Step 5 (optional): Convert RAPID output to be CF Compliant
This will convert the RAPID output to be CF compliant. This will require a comid_lat_lon_z file.
Additionally, it prepends time zero to you simulation. If no qinit file is given, a value of zero is added.

Example:
```python
rapid_manager.make_output_CF_compliant(simulation_start_datetime=datetime.datetime(1980, 1, 1),
                                       comid_lat_lon_z_file='../rapid_input_directory/comid_lat_lon_z.csv',
                                       project_name="ERA Interim Historical flows by US Army ERDC") 
```
##Getting USGS Daily Gage Data
You can use this to download USGS daily avg flow data. BE CAREFUL because overuse will get you blocked from downloading data.


Example reach_id_gage_id_file:
```
COMID, USGS_GAGE_ID
2000, 503944
...
```

Example:
```python
from os.path import join
main = "/home/username/data"
rapid_manager.generate_usgs_avg_daily_flows_opt(reach_id_gage_id_file=join(main,"mississippi_usgsgage_id_comid.csv"),
												start_datetime=datetime.datetime(2000,1,1),
												end_datetime=datetime.datetime(2014,12,31),
												out_streamflow_file=join(main,"streamflow_2000_2014.csv"), 
												out_stream_id_file=join(main,"streamid_2000_2014.csv"))
```
												
##Merge RAPID output
You can use this to combine consecutive RAPID output files into one file. WARNING: This code replaces the first file with the combined output and deletes the second file. BACK UP YOUR FILES!!!!

Example:
```python
import datetime
from RAPIDpy.make_CF_RAPID_output import ConvertRAPIDOutputToCF
file1 = ""
file2 = ""
cv = ConvertRAPIDOutputToCF(rapid_output_file=[file1, file2],
                            start_datetime=datetime.datetime(2005,1,1),
                            time_step=[3*3600, 3*3600],
                            qinit_file="",
                            comid_lat_lon_z_file="",
                            rapid_connect_file="",
                            project_name="NLDAS(VIC)-RAPID historical flows by US Army ERDC",
                            output_id_dim_name='COMID',
                            output_flow_var_name='Qout',
                            print_debug=False)
cv.convert()
```

##Write Qout timeseries to csv file
This function simplifies writing time series for each stream reach to a csv file.

```python
from RAPIDpy.rapid import write_flows_to_csv
path_to_file = '/p/work1/u4hfraat/rapid_projects/output_mississippi-nfie/Qout_k3v1_2005to2009.nc'
#for writing entire time series to file
write_flows_to_csv(path_to_file, 
		   reach_id=3624735) #COMID or rivid
		   
#if file is CF compliant, you can write out daily average
write_flows_to_csv(path_to_file, 
		   ind=20, #index of COMID or rivid (example if you already know index instead if reach_id)
		   daily=True) #if file is CF compliant, write out daily flows
```

##Generate qinit from past qout
RAPIDpy also creates a qinit file from a RAPID qout file. This example shows how.
```python
from RAPIDpy.rapid import RAPID
rapid_manager = RAPID(rapid_executable_location='/p/u4hfraat/rapid/src/rapid',
                      Qout_file='/rapid_projects/output_mississippi-nfie/Qout_k2v1_2005to2009.nc', 
		      rapid_connect_file='/autorapid/input_mississippi_nfie/rapid_connect_ECMWF.csv'
                     )
rapid_manager.generate_qinit_from_past_qout(qinit_file='/autorapid/input_mississippi_nfie/Qinit_2008_flood.csv',
					    time_index=10162) #time_index is optional, if not included it will be last time step
```
##Goodness of Fit
To check how well your simulation performed versus observations, this function can help you.
```python
from RAPIDpy.goodness_of_fit import find_goodness_of_fit_csv, find_goodness_of_fit

#if you have observations and model results next to each other in columns, this works for you
find_goodness_of_fit_csv('/Users/rdchlads/autorapid/rapid-io/output/united_kingdom-thames/flows_kingston_gage_noah.csv')

#if you have gage files in the format RAPID needs for ccalibration, this will generate a file for you with results
find_goodness_of_fit(reach_id_file='/p/work1/u4hfraat/rapid_projects/input_mississippi-nfie/all_gauges_nfie_id.csv',
                     rapid_qout_file='/p/work1/u4hfraat/rapid_projects/output_mississippi-nfie/Qout_k4v2_2005to2014.nc',
                     observed_file='/p/work1/u4hfraat/rapid_projects/input_mississippi-nfie/obs_all_gauges_nfie.csv',
                     out_analysis_file='/p/work1/u4hfraat/rapid_projects/output_mississippi-nfie/flows_analysis_vic_k4v2_2005-14.csv',
                     #steps_per_group=8, #for raw rapid output (8 is produces daily flows for 3-hr timesteps)
                     daily=True) #use this option for CF compliant files
