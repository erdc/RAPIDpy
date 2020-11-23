#!/usr/bin/env python

from RAPIDpy.inflow import run_lsm_rapid_process
from datetime import datetime
import yaml
import sys

def get_params(filename):
    with open(filename, 'r') as f:
        p = yaml.load(f)
    return p

if __name__=='__main__':
    yaml_file = sys.argv[1]
    p = get_params(yaml_file)
   
    simulation_start_datetime = datetime.strptime(
        p['simulation_start_timestamp'], '%Y%m%d')
    simulation_end_datetime = datetime.strptime(
        p['simulation_end_timestamp'], '%Y%m%d')

    run_lsm_rapid_process(
	rapid_executable_location=p['rapid_executable_location'],
	lsm_data_location=p['lsm_data_location'],
	rapid_io_files_location=p['rapid_io_files_location'],
	rapid_input_location=p['rapid_input_location'],
	rapid_output_location=p['rapid_output_location'],
	simulation_start_datetime=simulation_start_datetime,
	simulation_end_datetime=simulation_end_datetime,
	file_datetime_pattern=p['file_datetime_pattern'],
	file_datetime_re_pattern=p['file_datetime_re_pattern'],
	initial_flows_file=p['initial_flows_file'],
	generate_rapid_namelist_file=p['generate_rapid_namelist_file'],
	run_rapid_simulation=p['run_rapid_simulation'],
	generate_return_periods_file=p['generate_return_periods_file'],
	return_period_method=p['return_period_method'],
	generate_seasonal_averages_file=p['generate_seasonal_averages_file'],
	generate_seasonal_initialization_file=p['generate_seasonal_averages_file'],
	generate_initialization_file=p['generate_initialization_file'],
	use_all_processors=p['use_all_processors'],
	num_processors=p['num_processors'],
	mpiexec_command=p['mpiexec_command'],
	cygwin_bin_location=p['cygwin_bin_location'],
	modeling_institution=p['modeling_institution'],
	convert_one_hour_to_three=p['convert_one_hour_to_three'],
	expected_time_step=p['expected_time_step'])

    print("completed.")
