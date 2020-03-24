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
        rapid_io_files_location=p['rapid_io_files_location'],
        lsm_data_location=p['lsm_data_location'],
        simulation_start_datetime=simulation_start_datetime,
        simulation_end_datetime=simulation_end_datetime,
        initial_flows_file=p['initial_flows_file'],
        run_rapid_simulation=p['run_rapid_simulation'],
        generate_initialization_file=p['generate_initialization_file'],
        file_datetime_pattern=p['file_datetime_pattern'],
        file_datetime_re_pattern=p['file_datetime_re_pattern'],
        generate_return_periods_file=p['generate_return_periods_file'],
        generate_seasonal_averages_file=p['generate_seasonal_averages_file'])

    print("completed.")
