import glob
import subproccess
import os
from RAPIDpy.gis import weight
import yaml
import sys
import copy
from RAPIDpy.inflow import run_lsm_rapid_process
from datetime import datetime

os.chdir('/p/home/nolsen2/RAPIDpy/workflow')


basins = glob.glob('/p/work/nolsen2/rapid-io/input/*')


for basin in basins:

    basinName = basin.split('/')[-1]
    
    # make output directory
    os.mkdir('/p/work/nolsen2/rapid-io/output/'+basinName)


    ##### WEIGHT TABLE
	#build weight table yaml
    
    weighttableyamlFilename = '/p/home/nolsen2/LIS/weight_lis.yml'
    
    in_catchment_shapefile = '/p/work/nolsen2/rapid-io/input/'+basinName+'/catchmentV.shp'
    in_connectivity_file = '/p/work/nolsen2/rapid-io/input/'+basinName+'/rapid_connect.csv'
    out_weight_table = '/p/home/nolsen2/LIS/weight_lis.csv'

    with open(weighttableyamlFilename, 'a') as the_file:
        the_file.write('lsm: \'LDAS\'')
        the_file.write('in_nc: \'/p/app/unsupported/LIS_CPRA/jules/SURFACEMODEL/allfiles/LIS_HIST_201511020000.d01.nc\'')
        the_file.write('in_nc_lon_var: \'lon\'')
        the_file.write('in_nc_lat_var: \'lat\'')
        the_file.write('in_catchment_shapefile: \''+in_catchment_shapefile+'\'')
        the_file.write('river_id: \'GHID_2\'')
        the_file.write('in_connectivity_file: \''+in_connectivity_file+'\'')
        the_file.write('out_weight_table: \''+out_weight_table+'\'')
        the_file.write('area_id: ')
        the_file.write('file_geodatabase: ')

    # # send command to rapidpy
    # bashCommand = "python generate_rapid_weight_table.py /p/home/nolsen2/LIS/weight_lis.yml"
    # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()
    
    def get_params(filename):
        with open(filename, 'r') as f:
            p = yaml.load(f)
        return p

    yaml_file = copy.copy(weighttableyamlFilename)
    p = get_params(yaml_file)
    lsm = p['lsm'].upper() 
    if lsm in ['LDAS', 'JULES', 'WSIM']:
        weight.CreateWeightTableLDAS(
            in_ldas_nc=p['in_nc'],
            in_nc_lon_var=p['in_nc_lon_var'],
            in_nc_lat_var=p['in_nc_lat_var'],
            in_catchment_shapefile=p['in_catchment_shapefile'],
            river_id=p['river_id'],
            in_connectivity_file=p['in_connectivity_file'],
            out_weight_table=p['out_weight_table'],
            area_id=p['area_id'],
            file_geodatabase=p['file_geodatabase'])
    elif lsm == 'ECMWF':
        weight.CreateWeightTableECMWF(
            in_ecmwf_nc=p['in_nc'],
            in_catchment_shapefile=p['in_catchment_shapefile'],
            river_id=p['river_id'],
            in_connectivity_file=p['in_connectivity_file'],
            out_weight_table=p['out_weight_table'],
            area_id=p['area_id'],
            file_geodatabase=p['file_geodatabase'])

    print("completed.")    
    
    
    
    
    
    ##### INFLOW
    
    #build inflow yaml
    inflowyamlFilename = '/p/home/nolsen2/LIS/inflow_first.yml'
    
    rapid_executable_location='/p/home/nolsen2/rapid/src/rapid'
    rapid_input_location='/p/work/nolsen2/rapid-io/input/'+basinName+'/'
    rapid_output_location='/p/work/nolsen2/rapid-io/output/'+basinName+'/'
    lsm_data_location='/p/app/unsupported/LIS_CPRA/jules/SURFACEMODEL/allfiles/'
    simulation_start_timestamp='20080101'
    simulation_end_timestamp='20200908'

    with open(inflowyamlFilename, 'a') as the_file:
        the_file.write('rapid_executable_location: \''+rapid_executable_location+'\'')
        the_file.write('rapid_io_files_location: ')
        the_file.write('rapid_input_location: \''+rapid_input_location+'\'')
        the_file.write('rapid_output_location: \''+rapid_output_location+'\'')
        the_file.write('lsm_data_location: \''+lsm_data_location+'\'')
        the_file.write('simulation_start_timestamp: \''+simulation_start_timestamp+'\'')
        the_file.write('simulation_end_timestamp: \''+simulation_end_timestamp+'\'')
        the_file.write('run_rapid_simulation: false')
        the_file.write('generate_initialization_file: false')
        the_file.write('file_datetime_pattern: \'%Y%m%d%H%M\'')
        the_file.write('file_datetime_re_pattern: \'\d{12}\'')
        the_file.write('generate_return_periods_file: false')
        the_file.write('generate_seasonal_averages_file: false')
        the_file.write('generate_seasonal_initialization_file: false')
        the_file.write('initial_flows_file:')
        the_file.write('generate_rapid_namelist_file: true')
        the_file.write('run_rapid_simulation: false')
        the_file.write('generate_return_periods_file: false')
        the_file.write('return_period_method: \'weibul\'')
        the_file.write('use_all_processors: true')
        the_file.write('num_processors: 1')
        the_file.write('mpiexec_command: \'mpiexec\'')
        the_file.write('cygwin_bin_location: \'\'')
        the_file.write('modeling_institution: \'US Army Engineer Research and Development Center\'')
        the_file.write('convert_one_hour_to_three: false')
        the_file.write('expected_time_step:')

    # # send command to rapidpy
    # bashCommand = "python generate_rapid_inflow_files.py /p/home/nolsen2/LIS/inflow_first.yml"
    # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()

    def get_params(filename):
        with open(filename, 'r') as f:
            p = yaml.load(f)
        return p

    yaml_file = copy.copy(inflowyamlFilename)
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




    
    
    # move all files
    
    b = glob.glob('/p/home/nolsen2/LIS/*')
    for bb in b:
        shutil.move(bb, '/p/work/nolsen2/rapid-io/input/'+basinName+'/'+bb.split('/')[-1])
    
    b = glob.glob('/p/work/nolsen2/rapid-io/output/'+basinName+'/*')
    for bb in b:
        shutil.move(bb, '/p/work/nolsen2/rapid-io/input/'+basinName+'/'+bb.split('/')[-1])
    
    
    
    ##### RAPID ??????
    
    
    
    
    
    
    
    
    
    
    