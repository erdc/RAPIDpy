#! /usr/bin/env python

from RAPIDpy.gis import weight
import yaml
import sys

def get_params(filename):
    with open(filename, 'r') as f:
        p = yaml.load(f)
    return p

if __name__=='__main__':
    yaml_file = sys.argv[1]
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
