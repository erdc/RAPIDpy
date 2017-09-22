# -*- coding: utf-8 -*-
"""
    muskingum.py
    RAPIDpy

    Created by Alan D Snow, 2016.
    License: BSD 3-Clause
"""
from csv import reader as csv_reader
from csv import writer as csv_writer

import numpy as np
from past.builtins import xrange  # pylint: disable=redefined-builtin
from osgeo import gdal

# local
from ..helper_functions import csv_to_list, log, open_csv
from . import open_shapefile

# Enable GDAL/OGR exceptions
gdal.UseExceptions()


def CreateMuskingumKfacFile(in_drainage_line,
                            river_id,
                            length_id,
                            slope_id,
                            celerity,
                            formula_type,
                            in_connectivity_file,
                            out_kfac_file,
                            length_units="km",
                            slope_percentage=False,
                            file_geodatabase=None):
    r"""
    Creates the Kfac file for calibration.

    The improved methods using slope to generate values
    for Kfac were used here:

    Tavakoly, A. A., A. D. Snow, C. H. David, M. L. Follum, D. R. Maidment,
    and Z.-L. Yang, (2016) "Continental-Scale River Flow Modeling of the
    Mississippi River Basin Using High-Resolution NHDPlus Dataset",
    Journal of the American Water Resources Association (JAWRA) 1-22.
    DOI: 10.1111/1752-1688.12456

    Formula Type Options:

    1. :math:`Kfac_n = \frac{RiverLength_n}{Celerity_n}`
    2. :math:`Kfac_n = \eta*\frac{RiverLength_n}{\sqrt{RiverSlope_n}}`
    3. :math:`Kfac_n = \eta*\frac{RiverLength_n}{\sqrt{RiverSlope_n}}\left[0.05, 0.95\right]`


    Where:

    :math:`a = \frac{\sum_{n=1}^{r} \frac{RiverLength_n}{Celerity_n}}{r}`

    :math:`b = \frac{\sum_{n=1}^{r} \frac{RiverLength_n}{\sqrt{RiverSlope_n}}}{r}`

    :math:`\eta = \frac{a}{b}`

    r = Number of river segments.


    Parameters
    ----------
    in_drainage_line: str
        Path to the stream network (i.e. Drainage Line) shapefile.
    river_id: str
        The name of the field with the river ID
        (Ex. 'HydroID', 'COMID', or 'LINKNO').
    length_id: str
        The field name containging the length of the river segment
        (Ex. 'LENGTHKM' or 'Length').
    slope_id: str
        The field name containging the slope of the river segment
        (Ex. 'Avg_Slope' or 'Slope').
    celerity: float
        The flow wave celerity for the watershed in meters per second.
        1 km/hr or 1000.0/3600.0 m/s is a reasonable value if unknown.
    formula_type: int
        An integer representing the formula type to use when calculating kfac.
    in_connectivity_file: str
        The path to the RAPID connectivity file.
    out_kfac_file: str
        The path to the output kfac file.
    length_units: str, optional
        The units for the length_id field. Supported types are "m" for meters
        and "km" for kilometers.
    slope_percentage: bool, optional
        If True, it assumes the slope given is in percentage and will
        divide by 100. Default is False.
    file_geodatabase: str, optional
        Path to the file geodatabase. If you use this option, in_drainage_line
         is the name of the stream network feature class
         (WARNING: Not always stable with GDAL).


    Example::

        from RAPIDpy.gis.muskingum import CreateMuskingumKfacFile

        CreateMuskingumKfacFile(
            in_drainage_line='/path/to/drainageline.shp',
            river_id='LINKNO',
            length_id='Length',
            slope_id='Slope',
            celerity=1000.0/3600.0,
            formula_type=3,
            in_connectivity_file='/path/to/rapid_connect.csv',
            out_kfac_file='/path/to/kfac.csv',
            length_units="m",
        )
    """  # noqa
    ogr_drainage_line_shapefile_lyr, ogr_drainage_line_shapefile = \
        open_shapefile(in_drainage_line, file_geodatabase)

    number_of_features = ogr_drainage_line_shapefile_lyr.GetFeatureCount()
    river_id_list = np.zeros(number_of_features, dtype=np.int32)

    length_list = \
        np.zeros(number_of_features, dtype=np.float32)
    slope_list = np.zeros(number_of_features, dtype=np.float32)
    for feature_idx, drainage_line_feature in \
            enumerate(ogr_drainage_line_shapefile_lyr):
        river_id_list[feature_idx] = drainage_line_feature.GetField(river_id)
        length = drainage_line_feature.GetField(length_id)
        if length is not None:
            length_list[feature_idx] = length
        slope = drainage_line_feature.GetField(slope_id)
        if slope is not None:
            slope_list[feature_idx] = slope

    del ogr_drainage_line_shapefile

    if slope_percentage:
        slope_list /= 100.0

    if length_units == "m":
        length_list /= 1000.0
    elif length_units != "km":
        raise Exception("Invalid length units supplied. "
                        "Supported units are m and km.")

    connectivity_table = np.loadtxt(in_connectivity_file,
                                    delimiter=",",
                                    ndmin=2,
                                    dtype=int)

    length_slope_array = []
    kfac2_array = []
    if formula_type == 1:
        log("River Length/Celerity")
    elif formula_type == 2:
        log("Eta*River Length/Sqrt(River Slope)")
    elif formula_type == 3:
        log("Eta*River Length/Sqrt(River Slope) [0.05, 0.95]")
    else:
        raise Exception("Invalid formula type. Valid range: 1-3 ...")

    with open_csv(out_kfac_file, 'w') as kfacfile:
        kfac_writer = csv_writer(kfacfile)
        for row in connectivity_table:
            stream_id = int(float(row[0]))

            stream_id_index = river_id_list == stream_id
            # find the length
            stream_length = length_list[stream_id_index] * 1000.0

            if formula_type >= 2:
                # find the slope
                stream_slope = slope_list[stream_id_index]

                if stream_slope <= 0:
                    # if no slope, take average of upstream
                    # and downstream to get it
                    next_down_id = int(float(row[1]))
                    next_down_slope = 0
                    try:
                        next_down_index = \
                            np.where(river_id_list == next_down_id)[0][0]
                        next_down_slope = slope_list[next_down_index]
                    except IndexError:
                        pass

                    next_up_id = int(float(row[3]))
                    next_up_slope = 0
                    try:
                        next_up_index = \
                            np.where(river_id_list == next_up_id)[0][0]
                        next_up_slope = slope_list[next_up_index]
                    except IndexError:
                        pass

                    stream_slope = (next_down_slope + next_up_slope) / 2.0
                    if stream_slope <= 0:
                        # if still no slope, set to 0.001
                        stream_slope = 0.001

                length_slope_array.append(stream_length / stream_slope**0.5)
                kfac2_array.append(stream_length / celerity)
            else:
                kfac = stream_length / celerity
                kfac_writer.writerow(kfac)

        if formula_type >= 2:
            if formula_type == 3:
                log("Filtering Data by 5th and 95th Percentiles ...")
                length_slope_array = np.array(length_slope_array)
                percentile_5 = np.percentile(length_slope_array, 5)
                percentile_95 = np.percentile(length_slope_array, 95)

                length_slope_array[length_slope_array < percentile_5] = \
                    percentile_5
                length_slope_array[length_slope_array > percentile_95] = \
                    percentile_95

            eta = np.mean(kfac2_array) / np.mean(length_slope_array)
            log("Kfac2_Avg {0}".format(np.mean(kfac2_array)))
            log("Length_Slope Avg {0}".format(np.mean(length_slope_array)))
            log("Eta {0}".format(eta))
            log("Writing Data ...")
            for len_slope in length_slope_array:
                kfac_writer.writerow(eta*len_slope)


def CreateMuskingumKFile(lambda_k,
                         in_kfac_file,
                         out_k_file):
    """
    Creates muskingum k file from kfac file.

    Parameters
    ----------
    lambda_k: float
        The value for lambda given from RAPID after the calibration process.
        If no calibration has been performed, 0.35 is reasonable.
    in_kfac_file: str
        The path to the input kfac file.
    out_k_file: str
        The path to the output k file.


    Example::

        from RAPIDpy.gis.muskingum import CreateMuskingumKFile

        CreateMuskingumKFile(
            lambda_k=0.35,
            in_kfac_file='/path/to/kfac.csv',
            out_k_file='/path/to/k.csv')

    """
    kfac_table = csv_to_list(in_kfac_file)

    with open_csv(out_k_file, 'w') as kfile:
        k_writer = csv_writer(kfile)
        for row in kfac_table:
            k_writer.writerow([lambda_k * float(row[0])])


def CreateMuskingumXFileFromDranageLine(in_drainage_line,
                                        x_id,
                                        out_x_file,
                                        file_geodatabase=None):
    """
    Create muskingum X file from drainage line.

    Parameters
    ----------
    in_drainage_line: str
        Path to the stream network (i.e. Drainage Line) shapefile.
    x_id: str
        The name of the muksingum X field (i.e. 'Musk_x').
    out_x_file: str
        The path to the output x file.
    file_geodatabase: str, optional
        Path to the file geodatabase. If you use this option,
        in_drainage_line is the name of the stream network feature class
        (WARNING: Not always stable with GDAL).


    Example::

        from RAPIDpy.gis.muskingum import CreateMuskingumXFileFromDranageLine

        CreateMuskingumXFileFromDranageLine(
            in_drainage_line='/path/to/drainageline.shp',
            x_id='Musk_x',
            out_x_file='/path/to/x.csv')

    """
    ogr_drainage_line_shapefile_lyr, ogr_drainage_line_shapefile = \
        open_shapefile(in_drainage_line, file_geodatabase)

    with open_csv(out_x_file, 'w') as kfile:
        x_writer = csv_writer(kfile)
        for drainage_line_feature in ogr_drainage_line_shapefile_lyr:
            x_writer.writerow([drainage_line_feature.GetField(x_id)])

    del ogr_drainage_line_shapefile


def CreateConstMuskingumXFile(x_value,
                              in_connectivity_file,
                              out_x_file):
    """
    Create muskingum X file from value that is constant all the way through
    for each river segment.

    Parameters
    ----------
    x_value: float
        Value for the muskingum X parameter [0-0.5].
    in_connectivity_file: str
        The path to the RAPID connectivity file.
    out_x_file: str
        The path to the output x file.


    Example::

        from RAPIDpy.gis.muskingum import CreateConstMuskingumXFile

        CreateConstMuskingumXFile(
            x_value=0.3,
            in_connectivity_file='/path/to/rapid_connect.csv',
            out_x_file='/path/to/x.csv')

    """
    num_rivers = 0
    with open_csv(in_connectivity_file, "r") as csvfile:
        reader = csv_reader(csvfile)
        for _ in reader:
            num_rivers += 1

    with open_csv(out_x_file, 'w') as kfile:
        x_writer = csv_writer(kfile)
        for _ in xrange(num_rivers):
            x_writer.writerow([x_value])
