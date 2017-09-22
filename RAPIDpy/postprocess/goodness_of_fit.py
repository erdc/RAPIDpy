# -*- coding: utf-8 -*-
"""
    goodness_of_fit.py
    RAPIDpy

    Created by Alan D Snow, 2015.
    Based on RAPID_Toolbox for ArcMap
    License: BSD 3-Clause
"""
from __future__ import print_function
from csv import writer as csvwriter
import numpy as np

from ..dataset import RAPIDDataset


# ------------------------------------------------------------------------------
# statistic functions
# ------------------------------------------------------------------------------
# FUNCTIONS FROM http://pydoc.net/Python/ambhas/0.4.0/ambhas.errlib/
def filter_nan(s, o):
    """
        this functions removed the data  from simulated and observed data
        whereever the observed data contains nan

        this is used by all other functions, otherwise they will produce nan as
        output
        """
    data = np.array([s.flatten(), o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    return data[:, 0], data[:, 1]


def pc_bias(s, o):
    """
        Percent Bias
        input:
        s: simulated
        o: observed
        output:
        pc_bias: percent bias
        """
    # s,o = filter_nan(s,o)
    return 100.0*np.sum(s-o)/np.sum(o)


def apb(s, o):
    """
        Absolute Percent Bias
        input:
        s: simulated
        o: observed
        output:
        apb_bias: absolute percent bias
        """
    # s,o = filter_nan(s,o)
    return 100.0*np.sum(np.abs(s-o))/np.sum(o)


def rmse(s, o):
    """
        Root Mean Squared Error
        input:
        s: simulated
        o: observed
        output:
        rmses: root mean squared error
        """
    # s,o = filter_nan(s,o)
    return np.sqrt(np.mean((s-o)**2))


def mae(s, o):
    """
        Mean Absolute Error
        input:
        s: simulated
        o: observed
        output:
        maes: mean absolute error
        """
    # s,o = filter_nan(s,o)
    return np.mean(np.abs(s-o))


def bias(s, o):
    """
        Bias
        input:
        s: simulated
        o: observed
        output:
        bias: bias
        """
    # s,o = filter_nan(s,o)
    return np.mean(s-o)


def NS(s, o):
    """
        Nash Sutcliffe efficiency coefficient
        input:
        s: simulated
        o: observed
        output:
        ns: Nash Sutcliffe efficient coefficient
        """
    # s,o = filter_nan(s,o)
    return 1 - np.sum((s-o)**2)/np.sum((o-np.mean(o))**2)


def L(s, o, N=5):
    """
        Likelihood
        input:
        s: simulated
        o: observed
        output:
        L: likelihood
        """
    # s,o = filter_nan(s,o)
    return np.exp(-N*np.sum((s-o)**2)/np.sum((o-np.mean(o))**2))


def correlation(s, o):
    """
        correlation coefficient
        input:
        s: simulated
        o: observed
        output:
        correlation: correlation coefficient
        """
    # s,o = filter_nan(s,o)
    if s.size == 0:
        corr = np.NaN
    else:
        corr = np.corrcoef(o, s)[0, 1]
    return corr


def index_agreement(s, o):
    """
        index of agreement
        input:
        s: simulated
        o: observed
        output:
        ia: index of agreement
        """
    # s,o = filter_nan(s,o)
    ia = 1 - (np.sum((o-s)**2)) /\
             (np.sum((np.abs(s-np.mean(o))+np.abs(o-np.mean(o)))**2))
    return ia


def KGE(s, o):
    """
        Kling-Gupta Efficiency
        input:
        s: simulated
        o: observed
        output:
        kge: Kling-Gupta Efficiency
        cc: correlation
        alpha: ratio of the standard deviation
        beta: ratio of the mean
        """
    # s,o = filter_nan(s, o)
    cc = correlation(s, o)
    alpha = np.std(s)/np.std(o)
    beta = np.sum(s)/np.sum(o)
    kge = 1 - np.sqrt((cc-1)**2 + (alpha-1)**2 + (beta-1)**2)
    return kge, cc, alpha, beta

# END FUNCTIONS FROM http://pydoc.net/Python/ambhas/0.4.0/ambhas.errlib/


# ------------------------------------------------------------------------------
# Time Series comparison functions
# ------------------------------------------------------------------------------
def find_goodness_of_fit(rapid_qout_file, reach_id_file, observed_file,
                         out_analysis_file, daily=False):
    """
    Finds the goodness of fit comparing observed streamflow in a rapid Qout
    file with simulated flows in a csv file.

    Parameters
    ----------
    rapid_qout_file: str
        Path to the RAPID Qout file.
    reach_id_file: str
        ath to file with river reach ID's associate with the RAPID Qout file.
        It is in the format of the RAPID observed flows reach ID file.
    observed_file: str
        Path to input csv with with observed flows corresponding to the
        RAPID Qout. It is in the format of the RAPID observed flows file.
    out_analysis_file: str
        Path to the analysis output csv file.
    daily: bool, optional
        If True and the file is CF-Compliant, it will compare the
        *observed_file* with daily average flow from Qout. Default is False.


    Example with CF-Compliant RAPID Qout file:

    .. code:: python

        import os
        from RAPIDpy.postprocess import find_goodness_of_fit

        INPUT_DATA_PATH = '/path/to/data'
        reach_id_file = os.path.join(INPUT_DATA_PATH, 'obs_reach_id.csv')
        observed_file = os.path.join(INPUT_DATA_PATH, 'obs_flow.csv')

        cf_input_qout_file = os.path.join(COMPARE_DATA_PATH,
                                         'Qout_nasa_lis_3hr_20020830_CF.nc')
        cf_out_analysis_file = \
            os.path.join(OUTPUT_DATA_PATH,
                        'cf_goodness_of_fit_results-daily.csv')
        find_goodness_of_fit(cf_input_qout_file,
                             reach_id_file,
                             observed_file,
                             cf_out_analysis_file,
                             daily=True)

    """
    reach_id_list = np.loadtxt(reach_id_file,
                               delimiter=",", usecols=(0,),
                               ndmin=1, dtype=np.int32)

    data_nc = RAPIDDataset(rapid_qout_file)

    # analyze and write
    observed_table = np.loadtxt(observed_file,
                                ndmin=2, delimiter=",",
                                usecols=tuple(range(reach_id_list.size)))
    with open(out_analysis_file, 'w') as outcsv:
        writer = csvwriter(outcsv)
        writer.writerow(["reach_id",
                         "percent_bias",
                         "abs_percent_bias",
                         "rmse",
                         "mae",
                         "bias",
                         "NSE",
                         "likelihood",
                         "correlation_coeff",
                         "index_agreement",
                         "KGE"])

        for index, reach_id in enumerate(reach_id_list):
            observed_array = observed_table[:, index]
            simulated_array = data_nc.get_qout(reach_id, daily=daily)
            # make sure they are the same length
            simulated_array = simulated_array[:len(observed_array)]
            observed_array = observed_array[:len(simulated_array)]
            simulated_array, observed_array = \
                filter_nan(simulated_array, observed_array)
            writer.writerow([reach_id,
                             pc_bias(simulated_array, observed_array),
                             apb(simulated_array, observed_array),
                             rmse(simulated_array, observed_array),
                             mae(simulated_array, observed_array),
                             bias(simulated_array, observed_array),
                             NS(simulated_array, observed_array),
                             L(simulated_array, observed_array),
                             correlation(simulated_array, observed_array),
                             index_agreement(simulated_array, observed_array),
                             KGE(simulated_array, observed_array)[0]])


def find_goodness_of_fit_csv(observed_simulated_file, out_file=None):
    """
    Finds the goodness of fit comparing observed and simulated flows
    In the file, the first column is the observed flows and the
    second column is the simulated flows.

    Example::

        33.5, 77.2
        34.7, 73.0

    Parameters
    ----------
    observed_simulated_file: str
        Path to the csv file with the observed and simulated flows.
    out_file: str, optional
        Path to output file. If not provided, it will print to console.


    Example:

    .. code:: python

        from RAPIDpy.postprocess import find_goodness_of_fit_csv

        find_goodness_of_fit_csv('
            /united_kingdom-thames/flows_kingston_gage_noah.csv')

    """
    observed_simulated_table = np.loadtxt(observed_simulated_file,
                                          ndmin=2, delimiter=",",
                                          usecols=(0, 1))

    observed_array, simulated_array = \
        filter_nan(observed_simulated_table[:, 0],
                   observed_simulated_table[:, 1])

    # print error indices
    if out_file:
        print_file = open(out_file, 'w')
    else:
        print_file = None

    print("\n".join([
               "Percent Bias: {0:.4f}"
               .format(pc_bias(simulated_array, observed_array)),
               "Absolute Percent Bias: {0:.4f}"
               .format(apb(simulated_array, observed_array)),
               "Root Mean Squared Error: {0:.4f}"
               .format(rmse(simulated_array, observed_array)),
               "Mean Absolute Error: {0:.4f}"
               .format(mae(simulated_array, observed_array)),
               "Bias: {0}".format(bias(simulated_array, observed_array)),
               "Nash Sutcliffe efficiency coefficient: {0:.4f}"
               .format(NS(simulated_array, observed_array)),
               "Likelihood: {0:.4f}"
               .format(L(simulated_array, observed_array)),
               "correlation coefficient: {0:.4f}"
               .format(correlation(simulated_array, observed_array)),
               "index of agreement: {0:.4f}"
               .format(index_agreement(simulated_array, observed_array)),
               "Kling-Gupta Efficiency: {0:.4f}"
               .format(KGE(simulated_array, observed_array)[0]),
               ]),
          file=print_file)

    if print_file:
        print_file.close()
