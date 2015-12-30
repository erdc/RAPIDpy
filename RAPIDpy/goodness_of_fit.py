# -*- coding: utf-8 -*-
##
##  goodnessOfFit.py
##  RAPIDpy
##
##  Created by Alan D Snow 2015.
##  Copyright © 2015 Alan D Snow. All rights reserved.
##
__author__ = 'Alan D. Snow'
from csv import writer as csvwriter
import datetime
from netCDF4 import Dataset
import numpy as np
from pytz import utc

from helper_functions import csv_to_list

#------------------------------------------------------------------------------
#statistic functions
#------------------------------------------------------------------------------
## FUNCTIONS FROM http://pydoc.net/Python/ambhas/0.4.0/ambhas.errlib/
def filter_nan(s,o):
    """
        this functions removed the data  from simulated and observed data
        whereever the observed data contains nan
        
        this is used by all other functions, otherwise they will produce nan as
        output
        """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    return data[:,0],data[:,1]

def pc_bias(s,o):
    """
        Percent Bias
        input:
        s: simulated
        o: observed
        output:
        pc_bias: percent bias
        """
    #s,o = filter_nan(s,o)
    return 100.0*sum(s-o)/sum(o)

def apb(s,o):
    """
        Absolute Percent Bias
        input:
        s: simulated
        o: observed
        output:
        apb_bias: absolute percent bias
        """
    #s,o = filter_nan(s,o)
    return 100.0*sum(abs(s-o))/sum(o)

def rmse(s,o):
    """
        Root Mean Squared Error
        input:
        s: simulated
        o: observed
        output:
        rmses: root mean squared error
        """
    #s,o = filter_nan(s,o)
    return np.sqrt(np.mean((s-o)**2))

def mae(s,o):
    """
        Mean Absolute Error
        input:
        s: simulated
        o: observed
        output:
        maes: mean absolute error
        """
    #s,o = filter_nan(s,o)
    return np.mean(abs(s-o))

def bias(s,o):
    """
        Bias
        input:
        s: simulated
        o: observed
        output:
        bias: bias
        """
    #s,o = filter_nan(s,o)
    return np.mean(s-o)

def NS(s,o):
    """
        Nash Sutcliffe efficiency coefficient
        input:
        s: simulated
        o: observed
        output:
        ns: Nash Sutcliffe efficient coefficient
        """
    #s,o = filter_nan(s,o)
    return 1 - sum((s-o)**2)/sum((o-np.mean(o))**2)

def L(s,o, N=5):
    """
        Likelihood
        input:
        s: simulated
        o: observed
        output:
        L: likelihood
        """
    #s,o = filter_nan(s,o)
    return np.exp(-N*sum((s-o)**2)/sum((o-np.mean(o))**2))

def correlation(s,o):
    """
        correlation coefficient
        input:
        s: simulated
        o: observed
        output:
        correlation: correlation coefficient
        """
    #s,o = filter_nan(s,o)
    if s.size == 0:
        corr = np.NaN
    else:
        corr = np.corrcoef(o, s)[0,1]
    
    return corr


def index_agreement(s,o):
    """
        index of agreement
        input:
        s: simulated
        o: observed
        output:
        ia: index of agreement
        """
    #s,o = filter_nan(s,o)
    ia = 1 -(np.sum((o-s)**2))/(np.sum((np.abs(s-np.mean(o))+np.abs(o-np.mean(o)))**2))
    return ia


def KGE(s,o):
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
    #s,o = filter_nan(s,o)
    cc = correlation(s,o)
    alpha = np.std(s)/np.std(o)
    beta = np.sum(s)/np.sum(o)
    kge = 1- np.sqrt( (cc-1)**2 + (alpha-1)**2 + (beta-1)**2 )
    return kge, cc, alpha, beta

def assimilation_eff(assimilated, simulated, observed):
    """
        Assimilation efficiency (Aubert et al., 2003)
        Input:
        assimilated: assimilated flow
        simulated: simulated flow
        observed: observed flow
        Output:
        Eff
        """
    s,o = filter_nan(simulated, observed)
    a,o = filter_nan(assimilated, observed)
    
    Eff = 100*(1 - np.sum((a-o)**2)/np.sum((s-o)**2))
    return Eff
## END FUNCTIONS FROM http://pydoc.net/Python/ambhas/0.4.0/ambhas.errlib/

#------------------------------------------------------------------------------
#Time Series comparison functions
#------------------------------------------------------------------------------
def find_goodness_of_fit(reach_id_file, rapid_qout_file, observed_file,
                         out_analysis_file, daily=False, steps_per_group=1):
    """
        Finds the goodness of fit comparing observed and simulated flows
    """
    reach_id_list = np.array([row[0] for row in csv_to_list(reach_id_file)])
   
    nc_file = Dataset(rapid_qout_file)
    nc_reach_id_list = nc_file.variables['COMID'][:]
    
    #analyze and write
    def get_daily_time_series(data_nc, reach_index, daily, steps_per_group):
        """
            Gets the daily time series from RAPID output
        """
        if 'time' in data_nc.variables.keys():
            if daily:
                current_day = datetime.datetime.fromtimestamp(data_nc.variables['time'][0], tz=utc)
                flow = 0.0
                num_days = 0
                qout_arr = data_nc.variables['Qout'][reach_index, :]
                daily_qout = []
                for idx, t in enumerate(data_nc.variables['time'][:]):
                    var_time = datetime.datetime.fromtimestamp(t, tz=utc)
                    if current_day.day == var_time.day:
                        flow += qout_arr[idx]
                        num_days += 1
                    else:
                        if num_days > 0:
                            #write last average
                            daily_qout.append(flow/num_days)
                        
                        #start new average
                        current_day = var_time
                        num_days = 1
                        flow = qout_arr[idx]
                
                return np.array(daily_qout, np.float32)
            return nc_file.variables['Qout'][reach_index,:]
        elif steps_per_group > 1:
            flow_data = data_nc.variables['Qout'][:, reach_index]
            daily_qout = []
            for step_index in xrange(0, len(flow_data), steps_per_group):
                flows_slice = flow_data[step_index:step_index + steps_per_group]
                daily_qout.append(np.mean(flows_slice))
            
            return np.array(daily_qout, np.float32)
        
        return nc_file.variables['Qout'][:,reach_index]


    observed_table = np.array(csv_to_list(observed_file), np.float32)
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
            reach_index = np.where(nc_reach_id_list == int(reach_id))[0][0]
            observed_array = observed_table[:, index]
            simulated_array = get_daily_time_series(nc_file, reach_index, daily, steps_per_group)
            #make sure they are the same length
            simulated_array = simulated_array[:len(observed_array)]
            observed_array = observed_array[:len(simulated_array)]
            simulated_array,observed_array = filter_nan(simulated_array,observed_array)
            writer.writerow([reach_id,
                             pc_bias(simulated_array,observed_array),
                             apb(simulated_array,observed_array),
                             rmse(simulated_array,observed_array),
                             mae(simulated_array,observed_array),
                             bias(simulated_array,observed_array),
                             NS(simulated_array,observed_array),
                             L(simulated_array,observed_array),
                             correlation(simulated_array,observed_array),
                             index_agreement(simulated_array,observed_array),
                             KGE(simulated_array,observed_array)[0]])

def find_goodness_of_fit_csv(observed_simulated_file):
    """
        Finds the goodness of fit comparing observed and simulated flows
        In file:
        observed_flow, simulated flow
        
        Ex.
        
        33.5, 77.2
        34.7, 73.0
    """
    flow_table = csv_to_list(observed_simulated_file)
    
    observed_array = []
    simulated_array = []
    for row in flow_table:
        observed_array.append(float(row[0]))
        simulated_array.append(float(row[1]))
    
    simulated_array,observed_array = filter_nan(np.array(simulated_array),
                                                np.array(observed_array))

    # print error indices
    print "Percent Bias:", pc_bias(simulated_array,observed_array)
    print "Absolute Percent Bias:", apb(simulated_array,observed_array)
    print "Root Mean Squared Error:", rmse(simulated_array,observed_array)
    print "Mean Absolute Error", mae(simulated_array,observed_array)
    print "Bias", bias(simulated_array,observed_array)
    print "Nash Sutcliffe efficiency coefficient", NS(simulated_array,observed_array)
    print "Likelihood", L(simulated_array,observed_array)
    print "correlation coefficient", correlation(simulated_array,observed_array)
    print "index of agreement", index_agreement(simulated_array,observed_array)
    print "Kling-Gupta Efficiency", KGE(simulated_array,observed_array)[0]



#------------------------------------------------------------------------------
#main
#------------------------------------------------------------------------------
"""
if __name__ == "__main__":
    find_goodness_of_fit_csv('/Users/rdchlads/autorapid/rapid-io/output/united_kingdom-thames/flows_kingston_gage.csv')
    find_goodness_of_fit(reach_id_file='/Users/rdchlads/autorapid/rapid-io/output/united_kingdom-thames/obs_reach_id.csv',
                         rapid_qout_file='/Users/rdchlads/autorapid/rapid-io/output/united_kingdom-thames/Qout_met_office_joules_3hr_20010101to20090101_ukv.nc',
                         observed_file='/Users/rdchlads/autorapid/rapid-io/output/united_kingdom-thames/flows_kingston_gage.csv',
                         out_analysis_file='/Users/rdchlads/autorapid/rapid-io/output/united_kingdom-thames/flows_kingston_analysis_ukv.csv',
                         daily=True)
"""