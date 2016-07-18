# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:03:39 2016
ncdump() was created by this dude:
http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html#code
JD_converter was created by this dude: (and modified by me)
https://asimpleweblog.wordpress.com/2010/06/20/julian-date-calculator/

@author: login
"""

import scipy.stats as st
import re
import mmap
import os
import numpy as np
import pandas as pd
from pandas.tseries.offsets import DateOffset
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import JD_converter as jd
import subprocess as sp
import copy
import time
from functools import wraps
import shutil
 
 
def fn_timer(function):
    """
    This decorator can be used to print the time it takes a fxn to execute
    """
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("Total time running %s: %s seconds" %
               (function.func_name, str(t1-t0))
               )
        return result
    return function_timer
    
def import_config(dl_default=False, verbose=True):
    # the glm contains the following blocks:
    # &glm_setup: 13 general simulation info and mixing parameters
    if dl_default:
        default_config = {'others':{}}
        default_order = {}
        block_order = []
        with open('../test_files/sample_glm.nml', 'r') as f:
            for line in f:
                key = line[0]
                line = line.strip()
                if key == '&':
                    new_block = line[1:]
                    orderedList = []
                    default_config[new_block] = {}
                elif key == '/':
                    default_order[new_block] = orderedList
                    block_order.append(new_block)
                    if verbose == True:
                        print "Section End"
                elif key != '!':
                    line_split = line.split('=')
                    line_ = [l.strip() for l in line_split]
                    orderedList.append(line_[0])
                    if verbose == True:
                        print line_
                    if line_[0][0:3] != 'aed':
                        default_config[new_block][line_[0]] = line_
                    
                else:
                    default_config['others'][line_split[0]] = line_split

    else:
        default_config, default_order, block_order = None, None, None
        
    glm_setup = {'max_layers' :200, 
                 'min_layer_vol' :0.025, 
                 'min_layer_thick' :0.25, 
                 'max_layer_thick' :0.500,  
                 'Kw' : 0.5, 
                 'coef_mix_conv' : 0.125, 
                 'coef_wind_stir' : 0.23, 
                 'coef_mix_shear' : 0.20, 
                 'coef_mix_turb' : 0.51, 
                 'coef_mix_KH' : 0.30, 
                 'coef_mix_hyp' : 0.5,
                 'deep_mixing':'.true.'}
                 
    wq_setup = {}
    
    # &morphometry: 8 vars
    morphometry = {'lake_name': "'UpperMysticLake'", 
                   "latitude": 42.4317,  
                   "longitude": -71.1483,
                   "bsn_len": 1073.637,
                   "bsn_wid": 632.60, 
                   "bsn_vals": None,
                   "H": [-23.384, -20.336, -17.288, -14.24, -11.192, -8.144, 
                         -5.096, -2.048, 1.00],
                   "A":[ 77373.80, 148475.73, 202472.97, 257818.95, 338552.69, 
                        397077.50, 460778.04, 524802.66, 560051.22]}
                        
    morphometry['bsn_vals'] = len(morphometry['H'])
    assert len(morphometry['H']) == len(morphometry['A'])
    
    # &time block 5 vars if 'timefmt' is 2
    time = {"timefmt" : 3, 
            "start" : "'2012-01-01 00:00:00'", 
            "stop" : "'2014-01-01 00:00:00'", 
            "dt" : 3600.0,
            "timezone" : 5.0 }
            
    start_t = pd.Timestamp(time["start"])
    end_t = pd.Timestamp(time["stop"])
    time['num_days'] = (end_t-start_t).days
    
    # &output 14 vars 
    output = {"out_dir" : "", 
              "out_fn" : "'output_'",
              "nsave" : 12, 
              "csv_lake_fname" : "'lake'",
              "csv_point_nlevs" : 0,
              "csv_point_fname" : "'WQ_'",
              "csv_point_at" : [1, 5, 9, 13, 17, 21],
              "csv_point_nvars" : 2,
              "csv_point_vars" : [ 'temp', 'salt'],
              "csv_outlet_allinone" : ".false.",
              "csv_outlet_fname" : "'outlet_'",
              "csv_outlet_nvars" : 3,
              "csv_outlet_vars" : ['flow', 'temp', 'salt'],
              "csv_ovrflw_fname" : "\"overflow\"" }
              
    output["csv_point_nlevs"] = len(output["csv_point_at"])
    
    #&init_profiles            
    init_profiles = {"lake_depth": 24.384, 
                     "num_depths": 0, 
                     "the_depths": [1, 5, 9, 13, 17, 21], 
                     "the_temps": [4.0, 4.0, 4.0, 4.0, 4.0, 4.0], 
                     "the_sals": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                     }
    init_profiles['num_depths'] = len(init_profiles["the_depths"])
    assert len(init_profiles["the_depths"]) == len(init_profiles["the_temps"])
    assert len(init_profiles["the_depths"]) == len(init_profiles["the_sals"])
    
    if wq_setup.keys() != []:        
        init_profiles["num_wq_vars"] = 0
        init_profiles["wq_names"] = []
        init_profiles["wq_init_vals"] = "" 
                     
    #&meteorology                 
    meteorology = {"met_sw" : ".true.",
                   "lw_type" : "'LW_IN'", 
                   "rain_sw" : ".false.",
                   "atm_stab" : ".false.",
                   "catchrain" : ".false.",
                   "rad_mode": 1,
                   "albedo_mode": 3, 
                   "cloud_mode": 3, 
                   "subdaily": '.false.', 
                   "meteo_fl": "'met_daily.csv'",
                   "wind_factor": 1.0,
                   "sw_factor":  1.0,
                   "lw_factor": 1.0, 
                   "at_factor": 1.0,
                   "rh_factor": 1.0,
                   "rain_factor": 1.0,
                   "ce":  0.0013, 
                   "ch": 0.0013, 
                   "cd": 0.0013, 
                   "rain_threshold": 0.01, 
                   "runoff_coef":  0.3}
    #&bird_model            
    bird = {"AP" : 973, 
            "Oz" : 0.279, 
            "WatVap" : 1.1, 
            "AOD500" : 0.033, 
            "AOD380" : 0.038, 
            "Albedo" : 0.2}
    
    #&inflow                
    inflow = {"num_inflows": 1,
              "names_of_strms": "'Aberjona'", 
              "subm_flag": ".false.", 
              "strm_hf_angle": 65.0, 
              "strmbd_slope": 2.0, 
              "strmbd_drag": 0.0160, 
              "inflow_factor": 1.0, 
              "inflow_fl": "'inflow.csv'", 
              "inflow_varnum": 3, 
              "inflow_vars": ['FLOW','TEMP','SALT'],
              "coef_inf_entrain": 0. }
              
    #&outflow
    outflow = {"num_outlet": 1,
               "flt_off_sw": ".false.",
               "outl_elvs": 0.00,
               "bsn_len_outl": 299, 
               "bsn_wid_outl" : 199,
               "outflow_fl" : "'outflow.csv'",
               "outflow_factor": 0.8,
               "seepage" : ".true.",
               "seepage_rate" : 0.0 }
                
    glm_config = { "glm_setup" : glm_setup, 
                   "wq_setup" : wq_setup, 
                   "morphometry" : morphometry, 
                   "time" : time,
                   "output" : output,
                   "init_profiles" : init_profiles,
                   "meteorology" : meteorology,
                   "bird_model" : bird, 
                   "outflow" : outflow,
                   "inflow" : inflow }

    return glm_config, default_config, default_order, block_order

def log_likelihood(err, t_obs):
    std=100    #EDIT THIS VALUE
    log_lik=0
    #mutiply prob. of observation at each time t
    for t in np.arange(t_obs):
        lik_this = st.norm.pdf(err[t].mean(), 0, std)
        log_lik += np.log(lik_this)
    return np.exp(log_lik)
    
def error_metrics(Obs, Sim):
    """ returns MSE, rMSE, NSE, Pearson r, and others if given two 
    numpy ndarrays"""
    t_obs, _ = Sim.shape
    Sim_f, Obs_f = Sim.flatten(), Obs.flatten()
    Err = Sim - Obs
    MSE = np.mean(Err**2)
    rMSE = np.sqrt(MSE)
    NSE = 1. - (Err.var()/np.var(Obs))
    r = st.pearsonr(Sim_f, Obs_f)        
    a = np.std(Sim) / np.std(Obs)
    b = np.mean(Sim) / np.mean(Obs)
    # Kling-Gupta Efficiency
    KGE = 1- np.sqrt((r[0]-1)**2 + (a-1)**2 + (b-1)**2)
    return {'MSE':MSE, 'rMSE':rMSE, 'NSE':NSE, 'r':r, 'a':a, 'b':b, 'KGE':KGE,
            'LOG': log_likelihood(Err, t_obs) }

def subsectbydate_2(df1, df2):
    lowEnd = np.max([df1.index[0],df2.index[0]])
    highEnd = np.min([df1.index[-1],df2.index[-1]])
    c1, c2 = df1.index >= lowEnd, df2.index >= lowEnd
    c3, c4 = df1.index <= highEnd, df2.index <= highEnd
    return df1[c1 & c3].copy(), df2[c2 & c4].copy()

def subsectbydate_1(largedf, smalldf):
    crit1 = largedf.index >= smalldf.index[0]
    crit2 = largedf.index <= smalldf.index[-1]
    bothCrit = crit1 & crit2
    return largedf[bothCrit].copy()

def addtodict(vals, vars, dict):
    for var, val in zip(vars, vals):
        dict[var] = val
    
def make_dir(s, verbose=True):
    if os.path.exists(s):
        if verbose==True:
            print "%s exists\n" % os.path.basename(s)
    else:
        os.mkdir( s, 0760)

def printNaNWarning(df, label, fill_method=None):
    print "\nNaN Values read into {0} ({1})".format(label, str(type(df)))
    print df.isnull().sum()
    print "Total measurements per variable: {}".format(df.shape[0])
    return None

def plotTemporalDist(df, fignum, clear=True, bins=None):
    days = df.index.dayofyear
    if bins == None:
        buckets = len(np.unique(days))
    else:
        buckets = bins
        
    plt.figure(fignum)
    
    if clear == True:
        plt.clf()
        
    n, bins, patches = plt.hist(days, buckets, facecolor='green', alpha=0.75)
    
    return n, bins, patches

def pandasfillwrap(seriesORdf, fill_method):
    if fill_method=="ZERO":
        print "\tSetting them equal to zero"
        seriesORdf.fillna(value=0, inplace=True)
    elif fill_method in ['backfill', 'bfill', 'pad', 'ffill']:
        seriesORdf.fillna(method=fill_method, inplace=True)
        print "\t using Pandas {} method".format(fill_method)
    else:
        print "\t FILL METHOD NOT RECOGNIZED"
        print "\tLeaving them alone"
    return None

def z_score(df):
    z_df = pd.DataFrame(index=df.index, columns=df.columns)    
    for col in df.columns:
        z_df[col] = (df[col] - np.mean(df[col])) / np.std(df[col], ddof=1)
    return z_df
    
def makeAggregations(df, indices, fxn):
    aggs = {}    
    for idx in indices:
        new_df = df.copy().groupby(df[idx]).agg(fxn)
        aggs[idx] = new_df.drop(indices, axis=1)
    
    return aggs
    
def show_agg_resolution(aggs, col_name):
    day_frame = pd.DataFrame(index=aggs['day_i'].index, columns=aggs.keys())
    for col in day_frame.columns:
        ser = aggs[col][col_name]
        if col == 'day_i':
            new_idx = ser.index
        elif col == 'week_i':
            new_idx = ser.index*7+2
        elif col == 'month_i':
            new_idx = ser.index*30+6
        elif col == 'season_i':
            new_idx = ser.index*(366/4)+2
        elif col == 'year_i':
            new_idx = pd.Index([366], dtype='Int64')
            
        if col == 'year_i':
            mod_ser = pd.Series(data=ser.values.mean(), index=new_idx)
            day_frame[col] = mod_ser.reindex(day_frame.index)
        else:
            mod_ser = pd.Series(data=ser.values, index=new_idx)
            day_frame[col] = mod_ser.reindex(day_frame.index)
        
    return day_frame


# TODO: def plotAggregations():
# TODO: redo this to apply to aggregations 

def printAutocorr(aggs, threshold=None):
    print "Autocorrelation peaks"    
#    data_cols = {}
#    ac_df = pd.DataFrame(index=df.index, columns=df.columns)
#    
#    for col in df.columns:
#        if df[col].dtype != '<M8[ns]':
#            print "\t{}:".format(col)
#            autocorrs = np.array([df[col].autocorr(i) for i in range(1, len(df[col]))])
#            data_cols[col] = autocorrs
#            peaks = peakutils.indexes(autocorrs, thres=0.5, min_dist=100)
#            #peakvals = {autocorrs[j]:j for j in peaks}
#            print "\t\t{} peaks detected".format(len(peaks))
#            print "\t\tTallest peak ({0}) @ lag {1}".format(autocorrs.max(),
#                                                            autocorrs.argmax())
#            
#    for key in data_cols.keys():
#            ac_df[key] = data_cols[key]
    for key in aggs.keys():
        df = aggs[key].dropna()
        print key
        for col in df.columns:
            if col[-1] != "i":
                autocorrs = np.array([df[col].autocorr(i) for i in range(1, len(df[col])/2)])
                print "\t{0} max autocorr @ {1} = {2:.2f}".format(col, np.argmax(autocorrs)+1, autocorrs.mean())
            
    return None


# TODO: def plotfrequencydomain():

def TimeIdx(df, dates=None, insertAt=0):
    """this inserts day of the year, month of the year, and season of the year
    columns into a dataframe for use by the groupby function"""

    if dates == None:
        dates = df.index
    else:
        pass
    
    dates = pd.date_range(dates[0], periods=len(dates))
    month_i = dates.month
    week_i = dates.week
    day_i = dates.dayofyear
    year_i = dates.year
    fall, winter, spring, summer = [10,11,12], [1,2,3], [4,5,6], [7,8,9]
    season_i = np.zeros(len(month_i))
    for idx in range(len(month_i)):
        if month_i[idx] in fall:
            season_i[idx] = 1
        elif month_i[idx] in winter:
            season_i[idx] = 2
        elif month_i[idx] in spring:
            season_i[idx] = 3
        elif month_i[idx] in summer:
            season_i[idx] = 4
        else:
            print "unexpected error"

    i_names = ["year_i", "month_i", "season_i", 'week_i', "day_i"]
    new_cols = [year_i, month_i, season_i, week_i, day_i]
    for ind, name, col in zip(range(insertAt,len(new_cols)), i_names, new_cols):
        df.insert(ind, name, col)
    return df, i_names

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print "NetCDF Global Attributes:"
        for nc_attr in nc_attrs:
            print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print "NetCDF dimension information:"
        for dim in nc_dims:
            print "\tName:", dim 
            print "\t\tsize:", len(nc_fid.dimensions[dim])
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print "NetCDF variable information:"
        for var in nc_vars:
           if var not in nc_dims:
                print '\tName:', var
                print "\t\tdimensions:", nc_fid.variables[var].dimensions
                print "\t\tsize:", nc_fid.variables[var].size
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars

            
class Lake(object):
    
    def __init__(self, name, dir_path):
        self.stdout, self.stderr, self.ran = None, None, False
        self.variant_var, self.variant_val = "baseline", "baseline"
        
        self.name = name
        self.dir_path = dir_path
        make_dir(dir_path, verbose=False)
        self.glm_path = os.path.join(self.dir_path,"glm2.nml")
        
        config_pack = import_config(dl_default=True, verbose=False)
        self.glm_config, self.default_config = config_pack[0], config_pack[1]
        self.default_order, self.block_order = config_pack[2], config_pack[3]
        self.baseline_configured=False
       
    def write_glm_config(self, verbose=True):
        self.glm_config['output']['out_dir'] = "'"+str(self.dir_path)+"'"
        self.glm_config['output']['out_fn'] = self.glm_config['output']['out_fn'][0:-1] + self.name+"'"
        self.glm_config['glm_setup']['sim_name'] = "'"+str(self.name)+"'"
        
        if verbose:
            print "Formatting Check"
            if self.default_config != None:
                for key in self.glm_config.keys():
                    if key in self.default_config.keys():
                        default = self.default_config[key].keys()
                        glm = self.glm_config[key].keys()
                        if verbose == True:
                            print key
                            print "\tonly in default: ", list(set(default)-set(glm))
                            print "\tonly in glm: ", list(set(glm)-set(default))
                    else:
                        print key, " is not in default config"


        # start writing config file
        if os.path.exists(self.glm_path):
            if verbose==True:
                print "\n\nOverwriting Existing `glm2.nml` Configuration File\n"
        glm_handle = open(self.glm_path, 'w+')
        non_empty = []
        for key in self.block_order:
            if self.glm_config[key].keys() != []:
                non_empty.append(key)
        self.non_empty = non_empty
        for block in self.non_empty:
            glm_handle.write("&"+block+"\n")
            
            for param in self.default_order[block]:
                if verbose == True:
                    print block, param
                value = self.glm_config[block][param]
                if type(value) == list:
                    glm_handle.write("   {0} = ".format(param))
                    for i, v in enumerate(value):
                        if (i+1) == len(value):
                            glm_handle.write("%r\n" % v)
                        else:
                            glm_handle.write("%r, " % v)
                else:
                    glm_handle.write("   {0} = {1}\n".format(param, value))
            glm_handle.write("/\n")
        glm_handle.close()
        
        self.baseline_configured = True    
        
    def read_GEODISC_cloud_data(self, fname=None, data_dir=None):
        if fname == None:
            self.geodisc_f = 'cloud_data.csv'
        else:
            self.geodisc_f = fname
        if data_dir == None:
            self.geodisc_d = os.path.join(os.getcwd(), "GEODESC")
        else:
            self.geodisc_d = data_dir

        self.geodisc_p = os.path.join(self.geodisc_d, self.geodisc_f)        
        
        try:
            self.geodisc_clouds = pd.read_csv(self.geodisc_p, index_col = 0,
                                              parse_dates = ['start', 'end'],
                                              infer_datetime_format = True)
        except IOError:
            print "no file detected"
            raise IOError
    
    def temporal_clipping(self, data_pack):
        
        self.firstDay = np.max([i.index[0] for i in data_pack])
        self.lastDay = np.min([i.index[-1] for i in data_pack])
        frontCrit = [i.index >= self.firstDay for i in data_pack]
        backCrit = [i.index <= self.lastDay for i in data_pack]
        self.aligned_columns = {}
        
        for i, j, k in zip(frontCrit, backCrit, data_pack):
            self.aligned_columns[k.name] = k[i & j].copy()

        self.date_range = pd.date_range(start=self.firstDay, end=self.lastDay)
    
    def create_metcsv(self, metPack):
        """
        This creates the meteorological variable CSV. It pulls the 
        the expected file name from the `expectations` object and writes a 
        csv to that location. The format is as follows:
        #time,ShortWave,LongWave, AirTemp,RelHum, WindSpeed,Rain,Snow
        #1997-01-01,128.1119583,273.7329402,15.48908333,76.07847214, 0.950309798, 0,0       
        """
        self.config_dir = os.path.dirname(self.glm_path)
        met_fn = self.glm_config['meteorology']['meteo_fl'][1:-1]
        self.metcsv_path = os.path.join(self.config_dir, met_fn)
        self.met_df = pd.DataFrame(index= self.date_range, columns= metPack)
        for i in metPack:
            self.met_df[i] = self.aligned_columns[i]
            
        print "\tWriting met.csv"
        self.met_df.to_csv(path_or_buf=self.metcsv_path, index_label='time')
        
    def create_flowcsvs(self, InPack, OutPack):
        """
        This fxn does temporal alignment of input datasets. It pulls the 
        the expected file name from the `expectations` object and writes a 
        csv to that location. 
        """
        out_fn = self.glm_config['outflow']['outflow_fl'][1:-1]
        in_fn = self.glm_config['inflow']['inflow_fl'][1:-1]
        self.incsv_path = os.path.join(self.config_dir, in_fn)
        self.outcsv_path = os.path.join(self.config_dir, out_fn)
        self.out_df = pd.DataFrame(index= self.date_range, columns= OutPack)
        self.in_df = pd.DataFrame(index= self.date_range, columns= InPack)
        
        for i in InPack:
            self.in_df[i] = self.aligned_columns[i]
        for i in OutPack:
            self.out_df[i] = self.aligned_columns[i]
            
        self.out_df.rename(columns = {'OUTFLOW':'FLOW'}, inplace = True)
        self.in_df.rename(columns = {'INFLOW':'FLOW'}, inplace = True)
        
        print "\tWriting inflow.csv"
        self.in_df.to_csv(path_or_buf=self.incsv_path, index_label='time')
        print "\tWriting outflow.csv"
        self.out_df.to_csv(path_or_buf=self.outcsv_path, index_label='time')

        
    def read_variants(self, path, verbose=True):
        meta_config = self.glm_config
        block_key_option = []
        variants = {'grid': {},
                    'prior': {},
                    'likelihood': {},
                    'posterior': {} }
        
        
        with open(path, 'r') as f:
            for line in f:
                if line[0] == '#':
                    pass
                elif line[0] == '&':
                    new_block = line[1:].strip()
                    assert new_block in meta_config.keys()
                    if verbose == True:
                        print "\n", new_block
                elif line[0] == '-':
                    splitted = line[2:].split('=')
                    spit_cleaned = [i.strip() for i in splitted]
                    sub_block = spit_cleaned[0]
                    variants['grid'][sub_block] = []
                    if verbose == True:
                        print "\t%r" % sub_block
                        
                    [val_space, dType] = spit_cleaned[1].split("],[")
                    val_space, dType = val_space[1:], dType[:-1]
                    val_split = [l.strip() for l in val_space.split(",")]
                    
                    if val_split[-1] == 'None':
                        val_split.pop()
                        
                    if dType == 'float':
                        val_split = [float(l) for l in val_split]
                    elif dType == 'int':
                        val_split = [int(l) for l in val_split]
                            
                    if dType == 'bool':
                        for idx, i in enumerate(val_split):
                            if verbose == True:
                                print "\t\t%r" % i
                            sim_name = sub_block+"_"+str(idx+1)
                            block_key_option.append((new_block, sub_block, i, 
                                                     sim_name))
                            variants['grid'][sub_block].append(i)
                                                     
                    elif len(val_split) != 3 and dType != 'bool':
                        for idx, i in enumerate(val_split):
                            sim_name = sub_block+"_"+str(idx+1)
                            block_key_option.append((new_block, sub_block, i, 
                                                     sim_name))
                            if verbose == True:
                                print "\t\t%r" % i
                            variants['grid'][sub_block].append(i)
                            
                    elif len(val_split) == 3:
                        arange = np.arange(val_split[0], val_split[1]+val_split[2],
                                           val_split[2])
                        for idx, i in enumerate(arange):
                            j = int(np.floor(i*1000))
                            j = j/1000.
                            
                            if verbose==True:
                                print "\t\t%r" % j

                            sim_name = sub_block+"_"+str(idx+1)
                            block_key_option.append((new_block, sub_block, j, 
                                                     sim_name))
                            variants['grid'][sub_block].append(j)
                    else:
                        print "\n\nunparsed option\n\n"
                                                
                else:
                    print "Unparsable line"
                    print line
                    
        print "\t%i variant cases specified" % len(block_key_option)
        
        self.variant_cases = block_key_option
        
        for k in variants['grid'].keys():
            len_ = len(variants['grid'][k])
            variants['grid'][k] = np.array(variants['grid'][k])
            variants['prior'][k] = np.ones(len_)
            variants['likelihood'][k] = np.zeros(len_)
            variants['posterior'][k] = np.zeros(len_)
            
        self.variants = variants

    def write_variant_configs(self, copycsvs=False, verbose=True):
        
        if self.glm_config and self.variant_cases:
            variant_lakes = []
            for i in self.variant_cases:
                variant_lake = copy.deepcopy(self)
                key1, key2, val, simname = i[0],i[1],i[2],i[3]
                variant_lake.glm_config[key1][key2] = val
                variant_lake.name = simname
                variant_lake.variant_var = key2
                variant_lake.variant_val = val
                
                variant_lake.ran = False
                
                variant_lake.dir_path = os.path.join(self.dir_path, simname)
                make_dir(variant_lake.dir_path, verbose)
                variant_lake.glm_path = os.path.join(variant_lake.dir_path,
                                                     "glm2.nml")
                variant_lake.write_glm_config(verbose)
                
                if copycsvs==True:
                    if verbose==True:
                        print "Copying original data csvs to variant folders"
                    shutil.copy(self.metcsv_path, variant_lake.dir_path)
                    shutil.copy(self.outcsv_path, variant_lake.dir_path)
                    shutil.copy(self.incsv_path, variant_lake.dir_path)
                else:
                    if verbose==True:
                        print "Original data csvs not moved into variant folders"
                    
                variant_lakes.append(variant_lake)
        return variant_lakes
    
    

    def score_variants(self, observations, obs_type, variant_list, lik_fxn='NSE'):
        
        assert self.ran == True
        
        for v_l in variant_list:
            # ensure input object is a Lake model that was run and has had its 
            # data collected. 

            assert type(v_l) == Lake
            assert v_l.ran == True
            assert v_l.csv_dict

            #is top bottom orientation correct??
        if obs_type in self.depth_dfs.keys():
            baseline = self.depth_dfs[obs_type].resample('D').mean()
            self.observed_dm = time_depth_df_filter(baseline, observations)
            self.error = error_metrics(observations.values, self.observed_dm)
            
            for variant in variant_list:
                # Step 2: estimate the liklihood 
                sub_block_ = variant.variant_var
                curr_val = variant.variant_val
                all_idxs = np.where(self.variants['grid'][sub_block_] == curr_val)
                idx = all_idxs[0]
                
                modelled = variant.depth_dfs[obs_type].resample('D').mean()
                variant.observed_dm = time_depth_df_filter(modelled, 
                                                           observations)
                
                err = error_metrics(observations.values, variant.observed_dm)
                self.variants['likelihood'][sub_block_][idx] = err[lik_fxn]
                
            for sb2 in self.variants['likelihood'].keys():
                if len(self.variants['grid'][sb2]) > 4:
                    raw_post=self.variants['prior'][sb2]*self.variants['likelihood'][sb2]
                    evidence=np.trapz(raw_post, self.variants['grid'][sb2])
                    self.variants['posterior'][sub_block_]=raw_post/evidence
                
            
def time_depth_df_filter(baseline, observations):
    #extract proper baseline here based on real world observations
    time_matched = baseline[baseline.index.isin(observations.index)]
    observed_dm = np.zeros(observations.shape)
    obs_days, depths = observations.shape
    for row in range(obs_days):
        this_row = time_matched.iloc[row, :].dropna()
        this_len = len(this_row)-1
        binsize = int(np.round(this_len/float(depths)))
        for depth, elem in zip(range(depths),range(0,this_len, binsize)):
            observed_dm[row, depth] = this_row[elem:elem+binsize].mean()
    return observed_dm

    
class USGS_water_data(object):

    def __init__(self, fn):
        self.fn = fn
    
    def preprocess(self):
        f = open(self.fn)    
        s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        self.alt1 = False  
        
        try:
            header = s.find("#    DD parameter   Description")
            s.seek(header)
        except ValueError:
            header = s.find("#    DD parameter statistic")
            s.seek(header)
            self.alt1=True            

        s.seek(header)
        s.readline()
        text = s.readline()
        words = [splits for splits in text.split(" ") if splits is not ""]
        self.words = words
        dds, params, descrips  = [], [], []
        
        while words[0] != "#\n":
            dds.append(words[1])
            
            if self.alt1 == False:
                params.append(words[2])
                descrips.append(" ".join(words[3:])[:-1])
            else:
                params.append(words[2]+"_"+words[3])
                descrips.append(" ".join(words[4:])[:-1])
                
            text = s.readline()
            words = [splits for splits in text.split(" ") if splits is not ""]
        s.seek(0)        
        line = s.readline()
        com_lines = 0
        while line[0] == '#':
            line = s.readline()
            com_lines +=1
            
        self.cols = line
        re.split(r'\t+', self.cols)
        self.skipRows = []
        self.com_lines, counter2 = com_lines, com_lines
        
        while line[:4]!='USGS':
            line = s.readline()
            counter2+=1
            self.skipRows.append(counter2)
            
        s.close()
        f.close()
        
        col_keys = [i+"_"+j for i, j in zip(dds, params)]            
        self.col_dict = {i:j for i, j in zip(col_keys, descrips)}
        self.qual_cols = [i+"_cd" for i in col_keys]
        
    def read(self):
        self.df = pd.read_csv(self.fn, sep="\t", header=self.com_lines, 
                                  skiprows=self.skipRows, index_col=[2], 
                                  parse_dates=[2], low_memory=False,
                                  infer_datetime_format=True)
        
        new_cols = []
        for vals in self.df.columns:
            if self.col_dict.has_key(vals):
                new_cols.append(self.col_dict[vals])
            else:
                new_cols.append(vals)
        self.df.columns = new_cols

    def print_metadata(self):
        print "\n Metadata for %r" %os.path.basename(self.fn)
        
        for i in self.df.columns:
            if i[-3:] == '_cd':
               qual_col = self.df[i]
               print i
               print "\t%r" % qual_col.value_counts().index
               print "\t%r" % qual_col.value_counts().values
        print ""
        
    def removeQuals(self, verbose=True):
        for column in self.df.columns:
            if '_cd' in column:
                self.df.drop(column, axis=1, inplace=True)
            elif 'site_' in column:
                self.df.drop(column, axis=1, inplace=True)
            else:
                if verbose==True:
                    print "{0} is a valid column".format(column)
                
class GHCN_weather_data(object):
    
    def __init__(self, met_path):
        """this creates an dataframe object called `raw met` that contains all
        information within the csv provided at `met_path`"""
        
        self.met_path = met_path
        self.raw_met = pd.read_csv(met_path, parse_dates=[2], 
                                   infer_datetime_format=True, 
                                   na_values=[-9999])

    def clean_columns(self):
        """
        This converts precipitation from tenths of mm to meters. Windspeed is
        converted from tenths of meters per second to meters per second. 
        Max & min temperatures are converted to celsius. Snowfall is converted 
        to meters. This also modifies the column names.
        """
        #PRCP - Precipitation (tenths of mm)        
        self.precip = self.raw_met.PRCP/10000
        
        #AWND - Average daily wind speed (tenths of meters per second)        
        self.wind_speed = self.raw_met.AWND/10
        
        #TMAX - Maximum temperature (tenths of degrees C)
        self.t_max = self.raw_met.TMAX/10
        self.t_min = self.raw_met.TMIN/10
        #SNOW - Snowfall (mm) -> meters
        self.snow_fall = self.raw_met.SNOW/1000
        #SNWD - Snow depth (mm) -> meters
        self.snow_depth = self.raw_met.SNWD/1000
        self.clean_cols = [self.precip, self.wind_speed, self.t_max, self.t_min, 
                      self.snow_fall, self.snow_depth]
        col_labels = ["precip", "wind_speed", "t_max", "t_min", "snow_fall",
                      "snow_depth"]
        self.df = pd.DataFrame(index=self.raw_met.DATE, columns=col_labels)

        for idx in range(len(self.clean_cols)):
            self.df.iloc[:, idx] = self.clean_cols[idx].values
        
        self.gchn_zs = z_score(self.df)
        #printNaNWarning(self.df, "GCHN data", None)

                                           
class CERES_nc(object):
    
    def __init__(self, ceres_path, type):
        if type == "3":
            rootgrp = Dataset(ceres_path, "r", format="NETCDF3_CLASSIC")
            nc_attrs, nc_dims, nc_vars = ncdump(rootgrp, False)
            # Temps in Kelvin        
            # Fluxes in u'Watts per square meter'
            # wind vectors in meters per second
            # humididity & cloud cover in ?
            toExtract = ["lon", "lat", "Time_of_observation", 
                         "Surface_skin_temperature", "Surface_wind___U_vector",
                         "Surface_wind___V_vector", 
                         "Column_averaged_relative_humidity", 
                         "Altitude_of_surface_above_sea_level",
                         "PSF_wtd_MOD04_cloud_fraction_land",
                         "CERES_net_LW_surface_flux___Model_B",
                         "CERES_net_SW_surface_flux___Model_B"]
        
            betterNames = ["lat", "lon", "time_J", "temp", "windU", "windV",
                       "humidity", "altitude", "cloud_frac", "LW_rad", "SW_rad"]
            self.ceres_ = {}
            
            for var, nam in zip(toExtract, betterNames):
                self.ceres_[nam] = pd.Series(rootgrp.variables[var][:])
            gDays = [jd.caldate(j) for j in self.ceres_["time_J"]]
            self.ceres_['date'] = pd.to_datetime(np.array(gDays))
            
            self.ceres_df = pd.DataFrame(index=self.ceres_['date'], 
                                         columns=betterNames)
    
            ignorable = ['lat', 'lon', 'time_J', 'altitude']
    
            for key in self.ceres_.keys():
                if key == 'temp':
                    self.ceres_df[key] = self.ceres_[key].values-273.15
                elif key in ignorable:
                    self.ceres_df = self.ceres_df.drop(key, 1)
                elif key not in self.ceres_df.columns.values:
                    pass
                else:
                    self.ceres_df[key] = self.ceres_[key].values
            
            self.ceres_zs = z_score(self.ceres_df)
            self.ceres_d_zs = self.ceres_zs.resample("D").mean()
            self.ceres_d_zs = self.ceres_d_zs.ix[:-1]        
            #printNaNWarning(self.ceres_df, "CERES data", None)
            self.ceres_reindex_zs, self.i_names = TimeIdx(self.ceres_d_zs)
            self.ceres_mean_aggs = makeAggregations(self.ceres_reindex_zs, 
                                                    self.i_names, np.mean)
                                                    
            self.humidity_scales = show_agg_resolution(self.ceres_mean_aggs, 
                                                       "humidity")
            self.cloud_scales = show_agg_resolution(self.ceres_mean_aggs, 
                                                       "cloud_frac")
            #printAutocorr(self.ceres_mean_aggs)
            #pass dataframe to insert new index columns
            rootgrp.close()
        elif type == "4":
            import xarray as xr
            data = xr.open_dataset(ceres_path)
            vars = data.keys()
            self.lats_n = data['lat'].values
            self.lons_n = data['lon'].values
            self.time_n = data['time'].values
            other_vars = vars[5:]
            time_ns, lats_ns, lons_ns = data[vars[3]].values.shape
            
            row1 = data[vars[3]].values[:,0,0]
            row2 = data[vars[4]].values[:,0,0]
            stack = np.vstack((row1, row2))
            for var in other_vars:
                row3 = data[var].values[:,0,0]    
                stack = np.vstack((stack, row3))
            
            self.df1 = pd.DataFrame(index=self.time_n, columns=vars[3:], 
                                    data = stack.T.copy())
            self.df1.aux_skint_daily = self.df1.aux_skint_daily - 273.15
                
            
            row1 = data[vars[3]].values[:,0,1]
            row2 = data[vars[4]].values[:,0,1]
            stack = np.vstack((row1, row2))
            for var in other_vars:
                row3 = data[var].values[:,0,1]
                stack = np.vstack((stack, row3))
            
            self.df2 = pd.DataFrame(index=self.time_n, columns=vars[3:], 
                                    data = stack.T.copy())
                
def run_model(Lake, ex_loc = None, verbose=True):
    """This function accesses the locations of each glm configuration, 
    provided in Lake.glm_path, and runs the model using the executable
    specefied in `ex_loc`
    """
    report_path = os.path.join(Lake.dir_path, "Run_report.txt")
    unique_divider = "\n*\n"
    
    if os.path.exists(report_path):
        if verbose:
            print "Prior run of {} model detected".format(Lake.name)        
        report_h = open(report_path, 'r')
        report_ = report_h.read().split(unique_divider)
        report_.pop()
        check3 = report_[-2] == Lake.variant_var
        check4 = report_[-1] == str(Lake.variant_val)
        
        if check3 and check4:
            Lake.ran, Lake.stderr, Lake.stdout = True, report_[-3], report_[-4]
            if verbose:
                print "Prior run of {} model accepted\n".format(Lake.name)
        else:
            report_h.close()
            if verbose:
                print "Prior run rejected. Re-running\n"
                print "Expected: %r, %r" % (Lake.variant_var, str(Lake.variant_val))
                print "Recieved: %r, %r" % (report_[-2], report_[-1])

            os.remove(report_path)
            Lake = run_model(Lake)  
    else: 
        if verbose:
            print "No prior run of model detected @ {}".format(report_path)
            
            
        if ex_loc == None:
            curr_dir = os.path.dirname(os.getcwd())
            ex_loc = 'GLM_Executables/glm.app/Contents/MacOS/glm'
            ex_loc = os.path.join(curr_dir, ex_loc)
            
        p = sp.Popen(ex_loc, cwd=Lake.dir_path, shell=True, stderr=sp.PIPE, 
                     stdout=sp.PIPE)
        p.wait()
        Lake.stdout, Lake.stderr = p.communicate()
        Lake.ran = True
        
        check1 = 'Simulation begins..' in Lake.stdout
        check2 = '100.00% of days complete' in Lake.stdout
        if check1 and check2:
            report_h = open(report_path, 'w')
            report_h.writelines("Run Report"+unique_divider)
            
            report_h.writelines(Lake.stdout+unique_divider)
            
            report_h.writelines(Lake.stderr+unique_divider)
            
            report_h.writelines(Lake.variant_var+unique_divider)
            
            report_h.writelines(str(Lake.variant_val)+unique_divider)
            
            if verbose:
                print "{} model run completed\n".format(Lake.name)
    
    return Lake

@fn_timer
def pull_output_nc(Lake, plot_var=[], verbose=False, plot=True):
    if Lake.ran == True:
        
        o_fn = Lake.glm_config['output']["out_fn"][1:-1] + '.nc'
        o_path = Lake.glm_config['output']["out_dir"][1:-1]
        
        o_mult = Lake.glm_config['output']['nsave']
        o_interval = Lake.glm_config['time']['dt']
        freq_S = str(int(o_mult*o_interval))+"S"
        if verbose == True:
            print "The time steps in the output netCDF are %s" % freq_S
        
        interval = pd.date_range(start=Lake.glm_config['time']['start'], 
                                 end=Lake.glm_config['time']['stop'], 
                                 freq=freq_S)
        Lake.nc_interval = interval[1:]                         
        Lake.output_nc_path = os.path.join(o_path, o_fn)
        output = Dataset(Lake.output_nc_path, "r")
        nc_vars = output.variables.keys()
        
        if 'time' in nc_vars:
            t_vector = output['time'][:]
            start_t = pd.to_datetime(output['time'].units[12:])
            d_vector = [start_t + DateOffset(hours=t) for t in t_vector]
        
        assert len(d_vector) == len(Lake.nc_interval)
        if verbose == True:
            print "THe netCDF time vector is the expected size"
        
        lake_cols, lake_vecs, depth_dfs = [], {}, {}
        
        for ncv in nc_vars:
            dimens = len(output[ncv].shape)
            if dimens == 3:
                lake_cols.append(ncv)
                lake_vecs[ncv] = output[ncv][:,0,0]
            if dimens == 4:
                data_buff = output[ncv][:,:,0,0]
                data_mat = data_buff.data
                fill_val = data_buff.fill_value
                data_mat[data_mat == fill_val] = np.nan
                col_n = range(data_mat.shape[1])
                this_df = pd.DataFrame(columns = col_n, index = d_vector, 
                                       data = data_mat)
                assert data_buff.mask.sum() == this_df.isnull().sum().sum()
                depth_dfs[ncv] = this_df

        Lake.depth_dfs = depth_dfs
        Lake.lake_df = pd.DataFrame(columns=lake_cols,
                                    index=d_vector,
                                    data=lake_vecs)
        if verbose == True:
            print "%r columns added to the lake_df from the netCDF" % lake_cols
        
        loc_files = [i for i in os.listdir(Lake.dir_path) if '.csv' in i]
        overfl_fname = Lake.glm_config['output']['csv_ovrflw_fname'][1:-1]
        lake_fname = Lake.glm_config['output']['csv_lake_fname'][1:-1]
        outlet_fname = Lake.glm_config['output']['csv_outlet_fname'][1:-1]
        
        good_files = []
        for file_ in loc_files:
            if overfl_fname in file_:
                good_files.append(file_)
            if lake_fname in file_:
                good_files.append(file_)
            if outlet_fname in file_:
                good_files.append(file_)
        #TODO: add WQ_ files here when ready
        if verbose == True:
            print "%r loaded from the configuration folder" % good_files
        
        hdls = [os.path.join(Lake.dir_path, i) for i in good_files]
        csv_dict = {}
        for p, f in zip(hdls, good_files):
            csv_dict[f] = pd.read_csv(p)

        for f in good_files:
            t_steps, lake_vs = csv_dict[f].shape
            csv_dict[f]['date'] = csv_dict[f].time.copy()
            csv_dict[f].time = csv_dict[f].time.apply(lambda x: x[-8:])
            csv_dict[f].date = csv_dict[f].date.apply(lambda x: x[:10])
            time_bool = csv_dict[f].loc[:, 'time'] == '24:00:00'
            csv_dict[f].loc[time_bool, 'time']  = '23:59:59'
            csv_dict[f].datetime = csv_dict[f].date + " "+ csv_dict[f].time
            csv_dict[f].index = pd.to_datetime(csv_dict[f].datetime)                
            csv_dict[f].drop(['date', 'time'], axis=1, inplace=True)

        for k in csv_dict.keys():
            old_cols = list(csv_dict[k].columns)
            if overfl_fname in k:
                new_cols = ["overflow_"+i for i in old_cols]
            elif outlet_fname in k:
                new_cols = ["outflow_"+i for i in old_cols]
            else:
                new_cols = old_cols
            csv_dict[k].columns = new_cols


        Lake.csv_dict = csv_dict

        if verbose == True:
            for k in output.variables.keys():
                print k
                
                if k != 'NS':
                    print "\t{}".format(output[k].units)
                    print "\t{}".format(output[k].shape)
                    
        if plot_var !=[] and plot==True:
            for i, _var in enumerate(plot_var):
                plt.figure()                    
                if len(output[_var].shape) == 4:
                    temp = output[_var][:,:,0,0].data
                    temp[temp == output[_var]._FillValue] = -1
                    plt.imshow(np.flipud(temp.T))
                    plt.colorbar()
                elif len(output[_var].shape) == 3:
                    temp = output[_var][:,0,0].astype(float)
                    temp[temp == output[_var]._FillValue] = np.nan
                    plt.plot(Lake.nc_interval, temp)
                
                plt.title(_var)                        
                plt.show()
        
                # the default value assigned to masked elements > 1e36, 
                # here we reassign them -1 for plotting purpose
                #np.flipud puts the lake bottom at the plot bottom 
                
                    
    else:
        print "Can't plot before running model"
    
    Lake.output_nc = output
    return Lake