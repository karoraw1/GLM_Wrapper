# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:03:39 2016
ncdump() was created by this dude:
http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html#code
JD_converter was created by this dude: (and modified by me)
https://asimpleweblog.wordpress.com/2010/06/20/julian-date-calculator/

@author: login
"""

from parse_aed_config import read_aed_nml, read_cell_params, read_zoo_params
from parse_aed_config import write_aed2_eco_configs, write_aed2_config
from parse_aed_config import load_aed_var_types

import scipy.stats as st
import re, mmap, cPickle, os, copy, time, sys, shutil
import numpy as np
import pandas as pd
from pandas.tseries.offsets import DateOffset
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import JD_converter as jd
import subprocess as sp
from functools import wraps
from mpl_toolkits.axes_grid1 import make_axes_locatable

def preheim_date_parser(date_str):
    date_part = date_str.split("_")[1]
    new_str = date_part[2:4] + "-" + date_part[4:6] + "-" + date_part[0:2]
    return pd.to_datetime(new_str)

def readInElisasAndPreheimData(chem_dir):
    do_eliza_fn = os.path.join(chem_dir, "Elizas_DO_data_mg.csv")
    ph_eliza_fn = os.path.join(chem_dir, "Elizas_PH_data.csv")
    do_preheim_fn = os.path.join(chem_dir, "Preheim_DO_mg.tsv")
    no3_preheim_fn = os.path.join(chem_dir, "Preheim_NO3_mg.tsv")
    temp_preheim_fn = os.path.join(chem_dir, "Preheim_TEMP_C.tsv")
    
    do_eliza_df = pd.read_csv(do_eliza_fn, index_col=0, dtype=float)
    do_eliza_df.columns = [pd.to_datetime(i) for i in do_eliza_df.columns]
    do_eliza_df = do_eliza_df.T
    
    ph_eliza_df = pd.read_csv(ph_eliza_fn, index_col=0, dtype=float)
    ph_eliza_df.columns = map(pd.to_datetime, ph_eliza_df.columns)
    ph_eliza_df = ph_eliza_df.T
    
    do_preheim_df = pd.read_csv(do_preheim_fn, sep="\t", index_col=0, dtype=float)
    do_preheim_df.columns = map(preheim_date_parser, do_preheim_df.columns)
    do_preheim_df.interpolate(axis=0, inplace=True)
    surf_null_mask_2 = do_preheim_df.ix[0, :].isnull().values
    do_preheim_df.ix[0, surf_null_mask_2] = do_preheim_df.ix[1, surf_null_mask_2]
    do_preheim_df = do_preheim_df.T
    idx_to_drop = do_preheim_df.index[5]

    print "dropping this {}".format(idx_to_drop)
    do_preheim_df.drop(idx_to_drop, axis=0, inplace=True)
    
    no3_preheim_df = pd.read_csv(no3_preheim_fn, sep="\t", index_col=0, dtype=float)
    no3_preheim_df.columns = map(preheim_date_parser, 
                                 no3_preheim_df.columns)
    no3_preheim_df.interpolate(axis=0, inplace=True)
    surf_null_mask_1 = no3_preheim_df.ix[0, :].isnull().values
    no3_preheim_df.ix[0, surf_null_mask_1] = no3_preheim_df.ix[1, surf_null_mask_1]
    no3_preheim_df = no3_preheim_df.T
    
    temp_preheim_df = pd.read_csv(temp_preheim_fn, sep="\t", index_col=0, 
                                  dtype=float)
    temp_preheim_df.columns = map(preheim_date_parser, 
                                  temp_preheim_df.columns)
    temp_preheim_df.ix[1, '2012-11-12'] = temp_preheim_df.ix[2, '2012-11-02']
    surf_null_mask_0 = temp_preheim_df.ix[0, :].isnull().values
    temp_preheim_df.ix[0, surf_null_mask_0] = temp_preheim_df.ix[1, surf_null_mask_0]
    temp_preheim_df.interpolate(axis=0, inplace=True)
    temp_preheim_df = temp_preheim_df.T
    
    new_cols = list(temp_preheim_df.columns)
    new_cols.reverse()
    temp_preheim_df.columns = new_cols
    temp_preheim_df.columns = [int(i) for i in temp_preheim_df.columns]
    no3_preheim_df.columns = new_cols
    no3_preheim_df.columns = [int(i) for i in no3_preheim_df.columns]
    do_preheim_df.columns = new_cols
    do_preheim_df.columns = [int(i) for i in do_preheim_df.columns]

    ph_eliza_df.columns = new_cols
    ph_eliza_df.columns = [int(i) for i in ph_eliza_df.columns]
                           
    do_eliza_df.columns = new_cols
    do_eliza_df.columns = [int(i) for i in do_eliza_df.columns]
                               
    obs_pack = [do_eliza_df, ph_eliza_df, do_preheim_df, no3_preheim_df, 
                temp_preheim_df]
    return obs_pack
    
    
def cond2sal(sal_ser):
    return 10**((np.log10(sal_ser.copy())*0.909)-2.823)
    
def CubicftPerS_to_MegalitersPerDay(df):
    return df*2.44658

def print_minimums(Lake):
    var_dict = Lake.variants
    for k in var_dict['grid'].keys():
        print "Variable:", k
        liks = var_dict['likelihood'][k]
        liks[np.isnan(liks)] = -1
        best_lik = liks.max()
        best_idxes = np.where(liks == best_lik)[0]
        best_idx = best_idxes[0]
        best_val = Lake.variants['grid'][k][best_idx]

        for k2, v2 in Lake.glm_config.items():
            if k in v2.keys():
                orig_val = v2[k]
        

        print "\t Original Value: %r" % orig_val
        if type(orig_val) == str:
            orig_idxes = np.where(var_dict['grid'][k] == orig_val)[0]
        elif type(orig_val) == list:
            orig_idxes = np.where(np.isclose(var_dict['grid'][k],1.0, rtol=5e-02))[0]
        else:                
            orig_idxes = np.where(np.isclose(var_dict['grid'][k],orig_val, rtol=5e-02))[0]
            
        if len(orig_idxes) == 0:
            orig_lik = liks[1:2].mean()
        else:
            orig_lik = liks[orig_idxes][0]
            
        if liks[best_idx] == orig_lik:
            print "Variable", k, "was set at thet best setting"            
        else:
            print "\t Best Value: %r" % best_val
            print "\t Likelihood: %.2f" % best_lik
            print "\t Relative diff:", abs(orig_lik-best_lik)/abs(orig_lik)
            if len(best_idxes) > 1:
                    print "\tMultiple peaks at:", best_idxes
                    for i in best_idxes:
                        print liks[i]

def plotLakeandError(scored_Lake, Lake_temps, fignum):
    modelled_dm = scored_Lake.modelled_mat
    plt.figure(fignum, figsize=(16,8))
    titles = ["Observed", "Model", "Error"]
    ylabs_ = [str(i.date()) for i in Lake_temps.index]
    for i, v in enumerate([Lake_temps, modelled_dm, abs(Lake_temps-modelled_dm)]):
        ax = plt.subplot(1, 3, i+1)
        ax1 = plt.imshow(v)
        plt.title(titles[i])
        if i == 0:
            plt.yticks(np.arange(-2,18,2), ylabs_)
        else:
            plt.yticks(np.arange(-2,18,2), [""]*len(ylabs_))
        if i == 1:
            plt.xlabel("Depth in meters")
        divider = make_axes_locatable(ax)
        cax1 = divider.append_axes("right", size="5%", pad="3%")
        plt.colorbar(ax1, cax=cax1)
    error = abs(Lake_temps-modelled_dm).sum().sum()
    denom = (Lake_temps.shape[0]*Lake_temps.shape[1])
    print error/denom, "average error per voxel"

def plotLakeProfiles(scored_lake, fignum):
    fig = plt.figure(fignum, figsize=(18,9))
    ax = fig.add_subplot(3,1,1)
    ax.set_title("Elevation", fontsize=18)
    cax = ax.imshow(np.flipud(scored_lake.depth_dfs['z'].iloc[:, :150].T))
    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes("right", size="5%", pad="3%")
    plt.colorbar(cax, cax=cax1)
    
    ax = fig.add_subplot(3,1,2)
    ax.set_title("Temperature", fontsize=18)
    ax.set_ylabel("Layer Number", fontsize=16)
    cax = ax.imshow(np.flipud(scored_lake.depth_dfs['temp'].iloc[:, :150].T))
    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes("right", size="5%", pad="3%")
    plt.colorbar(cax, cax=cax1)
    
    ax = fig.add_subplot(3,1,3)
    ax.set_title("Salinity", fontsize=18)
    ax.set_xlabel("Time Steps (12 hour)", fontsize=16)
    cax = ax.imshow(np.flipud(scored_lake.depth_dfs['salt'].iloc[:, :150].T))
    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes("right", size="5%", pad="3%")
    plt.colorbar(cax, cax=cax1)
 
def hard_round(x):
    return int(x*1000.0)/1000.0
    
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
    
def import_glm_config(dl_default=False, verbose=False):
    # the glm contains the following blocks:
    # &glm_setup: 13 general simulation info and mixing parameters
    print "Importing AED baseline config"
    if dl_default:
        default_config = {'others':{}}
        default_order = {}
        block_order = []
        with open('../test_files/sample_aed_glm.nml', 'r') as f:
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
        
    glm_setup = {'max_layers' :300, 
                 'min_layer_vol' :0.02, 
                 'min_layer_thick' :0.32, 
                 'max_layer_thick' :1.0,  
                 'Kw' : 0.4, 
                 'coef_mix_conv' : 0.1125, 
                 'coef_wind_stir' : 0.1625, 
                 'coef_mix_shear' : 0.21, 
                 'coef_mix_turb' : 0.6375, 
                 'coef_mix_KH' : 0.315, 
                 'coef_mix_hyp' : 0.675,
                 'deep_mixing':'.false.'}
                 
    wq_setup = {'wq_lib': "'aed2'",
                'wq_nml_file': "'aed2.nml'",
                'ode_method': 11,
                'split_factor': 1,
                'bioshade_feedback': '.true.',
                'repair_state': '.true.',
                'benthic_mode': 1}
    
    # &morphometry: 8 vars
    morphometry = {'lake_name': "'UpperMysticLake'", 
                   "latitude": 42,  
                   "longitude": -71,
                   "bsn_len": 1595.0,
                   "bsn_wid": 570.0, 
                   "bsn_vals": None,
                   #"A":[ 77373.80, 148475.73, 202472.97, 257818.95, 338552.69,
                   #    397077.50, 460778.04, 524802.66, 560051.22]
                   "A":[ 77373.80, 148475.73, 202472.97, 257818.95, 338552.69,
                       397077.50, 460778.04, 629763.19, 784071.71]
                   
                        }
    depths_ = [ 0., 3.048, 6.096, 9.144, 12.192, 15.24, 18.288, 21.336,  24.384]
    depths_.reverse()
    morphometry['H'] = [i*-1. for i in depths_]

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
    time['num_days'] = (end_t-start_t).days + 1
    
    # &output 14 vars 
    output = {"out_dir" : None, # set at write time
              "out_fn" : None,  # set at write time 
              "nsave" : 24, 
              "csv_lake_fname" : "'lake'",
              "csv_point_nlevs" : None,
              "csv_point_fname" : "'WQ_'",
              "csv_point_at" :  list(np.arange(1,21,3)) + [21],
              "csv_point_nvars" : None,
              "csv_point_vars" : ['temp', 'salt', 'OXY_oxy', 'CAR_dic', 'CAR_pH',
                                  'CAR_ch4', 'SIL_rsi', 'NIT_amm', 'NIT_nit',
                                  'PHS_frp', 'OGM_don', 'OGM_pon', 'OGM_dop',
                                  'OGM_pop','OGM_doc','OGM_poc','PHY_green'],
              "csv_outlet_allinone" : ".false.",
              "csv_outlet_fname" : "'outlet_'",
              "csv_outlet_nvars" : None,
              "csv_outlet_vars" : None,
              "csv_ovrflw_fname" : "\"overflow\"" }

    output["csv_outlet_vars"] = ['flow']+output['csv_point_vars']
    output["csv_outlet_nvars"] = len(output["csv_outlet_vars"])
    output["csv_point_nlevs"] = len(output["csv_point_at"])
    output["csv_point_nvars"] = len(output["csv_point_vars"])
    
    
    #&init_profiles            
    init_profiles = {"lake_depth": 24.384, 
                     "num_depths": 6, 
                     "the_depths": [1, 5, 9, 13, 17, 21], 
                     "the_temps": [4.0, 4.0, 4.0, 4.0, 4.0, 4.0], 
                     "the_sals": [0.549, 1.10, 1.492, 1.937, 2.373, 2.801]}
    
    init_profiles['num_depths'] = len(init_profiles["the_depths"])
    init_profiles['the_sals'] = [(i*1.8*1.945*1.91) for i in init_profiles['the_sals']]
    assert len(init_profiles["the_depths"]) == len(init_profiles["the_temps"])
    assert len(init_profiles["the_depths"]) == len(init_profiles["the_sals"])
    
    init_profiles["wq_names"] = ['OGM_don', 'OGM_pon', 'OGM_dop', 'OGM_pop',
                                 'OGM_doc', 'OGM_poc']
    init_profiles["num_wq_vars"] = len(init_profiles["wq_names"])
    #init_profiles["wq_init_vals"] = [.5, .5, .5, .5, .5, .5]*6
    init_profiles["wq_init_vals"] = [.5, .5, .5, .5, .5, .5]*6

    #&meteorology                 
    meteorology = {"met_sw" : ".true.",
                   "lw_type" : "'LW_IN'", 
                   "rain_sw" : ".false.",
                   "atm_stab" : ".false.", 
                   "catchrain" : ".false.",
                   "rad_mode": 0,
                   "albedo_mode": 3, 
                   "cloud_mode": 3, 
                   "subdaily": '.false.', 
                   "meteo_fl": "'met_daily.csv'",
                   "wind_factor": 0.7,
                   "sw_factor":  1.0,
                   "lw_factor": 1.1, 
                   "at_factor": 1.0,
                   "rh_factor": 1.0,
                   "rain_factor": 1.00,
                   "ce":  0.001, 
                   "ch": 0.0011, 
                   "cd": 0.0015, 
                   "rain_threshold": 0.03, 
                   "runoff_coef":  0.1}
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
              "strm_hf_angle": 1.1, 
              "strmbd_slope": 0.5, 
              "strmbd_drag": 0.008, 
              "inflow_factor": 0.55, 
              "inflow_fl": "'inflow.csv'", 
              "inflow_varnum": None, 
              "inflow_vars": ['FLOW','TEMP','SALT','OXY_oxy','SIL_rsi','NIT_amm',
                              'NIT_nit','PHS_frp','OGM_don','OGM_pon','OGM_dop',
                              'OGM_pop','OGM_doc','OGM_poc','PHY_green',
                              'PHY_crypto','PHY_diatom', 'CAR_dic', 'CAR_pH'],
              "coef_inf_entrain": 0.7 }
    inflow["inflow_varnum"] = len(inflow['inflow_vars'])
    #&outflow
    outflow = {"num_outlet": 1,
               "flt_off_sw": ".false.",
               "outl_elvs": 0.,
               "bsn_len_outl": 100.0, 
               "bsn_wid_outl" : 274.0,
               "outflow_fl" : "'outflow.csv'",
               "outflow_factor": 0.5,
               "seepage" : ".false.",
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

def calculate_vector_space_by_midpoint(midpoints, increment_n):
    """
    This function takes a list of values and prints the inputs for np.arange()
    such that a vector of equidistant points is created around each value. 
    If the boundaries of the space fall into negative territory, the entire
    space is shifted in the positive direction until the lower limit is 
    equal to the midpoint[x]*(1/increment_n)
    """
    #midpoints = [0.4, 0.1125, 0.1625, 0.21, 0.6375, 0.315, 0.675, 1595.0, 
    #             570.0, 1.0, 1.0, 100.0, 274.0, 0.7, 0.008, 0.5, 1.1, 0.55, 
    #             0.03, 0.7, 1.00, 1.0, 1.0, 1.0, 1.1, 0.0015, 0.001, 0.0011]
    #opt_space = 200
    opt_space = increment_n
    pct_inc = (1/float(opt_space))
    for i in midpoints:
        increment = i*pct_inc
        low_point = i - increment*(opt_space/2)
        high_point = i + increment*(opt_space/2 + 1)
        if low_point < 0:
            neg_dist = 0 - low_point
            low_point += neg_dist+increment
            high_point += neg_dist+increment
        print ""
        print "{}, {}, {}".format(low_point, high_point, increment)
        print "{0:.2f}\t{1:.2f}".format(i, np.arange(low_point, high_point, increment).mean())
        print "num points:", len(np.arange(low_point, high_point, increment))
    

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
        # modelled matrices subsected to match observation matrices
        self.modelled_sub = {}
        # corresponding (posibly) subsected observation matrices
        self.observations_sub = {}        
        # Dicts for raw observation matrices and liklihood scores
        self.observations = {}
        self.liklihood = {}
        # Dict of modelled vars & corresponding units
        self.model_var_units = {}
        # Flag to flip after execution and objects for output pipes to fill
        self.stdout, self.stderr, self.ran = None, None, False
        # Identifier for optimization variants
        self.variant_var, self.variant_val = "baseline", "baseline"

        self.name = name
        self.dir_path = dir_path
        make_dir(dir_path, verbose=False)
        self.glm_path = os.path.join(self.dir_path,"glm2.nml")
        
        config_pack = import_glm_config(dl_default=True, verbose=False)
        self.glm_config, self.default_config = config_pack[0], config_pack[1]
        self.default_order, self.block_order = config_pack[2], config_pack[3]
        self.baseline_configured=False
        self.aed_config = read_aed_nml()
        self.pathogen_config = read_cell_params('pathogen')
        self.phyto_config = read_cell_params('phyto')
        self.zoo_config = read_zoo_params()
        self.patho_params = load_aed_var_types('pathogen')
        self.phyto_params = load_aed_var_types('phyto')
        self.zoo_params = load_aed_var_types('zoo')
       
    def cleanup(self, verbose=True):
        shutil.rmtree(self.dir_path)
        
    def write_aed_files(self, aed_dict=None, phyto_dict=None, patho_dict=None, 
                        zoo_dict=None):
        """
        This function writes the values in the configuration dictionaries
        to the appropriate file names as specified in the `glm_config` object.
        If None is provided for a given dictionary, the default values are 
        written out.
        """
        aed_fn = self.glm_config['wq_setup']['wq_nml_file'][1:-1]
        self.aed_path = os.path.join(self.dir_path, aed_fn)
        phyto_fn = 'aed2_phyto_pars.nml'
        patho_fn = 'aed2_pathogen_pars.nml'
        zoo_fn = 'aed2_zoop_pars.nml'
        self.phyto_path = os.path.join(self.dir_path, phyto_fn)
        self.patho_path = os.path.join(self.dir_path, patho_fn) 
        self.zoo_path =  os.path.join(self.dir_path, zoo_fn)
        
        if not aed_dict:
            _ = write_aed2_config(self.aed_config, self.aed_path)
        else:
            _ = write_aed2_config(aed_dict, self.aed_path)
            
        if not phyto_dict:
            _ = write_aed2_eco_configs(self.phyto_config, self.phyto_path, 
                                       'phyto')
        else:
            _ = write_aed2_eco_configs(phyto_dict, self.phyto_path, 'phyto')
            
        if not patho_dict:
            _ = write_aed2_eco_configs(self.pathogen_config, self.patho_path, 
                                       'pathogen')
        else:
            _ = write_aed2_eco_configs(patho_dict, self.patho_path, 'pathogen')
            
        if not zoo_dict:
            _ = write_aed2_eco_configs(self.zoo_config, self.zoo_path, 'zoo')
        else:
            _ = write_aed2_eco_configs(zoo_dict, self.zoo_path, 'zoo')
        del _
        
    def write_glm_config(self, verbose=True):
        self.glm_config['output']['out_dir'] = "'"+self.dir_path+"'"
        self.glm_config['output']['out_fn'] = "'full_output'"
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
                    glm_handle.write("{0} = ".format(param))
                    for i, v in enumerate(value):
                        if (i+1) == len(value):
                            glm_handle.write("%r\n" % v)
                        else:
                            glm_handle.write("%r," % v)
                else:
                    glm_handle.write("{0} = {1}\n".format(param, value))
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
        """
        This fxn does temporal alignment of input datasets
        """
        self.firstDay = pd.to_datetime(self.glm_config['time']['start'][1:-1])
        self.lastDay = pd.to_datetime(self.glm_config['time']['stop'][1:-1])
        data_first_day = np.max([i.index[0] for i in data_pack])
        data_last_day = np.min([i.index[-1] for i in data_pack])
        
        if data_first_day > self.firstDay:
            sys.exit("Simulation begins before input data")
        elif data_last_day < self.lastDay:
            sys.exit("Prescribed simulation ends before input data")
        
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
            
        print "\tWriting met.csv for:"
        print self.met_df.index[0], "to", self.met_df.index[-1]
        self.met_df.to_csv(path_or_buf=self.metcsv_path, index_label='time')
        
    def create_flowcsvs(self, InPack, OutPack, verbose=False):
        """
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
        if verbose:
            print "Adding AED columns to a matrix of size:"
            print self.in_df.shape
        self.in_df_aed = add_aed_columns_to_input(self.in_df)
        if verbose:
            print "\tWriting inflow.csv"
        self.in_df_aed.to_csv(path_or_buf=self.incsv_path, index_label='time')
        if verbose:
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
                    print sub_block
                    print val_space
                    print dType
                    
                    if val_split[-1] == 'None':
                        val_split.pop()
                        
                    if dType == 'float':
                        val_split = [float(l) for l in val_split]
                    elif dType == 'list':
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
                    elif len(val_split) == 3 and dType == 'list':
                        old_val = self.glm_config[new_block][sub_block]
                        arange = np.arange(val_split[0], val_split[1]+val_split[2],
                                           val_split[2])
                        for idx, i in enumerate(arange):
                            new_val = np.array(old_val)*i
                            new_val = [int(np.floor(j*10000))/10000. for j in new_val]
                            sim_name = sub_block+"_"+str(idx+1)
                            block_key_option.append((new_block, sub_block, 
                                                     new_val, sim_name))
                            variants['grid'][sub_block].append(i)

                    elif len(val_split) == 3:
                        arange = np.arange(val_split[0], val_split[1]+val_split[2],
                                           val_split[2])
                        for idx, i in enumerate(arange):
                            j = int(np.floor(i*10000))
                            j = j/10000.
                            
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
        combination_max = 1
        for k in variants['grid'].keys():
            len_ = len(variants['grid'][k])
            combination_max *= len_
            variants['grid'][k] = np.array(variants['grid'][k])
            variants['prior'][k] = np.ones(len_)
            variants['likelihood'][k] = np.zeros(len_)
            variants['posterior'][k] = np.zeros(len_)

        print "\t%2.f possible models in this parameter space" % combination_max
            
        self.variants = variants

    def randomize(self, path, verbose=True):
        random_lake = copy.deepcopy(self)
        random_lake.ran = False
        random_lake.name = path
        random_lake.dir_path = os.path.join(self.dir_path, path)
        if not os.path.exists(random_lake.dir_path):
            os.mkdir(random_lake.dir_path)
        random_lake.glm_path = os.path.join(random_lake.dir_path, "glm2.nml")
        old_csvs = [self.metcsv_path, self.outcsv_path, self.incsv_path]
        new_csvs = [random_lake.metcsv_path, random_lake.outcsv_path, 
                    random_lake.incsv_path]
        for csv_old, csv_new in zip(old_csvs, new_csvs):
            csv_bn = os.path.basename(csv_old)
            csv_new = os.path.join(random_lake.dir_path, csv_bn)
            if not os.path.exists(csv_new):
                os.symlink(csv_old, csv_new)
        
        p_count = 0
        p_found = 0
        if self.variants:
            for param, vals in self.variants['grid'].items():
                p_count+=1
                p_flag = False
                randIdx = np.random.randint(len(vals))
                for block, sub_blocks in random_lake.glm_config.items():
                    if param in sub_blocks.keys():
                        p_found+=1
                        p_flag = True
                        if type(sub_blocks[param]) == list:
                            new_vals = [hard_round(i*vals[randIdx]) for i in sub_blocks[param]]
                        else:
                            new_vals = vals[randIdx]
                        print param, "changed from", sub_blocks[param], "to", new_vals
                        sub_blocks[param] = new_vals
                                                
            if p_flag == False:
                print param, "not found"
                
            print p_count, "external parameter spaces defined"
            print p_found, "of them matched to internal declarations"
            out_el = random_lake.glm_config['outflow']['outl_elvs']
            lake_base = random_lake.glm_config['morphometry']['H'][0]
            lake_crest = random_lake.glm_config['morphometry']['H'][-1]            
            if out_el < lake_base:
                random_lake.glm_config['outflow']['outl_elvs'] = lake_base
            elif out_el > lake_crest:
                random_lake.glm_config['outflow']['outl_elvs'] = lake_crest
                    
            random_lake.write_glm_config(False)
            
            aed_ps = [random_lake.aed_path, random_lake.phyto_path, 
                      random_lake.pathogen_path, random_lake.zoo_path]
                      
            aed_fs = [os.path.basename(i) for i in aed_ps]
            random_lake.write_aed_files(aed_fs[0],aed_fs[1],
                                            aed_fs[2],aed_fs[3])
        else:
            sys.exit("No variants detected")

            
        return random_lake
        
    def test_variant_configs(self, Lake_temps, copycsvs=False, verbose=True):
        
        if self.glm_config and self.variant_cases:
            variant_lakes = []
            lakes_n = len(self.variant_cases)
            for n_var, i in enumerate(self.variant_cases):
                variant_lake = copy.deepcopy(self)
                key1, key2, val, simname = i[0],i[1],i[2],i[3]
                listVars = ['the_temps', 'the_sals', 'A', 'H' ]
                    
                if key2 == "H":
                    old_elevs = self.glm_config['morphometry']['H'][1]
                    multiplier = val[1]/old_elevs
                    old_depths = self.glm_config['init_profiles']['the_depths']
                    new_depths = [i*multiplier for i in old_depths]
                    # change the lake depth to the new deepest elevation
                    variant_lake.glm_config['init_profiles']['the_depths'] = new_depths
                    variant_lake.glm_config['init_profiles']['lake_depth'] = np.array(val).max()
                    variant_lake.glm_config['outflow']['outl_elvs'] = np.array(val).max()
                    
                if key2 == "A":
                    old_areas = self.glm_config['morphometry']['A'][1]
                    multiplier = val[1]/old_areas                    
                    old_minV = self.glm_config['glm_setup']['min_layer_vol']
                    variant_lake.glm_config['glm_setup']['min_layer_vol'] = old_minV*multiplier
                
                variant_lake.glm_config[key1][key2] = val
                variant_lake.name = simname
                variant_lake.variant_var = key2
                if key2 in listVars:
                    old_list = self.glm_config[key1][key2][1]
                    variant_lake.variant_val = val[1]/old_list
                else:
                    variant_lake.variant_val = val
                    
                variant_lake.ran = False
                
                variant_lake.dir_path = os.path.join(self.dir_path, simname)
                make_dir(variant_lake.dir_path, verbose)
                variant_lake.glm_path = os.path.join(variant_lake.dir_path,
                                                     "glm2.nml")
                variant_lake.write_glm_config(verbose)
                aed_ps = [variant_lake.aed_path, variant_lake.phyto_path, 
                          variant_lake.pathogen_path, variant_lake.zoo_path]
                aed_fs = [os.path.basename(i) for i in aed_ps]
                variant_lake.write_aed_files(aed_fs[0],aed_fs[1],
                                                 aed_fs[2],aed_fs[3])
                
                if copycsvs==True:
                    if verbose==True:
                        print "Copying original data csvs to variant folders"
                    csv_files = [self.metcsv_path, self.outcsv_path,
                                 self.incsv_path]
                    obj_names = [os.path.basename(i) for i in csv_files]
                    obj_paths = [os.path.join(variant_lake.dir_path, i) 
                                                            for i in obj_names]
                    
                    for csV, objN in zip(csv_files, obj_paths):
                        if not os.path.exists(objN):
                            os.symlink(csV, objN)
                        else:
                            pass
                        
                    variant_lake.metcsv_path = obj_paths[0]
                    variant_lake.outcsv_path = obj_paths[1]
                    variant_lake.incsv_path = obj_paths[2]
                        

                else:
                    if verbose==True:
                        print "Original data csvs not moved into variant folders"
                    
                variant_lake.assign_observation_matrix(Lake_temps, 'temp', 'celsius')
                #Run them all & score them all
                score_list = run_test_variant((n_var, variant_lake, 
                                               lakes_n, self))
                del variant_lake
                variant_lakes.append(score_list)
        return variant_lakes
    
    def assign_observation_matrix(self, observations, var_type, units):
        if units:
            if units.lower() == 'mg/l':
                if var_type == 'OXY_oxy':
                    print "Converting DO obs. units"
                    mg_DO_to_uM = (1./0.0319988)
                    observations *= mg_DO_to_uM
                elif var_type == 'NIT_nit':
                    mg_NO3_to_uM = 16.1278
                    observations *= mg_NO3_to_uM
                    print "Converting NO3 obs. units"    
                units = 'mmol/m**3'
        
        if var_type in self.observations.keys():
            old_obs, unit_type = self.observations[var_type]
            if unit_type == units:
                new_obs = old_obs.append(observations)
                self.observations[var_type] = (new_obs, units)
            else:
                sys.exit('convert observations to compatible units')
        else:
            self.observations[var_type] = (observations, units)
            
            
        
        
    def score_variant(self, base, obs_type, lik_fxn='NSE'):
        # check if model was ran
        if self.ran == False:
            print "Scoring attempt on Bad Run Skipped"
            return np.nan
        
        # check if data was collected
        assert self.csv_dict

        #is top bottom orientation correct??
        # Depth DFs is a dict of all the depth & time indexed measurements
        # This is developed for 'temp' but can be expanded to check other 
        # locations for outputs to score
        if obs_type in self.depth_dfs.keys():
           
            
            # Now repeat the last 3 steps for all variant runs in variant_list
            # Pull out each tweaked parameter and the corresponding value
            sub_block_ = self.variant_var
            curr_val = self.variant_val
            # Pull out the value's position in the vector of all tested values
            modelled = self.depth_dfs[obs_type].resample('D').mean()
            self.z_index = self.depth_dfs['z'].resample('D').mean()
            
            try:
                observations, _ = self.observations[obs_type]  
            except NameError:
                sys.exit("No observations assigned")
                
            modelled_sub, observations_sub = time_depth_df_filter(modelled, 
                                                                  observations,
                                                                  self.z_index)
            
            self.modelled_sub[obs_type] = modelled_sub
            self.observations_sub[obs_type] = observations_sub

            err = error_metrics(observations_sub.values, modelled_sub)
            self.liklihood[obs_type] = err[lik_fxn]
                
            if sub_block_ != 'baseline':
                all_idxs = np.where(base.variants['grid'][sub_block_] == curr_val)
                idx = all_idxs[0]
                
                # Steps repeated as above, with the specefied error parameter
                # pulled out of the dictionary returned from `error_metrics`
                # and assigned the the original LakeModel Lake within the
                # variant dictionary in the 'likelihood' section at the same 
                # position as the input was in the 'grid' section
                base.variants['likelihood'][sub_block_][idx] = err[lik_fxn]
            

                # We then go back into the new vectors containing likelihood values
                for sb2 in base.variants['likelihood'].keys():
                    # We exclude parameters that are booleans or categorical
                    if len(base.variants['grid'][sb2]) > 4:
                        # Multiply the posterior by 1 because we assume a uniform
                        # distribution for the prior
                        raw_post=base.variants['prior'][sb2]*base.variants['likelihood'][sb2]
                    
                        # We integrate the raw posterior distribution along the 
                        # input grid of parameter values and divide the raw
                        # posterior by that value to get the real posterior
                        evidence=np.trapz(raw_post, base.variants['grid'][sb2])
                        base.variants['posterior'][sb2]=raw_post/evidence
        
        return err[lik_fxn]
    
def run_test_variant(variant_lake_tuple):
    i = variant_lake_tuple[0]
    l = variant_lake_tuple[1]
    mod_fam_size = variant_lake_tuple[2]
    baseline_case = variant_lake_tuple[3]
    
    # test disk size 
    statvfs = os.statvfs(os.getcwd())
    bytes_avail = statvfs.f_frsize * statvfs.f_bavail
    if bytes_avail < 104857600/2:
        sys.exit("Disk space less than 100 Mb, aborting")
    #print "~{} Mb of disk space remaining".format(bytes_avail/1024/1024)
    
    print "Running #", i+1, "out of", mod_fam_size
    #print "%r = %r" % (l.variant_var, l.variant_val)
    print l.glm_config['output']["out_fn"]
    ran_lake = run_model(l, verbose=False, forceRerun=True)
    print ran_lake.variant_val
    print ran_lake.variant_var
    data_lake = pull_output_nc(ran_lake, force=True)
    if data_lake.ran == True:
        try: 
            NSE = data_lake.score_variant(baseline_case, 'temp', lik_fxn='NSE')
        except ValueError:
            sys.exit()
        data_lake.cleanup(verbose=True)
        to_return = [data_lake.variant_var, data_lake.variant_val, NSE]
        del data_lake
        return to_return
    else:
        print "Variant %s failed and was removed" % data_lake.glm_config['glm_setup']['sim_name']
        return None

def time_depth_df_filter(baseline, observations, depths):
    #only extract the rows where we have corresponding dates
    time_matched_vals = baseline[baseline.index.isin(observations.index)]
    time_matched_heights = depths[depths.index.isin(observations.index)]
    observations_mod, _ = subsectbydate_2(observations, time_matched_vals)
    
    observed_dm = np.zeros(observations_mod.shape)
    obs_days, e_cols = observations_mod.shape
    for row in range(obs_days):
        this_row = time_matched_vals.iloc[row, :].dropna(axis=0)
        this_z = time_matched_heights.iloc[row, :].dropna(axis=0)
        
        elevations = range(e_cols)
        elevations.reverse()
        elevations = np.array(elevations) + 1
        depths_ = range(e_cols)
        
        for m,e in zip(depths_, elevations):
            if m != range(e_cols)[0]:
                up_elev_bool = (this_z > e-1)
                low_elev_bool = (this_z < e)
                sweet_spot_bool = up_elev_bool & low_elev_bool
                observed_dm[row, m] = this_row[sweet_spot_bool].mean()
            else:
                observed_dm[row, m] = this_row[(this_z > e)].mean()
                
    return observed_dm, observations_mod

    
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
                
def run_model(Lake, ex_loc=None, verbose=True, forceRerun=False):
    """This function accesses the locations of each glm configuration, 
    provided in Lake.glm_path, and runs the model using the executable
    specefied in `ex_loc`
    """
    report_path = os.path.join(Lake.dir_path, "Run_report.txt")
    print report_path
    unique_divider = "\n*\n"
    
    if os.path.exists(report_path) and forceRerun == False:
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
        if verbose and forceRerun == False:
            print "No prior run of model detected @ {}".format(report_path)
        else:
            print "Prior Run detected @ {}, and overwritten".format(report_path)
            
            
        if ex_loc == None:
            curr_dir = os.path.dirname(os.getcwd())
            ex_loc = 'GLM_Executables/glm.app/Contents/MacOS/glm'
            ex_loc = os.path.join(curr_dir, ex_loc)
            
        p = sp.Popen(ex_loc, cwd=Lake.dir_path, shell=True, stderr=sp.PIPE, 
                     stdout=sp.PIPE)
        #Lake.ret_code = p.wait()
        Lake.stdout, Lake.stderr = p.communicate()
        

        
        check1 = 'Simulation begins..' in Lake.stdout
        check2 = '100.00% of days complete' in Lake.stdout
        
        if check1 and check2:
            Lake.ran = True
            report_h = open(report_path, 'w')
            report_h.writelines("Run Report"+unique_divider)
            
            report_h.writelines(Lake.stdout+unique_divider)
            
            report_h.writelines(Lake.stderr+unique_divider)
            
            report_h.writelines(Lake.variant_var+unique_divider)
            
            report_h.writelines(str(Lake.variant_val)+unique_divider)
        
            if verbose:
                print "{} model run completed\n".format(Lake.name)
        else:
            print Lake.stdout
            print Lake.stderr
            print "Model Run Terminated Prematurely"
            Lake.ran = False
        
        return Lake
    


def pull_output_nc(Lake, force, verbose=False, plot=True, plot_var=[],
                   pickle=True):
    
    if Lake.ran == True:
        
        # This section calculates the expected time interval for saved time
        # steps in the netCDF archive
        o_mult = Lake.glm_config['output']['nsave']
        o_interval = Lake.glm_config['time']['dt']
        freq_S = str(int(o_mult*o_interval))+"S"        
        interval = pd.date_range(start=Lake.glm_config['time']['start'], 
                                 end=Lake.glm_config['time']['stop'], 
                                 freq=freq_S)
        Lake.nc_interval = interval
        if verbose == True:
            print "The time steps in the output netCDF are %s" % freq_S
        
        ## This section pulls out the expected path for the archive and also
        # checks to see if a pickled version is stored first
        
        o_fn = Lake.glm_config['output']["out_fn"][1:-1] + '.nc'
        o_path = Lake.glm_config['output']["out_dir"][1:-1]
        Lake.output_nc_path = os.path.join(o_path, o_fn)
        expected_pickle = Lake.output_nc_path[:-3]+".pickle"
        
        if os.path.exists(expected_pickle) and force == False:
            print "\tNC Pickle detected"
            print "\t\t", expected_pickle
            f = open(expected_pickle, 'rb')
            unpickled = cPickle.load(f)       
            f.close()
            Lake.depth_dfs, Lake.csv_dict = unpickled[0], unpickled[2]
            Lake.lake_df = unpickled[1]
            Lake.output_nc = None
        else:
            print Lake.output_nc_path
            output = Dataset(Lake.output_nc_path, "r")
            nc_vars = output.variables.keys()
            for ncv in nc_vars:
                try:
                    Lake.model_var_units[ncv] = output[ncv].units
                except AttributeError:
                    Lake.model_var_units[ncv] = None

            if 'time' in nc_vars:
                t_vector = output['time'][:]
                start_t = pd.to_datetime(output['time'].units[12:])
                d_vector = [start_t + DateOffset(hours=t) for t in t_vector]
            #print len(d_vector)
            #print len(Lake.nc_interval)
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
                    #print ncv
                    #print data_buff.mask.sum(), this_df.isnull().sum().sum()
                    if not (data_buff.mask.sum() == this_df.isnull().sum().sum()):
                        print "Warning: illegal variable size"
                    depth_dfs[ncv] = this_df
    
            Lake.depth_dfs = depth_dfs
            lake_df = pd.DataFrame(columns=lake_cols, index=d_vector, 
                                   data=lake_vecs)
            Lake.lake_df = lake_df
            
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
            
            if pickle == True:
                to_be_pickled = [depth_dfs, lake_df, csv_dict]
                Lake.output_pickle = expected_pickle
                f = open(Lake.output_pickle, 'wb')   # 'wb' instead 'w' for binary file
                cPickle.dump(to_be_pickled, f, -1)       # -1 specifies highest binary protocol
                f.close()
        
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
            Lake.output_nc = output
            output.close()
    else:
        print "Model not executed"
    return Lake

    
from collections import OrderedDict
from ChemDataVectors import LoadChemDFs


def add_aed_columns_to_input(inflow_orig):
    # Pull full inflow, meteorology, and outflow time series
    # Pull in surface & depth profile chem data
    depthGrad_df, surfaceM_df = LoadChemDFs()

    # Drop outflow & lake surface columns & trim column names
    extra_cols = [col for col in list(surfaceM_df.columns) 
                          if 'AberjonaRiver(Lower)' not in col]
    inflow_chem = surfaceM_df.drop(extra_cols, axis=1, inplace=False)
    new_columns = [i[:-21] for i in inflow_chem.columns]
    inflow_chem.columns = new_columns # concordance checked by mean
    
    ## Convert Units of Inflow vars
    
    # input O2 units are mg/L -> Multiply by (mmol / g of o2) * (1 L / 0.001 m^3)
    DO_multiplier = (3.12512*10**-2)*(1/0.001)
    depthGrad_df.DO = depthGrad_df.DO * DO_multiplier
    inflow_chem.DO = inflow_chem.DO * DO_multiplier
    
    # convert specific conductance to freshwater salinity 
    inflow_chem['Salinity'] = cond2sal(inflow_chem.SpecCond)
    
    # convert NH4 from mg/L to mmol / m^3
    NH4_multiplier = (1000/18.0385)
    inflow_chem.NH4 = inflow_chem.NH4 * NH4_multiplier
    TP_multiplier = 10.5294857189
    inflow_chem.TotalP = inflow_chem.TotalP*TP_multiplier
    NO_multiplier = 16.127757645
    inflow_chem.NO3NO2 = inflow_chem.NO3NO2*NO_multiplier
    
    # convert period index to datetime index centered on the 1st of each month
    new_index = []
    for period in inflow_chem.index:
        ide = pd.to_datetime(str(period)+'-01')
        new_index.append(ide)
    inflow_chem.index = new_index

    # Trim extratemporaneous time points from MRWA dataset
    new_index = pd.date_range((pd.Timestamp('20121231')+pd.DateOffset(months=-12)), 
                              periods=37, freq='MS')
    alt_chem = pd.DataFrame(index=new_index, columns=inflow_chem.columns,
                            data=inflow_chem.values)
    inflow_mod1, mrwa_chem = subsectbydate_2(inflow_orig, alt_chem)

    # Fill in the rest of the month using a second order spline
    mrwa_upsample = mrwa_chem.reindex(inflow_mod1.index)
    mrwa_up_splined = mrwa_upsample.interpolate(method='pchip')
    
    # TODO: Convert Units of Depth Profiles for grading purposes
    
    aed_inflow_params = OrderedDict()
    aed_inflow_params['OXY_oxy'] = None # 225 # None
    aed_inflow_params['SIL_rsi'] = 12.5
    aed_inflow_params['NIT_amm'] = None # 12.5 # None
    aed_inflow_params['NIT_nit'] = None # 27.6 # None
    aed_inflow_params['PHS_frp'] = None # 0.25 # None
    aed_inflow_params['OGM_don'] = 70.
    aed_inflow_params['OGM_pon'] = 5.
    aed_inflow_params['OGM_dop'] = 10.
    aed_inflow_params['OGM_pop'] = 10.
    aed_inflow_params['OGM_doc'] = 80.
    aed_inflow_params['OGM_poc'] = 700.
    aed_inflow_params['PHY_green'] = 40.
    aed_inflow_params['PHY_crypto'] = 40.
    aed_inflow_params['PHY_diatom'] = 40. 
    aed_inflow_params['CAR_dic'] = 1600.
    aed_inflow_params['CAR_pH'] = 7.5
    
    inflow_mod2 = inflow_mod1.copy()

    for name, value in aed_inflow_params.items():
        if value:
            inflow_mod2[name] = np.ones((inflow_mod2.shape[0],))*value
        elif name == 'OXY_oxy':
            inflow_mod2[name] = mrwa_up_splined['DO']
        elif name == 'NIT_amm':
            inflow_mod2[name] = mrwa_up_splined['NH4']
        elif name == 'NIT_nit':
            inflow_mod2[name] = mrwa_up_splined['NO3NO2']
        elif name == 'PHS_frp':
            inflow_mod2[name] = mrwa_up_splined['TotalP']
        else:
            sys.exit("Illegal Inflow Parameter")

    return inflow_mod2
    
    
def extractDepthDf(lake, var_string, time_depth_axis):
    var_df = lake.depth_dfs[var_string]
    print "Extracting", var_string
    var_meters, time_depth_axis = time_depth_df_filter(var_df, time_depth_axis, 
                                                       lake.z_index)
    var_df_clean = pd.DataFrame(index=var_df.index,
                                columns=range(var_meters.shape[1]),
                                data=var_meters)
    var_df_clean_2 = var_df_clean.interpolate(axis=1)
    return var_df_clean_2