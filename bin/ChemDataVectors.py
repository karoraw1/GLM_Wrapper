#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 14:48:42 2017



@author: login
"""

import os, sys
import numpy as np
import pandas as pd

def readChemData(chem_path, units, ftype, plotbool=False, verbose=False):
    if verbose:
        print os.path.basename(chem_path)
        print units
    
    if ftype == 'depth_profile':
        site_chem_spec_df = pd.read_csv(chem_path, index_col=0, 
                                        parse_dates=True, 
                                        infer_datetime_format=True)
        new_idx = []
        for i in site_chem_spec_df.index:
            new_idx.append(pd.Period(i, 'M'))
        site_chem_spec_df.index = new_idx
        
    elif ftype == 'surface_measure':
        chem_spec_csv = pd.read_csv(chem_path, sep=",", index_col=0)
        if verbose:
            print "Null Values: {}".format(chem_spec_csv.isnull().sum().sum())
            print "Database Shape: {}".format(chem_spec_csv.shape)
        new_cols = [pd.Period(i) for i in chem_spec_csv.columns]
        chem_spec_csv.columns = new_cols
        site_chem_spec_df = chem_spec_csv.T.interpolate()
    else:
        sys.exit("invalid ftype")
        
    if plotbool == True:
        site_chem_spec_df.plot()
    return site_chem_spec_df
    
def LoadChemDFs(verbose=False):
    if verbose:
        print "Reading in all Chemical Data"
    
    chem_path = '/Users/login/Documents/GLM_Wrapper/OTU_Time_Series/ChemData'
    chem_csvs = sorted([os.path.join(chem_path, csv) for csv in os.listdir(chem_path) 
                                              if csv[-4:] == '.csv'])
    
    
    depth_grads = []
    surface_meas = []
    unit_lookup = {'Ammonia.csv': 'mg/L',
                   'DO_Depths.csv': 'mg/L',
                   'DissolvedOxygen.csv': 'mg/L',
                   'InitialValues.csv': 'mg/L',
                   'NitrateNitrite.csv': 'mg/L',
                   'Nitrate_Depths.csv': 'mg/L',
                   'SpecConductance.csv': 'uS/cm',
                   'Sulfate_Depths.csv': 'mg/L',
                   'TSS.csv': 'mg/L',
                   'TotalPhosphorus.csv': 'mg/L',
                   'EColi.csv':'MPN/100mL'}
                   
    
    for f in chem_csvs:
        if "Depths" in f:
            ftype = 'depth_profile'
        else:
            ftype = 'surface_measure'
        any_csv = readChemData(f, unit_lookup[os.path.basename(f)], ftype, False)
        if any_csv.shape[0] > any_csv.shape[1]:
            surface_meas.append(any_csv)
        elif any_csv.shape[0] < any_csv.shape[1]:
            depth_grads.append(any_csv)
        
    our_sites = ['AberjonaRiver(Lower)','MysticRiver(Upper)', 'UpperMysticLake']
                 
    
    column_prefixes = ['NH4', 'DO', 'EColi', 'NO3NO2', 'SpecCond', 'TSS', 'TotalP']
    for csv, pfx in zip(surface_meas, column_prefixes):
        new_cols = [pfx+"_"+i for i in csv.columns]
        csv.columns = new_cols
    
    surfaceM_df = pd.concat(surface_meas, axis=1, verify_integrity=True)
    
    def flatten_dg_dfs(df, varName):
        new_dates = []
        new_depths = []
        flat_vals = []
        for da in df.index:
            for de in df.columns:
                new_dates.append(da)
                new_depths.append(de)
                flat_vals.append(df.ix[da,de])
        new_cols = ['date', 'depth', varName]
        col_vals = [new_dates, new_depths, flat_vals]
        df_dict = {i:j for i, j in zip(new_cols, col_vals)}
        new_df = pd.DataFrame.from_dict(df_dict)
        return new_df
    
    
    varNames = ['DO', "Nitrate", "Sulfate"]
    dg_dfs = [flatten_dg_dfs(df, varName) for df, varName in zip(depth_grads, varNames)]
    
    depthGrad_df = pd.concat(dg_dfs, axis=1).T.drop_duplicates().T
    
    return (depthGrad_df, surfaceM_df)
    













