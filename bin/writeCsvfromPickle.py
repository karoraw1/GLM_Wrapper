#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 19:44:48 2017

@author: login
"""
import os, sys, cPickle

def parseConfigPickles(dir_list, prefix_list, verbose=True):
    rootPath = os.path.dirname(os.getcwd())
    caseBase = os.path.join(rootPath, 'glm_case_folders')

    if os.path.exists(caseBase) and verbose:
        print "\nBase Directory Detected"
    else:
        sys.exit("\nBase Directory Not Detected")
    
    caseDirs = [os.path.join(caseBase, i) for i in dir_list]
    
    cases = []
    for i in caseDirs:            
        if os.path.exists(i) and verbose:
            print "\nCase Directory {} Detected".format(i)
        else:
            print("\nCase Directory {} Not Detected".format(i))
        counter = 0    
        for j in os.listdir(i):
            for k in prefix_list:
                if k in j and 'pickle' in j: 
                    counter+=1
                    pickle_path = os.path.join(i, j)
                    cases.append(pickle_path)
                    if verbose:
                        
                        print "\t{}, ({}) pickle detected".format(j, os.path.getsize(pickle_path))
    
    uniq_cases = list(set(cases))
    if verbose:
        print "\n{} total cases detected".format(len(cases))
        print "Reduced to {} after dereplicaion".format(len(uniq_cases))
    
    compiled_runs = {}
    all_bad_lakes = []
    master_counter = 0
    for l in uniq_cases:
        f = open(l, 'rb')
        bstps_results, bad_lakes = cPickle.load(f)
        all_bad_lakes+= bad_lakes
        f.close()
        for m in bstps_results['config'].keys():
            master_counter += 1 
            compiled_runs[master_counter] = (bstps_results['config'][m],
                                             bstps_results['error'][m])
    return compiled_runs, all_bad_lakes
    


dir_list = ['/Users/login/Documents/GLM_Wrapper/glm_case_folders/testcase_5/interrupted_pickles']
prefix_list = ['randomize_5_res_']

compiled_runs, all_bad_lakes = parseConfigPickles(dir_list, prefix_list, verbose=True)    

opted_cols = ["coef_wind_stir", "coef_mix_conv", "coef_mix_shear", "coef_mix_hyp", 
              "coef_mix_turb", "coef_mix_KH", "Kw", "bsn_len", "bsn_wid", 
              "the_sals", "bsn_len_outl", "bsn_wid_outl", "coef_inf_entrain",
              "strmbd_drag", "strmbd_slope", "strm_hf_angle", "inflow_factor",
              "rain_threshold", "wind_factor", "rain_factor", "at_factor",
              "sw_factor", "lw_factor", "cd", "ce", "ch", "the_temps",
              "rh_factor"]
        
trials = compiled_runs.keys()

def pull_out_column(trial_no):
    var_values = []
    trial_data = compiled_runs[trial_no]
    var_values.append(trial_data[1])
    remaining_data = trial_data[0]
    heads = sorted(remaining_data.keys())
    for header in heads:
        subsection = remaining_data[header]
        subheads = sorted(subsection.keys())
        for subhead in subheads:
            if subhead in opted_cols:
                if subhead == 'the_temps' or subhead == 'the_sals':
                    var_values.append(subsection[subhead][0])
                else:
                    var_values.append(subsection[subhead])
    return var_values

columns = map(pull_out_column, trials)

var_headers = []
trial_data = compiled_runs[1]
var_headers.append("\"NSE\"")
remaining_data = trial_data[0]
heads = sorted(remaining_data.keys())
for header in heads:
    subsection = remaining_data[header]
    subheads = sorted(subsection.keys())
    for subhead in subheads:
        if subhead in opted_cols:
            var_headers.append('"'+subhead+'"')
            
import pandas as pd

opt_df = pd.DataFrame(index=trials,
                      columns=var_headers,
                      data = columns)

opt_df.to_csv(path_or_buf="/Users/login/Documents/DataAnalytics/Homework/HW1_Data1.csv",
              index_label=False, index=False, quotechar="\"", quoting=0)
            
            
            
            
            