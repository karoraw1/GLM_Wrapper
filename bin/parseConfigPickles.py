# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 12:28:36 2016

@author: login
"""

import cPickle, os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
import pandas as pd

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


optimized = ["Kw", "min_layer_vol", "min_layer_thick", "max_layer_thick", 
             "coef_mix_conv", "coef_wind_stir", "coef_mix_shear", 
             "coef_mix_turb", "coef_mix_KH", "coef_mix_hyp", "bsn_len", 
             "bsn_wid", "H", "the_sals", "bsn_len_outl", "bsn_wid_outl", 
             "outflow_factor", "coef_inf_entrain", "strmbd_drag", 
             "strmbd_slope", "strm_hf_angle", "rain_threshold", "runoff_coef",
             "wind_factor", "rain_factor", "at_factor",  "rh_factor", 
             "sw_factor", "lw_factor", "cd", "ce", "ch"]

prefix_list = ['run3_results']

dir_list = ['randomize_2']
result, bad_lakes = parseConfigPickles(dir_list, prefix_list)
errorList = np.array([n[1] for n in result.values()])
print "{} individual simulations detected".format(len(errorList))
configList = [o[0] for o in result.values()]
plt.figure(1)
plt.xlabel("Nash Sutcliffe Efficiency")
plt.ylabel('No. simulations (n)')
plt.hist(errorList, bins=80, color='g')

sortedErr = sorted(errorList)
pct_denom = len(sortedErr)
idx_cutoffs = [int(pct_denom*0.95), int(pct_denom*0.99)]
err_cutoffs = [sortedErr[i] for i in idx_cutoffs]

for i in err_cutoffs:
    denom = len(errorList)
    numer = float((errorList > i).sum())
    frac_above = (int((numer/denom)*1000.)/1000.)*100. 
    print "{}% of runs had NSE above {}".format(frac_above, i)

error_value_rels = {}
hsy_parsed = {}
for _, p in configList[0].items():
    for q in p.keys():
        hsy_parsed[q] = {'behavioural':[], 'non-behavioural':[], 'breakers':[]}
        error_value_rels[q] = {'value':[], 'error':[]}

bl_A = np.array([85110.3, 163322.5, 222719.2, 283600.9, 372408.3, 
                        436785.8, 506855.8, 577283.3, 616056.1])
bl_H = np.array([0.132, 4.055, 8.112, 12.169, 16.227, 20.284, 24.339, 
                 28.396, 32.454])
bl_sals = np.array([197.6, 395.16, 592.74, 790.32, 988., 1185.48])
bl_temps = np.array([5.04, 5.04, 5.04, 5.04, 5.04, 5.04])

listVars = ['the_temps', 'the_sals', 'A', 'H' ]
listVals = [bl_temps, bl_sals, bl_A, bl_H]

baselineLists = {t:u for t, u in zip(listVars, listVals)}
sig_params = ['A', 'max_layer_thick', 'the_sals', 'H']
                     
run_id = range(len(errorList))
x = np.zeros((len(errorList), len(optimized)))
y = np.zeros((len(errorList),))
 
for config, error, row_l in zip(configList, errorList, run_id):
    #load config & error
    for r in config.keys():
        param_block = config[r]
        for param, value in param_block.items():
            if param in listVars:
                multArray = (np.array(value)/baselineLists[param]).round(decimals=1)
                this_val = round(np.median(multArray),1)
            else:
                this_val = value
            
            if error > err_cutoffs[0]:
                hsy_parsed[param]['behavioural'].append(this_val)
                if param in sig_params:
                    error_value_rels[param]['error'].append(error)
                    error_value_rels[param]['value'].append(this_val)
            elif error <= err_cutoffs[0]:
                hsy_parsed[param]['non-behavioural'].append(this_val)
        
    for idx, param2 in enumerate(optimized):
        
        for r in config.keys():
            param_block = config[r]    
            if param2 in param_block.keys():
                new_val = param_block[param2]
                if param2 in listVars:
                    multArray = (np.array(new_val)/baselineLists[param2]).round(decimals=1)
                    this_val = round(np.median(multArray),1)
                else:
                    this_val = new_val
                
                x[row_l, idx] = this_val
                y[row_l] = error
            else:
                pass
            

from sklearn import preprocessing
from sklearn import linear_model
x_scaled = preprocessing.scale(x)
clf = linear_model.Ridge()
clf.fit(x_scaled, y*100)
coeffs = pd.DataFrame(index = optimized, data = clf.coef_)
coeffs.sort_values(0, inplace=True, ascending=False)

for bL in bad_lakes:
    for r in bL.keys():
        param_block = bL[r]
        for param, value in param_block.items():
            if param in listVars:
                multArray = (np.array(value)/baselineLists[param]).round(decimals=1)
                this_val = round(np.median(multArray),1)
            else:
                this_val = value

            hsy_parsed[param]['breakers'].append(this_val)

def HSY_Sensitivity(group1, group2):
    full_set = group1 + group2
    x = np.linspace(min(full_set), max(full_set))
    ecdf_1 = sm.distributions.ECDF(group1)
    ecdf_2 = sm.distributions.ECDF(group2)
    cumdist_1 = ecdf_1(x)
    cumdist_2 = ecdf_2(x)
    ks_stat, p_val = stats.ks_2samp(cumdist_1, cumdist_2)
    return x, cumdist_2, cumdist_1, p_val


breaker_ps = {}
p_vals = {}
sig_distributions = {}
breakers_sigs = {}
for param in hsy_parsed.keys():
    if param in optimized:
        behavioural = hsy_parsed[param]['behavioural']
        non_behavioural = hsy_parsed[param]['non-behavioural']
        breakers = hsy_parsed[param]['breakers']
        p_vals[param] = HSY_Sensitivity(behavioural, non_behavioural)
        breaker_ps[param]= HSY_Sensitivity(non_behavioural+non_behavioural, breakers)
        
        if p_vals[param][3] < 0.05:
            sig_distributions[param] = (p_vals[param][3],
                                        p_vals[param][0],
                                        p_vals[param][1],
                                        p_vals[param][2])
                                        
        if breaker_ps[param][3] < 0.05:
            breakers_sigs[param] = (breaker_ps[param][3],
                                    breaker_ps[param][0],
                                    breaker_ps[param][1],
                                    breaker_ps[param][2])
        


"""
plt.figure(2)
plt.ylabel("Probability")
ax1 = plt.subplot(211)
plt.title("Depth Profile Scalar")
plt.plot(sig_distributions['H'][1], sig_distributions['H'][2],
         label="non-behavioural")
plt.plot(sig_distributions['H'][1], sig_distributions['H'][3],
         label="behavioural")
plt.legend(loc='upper left')

# share x and y
ax2 = plt.subplot(212)
plt.title("Salinitiy Profile Scalar")
plt.plot(sig_distributions['the_sals'][1], sig_distributions['the_sals'][2])
plt.plot(sig_distributions['the_sals'][1], sig_distributions['the_sals'][3])
plt.show()
sys.exit()
plt.figure(3, figsize=(9,9))
ax1 = plt.subplot(211)
plt.plot(breakers_sigs['max_layer_thick'][1], breakers_sigs['max_layer_thick'][2],
         label="Invalid")
plt.plot(breakers_sigs['max_layer_thick'][1], breakers_sigs['max_layer_thick'][3],
         label="Valid")
plt.legend(loc='upper left')
plt.title("Maximum Layer Thickness")
# share x and y
ax2 = plt.subplot(212)
plt.plot(breakers_sigs['A'][1], breakers_sigs['A'][2])
plt.plot(breakers_sigs['A'][1], breakers_sigs['A'][3])
plt.xlim([0.5, 1.6])
plt.ylabel("Probability")
plt.title("Area Profile Scalar")
ax3 = plt.subplot(313)
plt.plot(breakers_sigs['H'][1], breakers_sigs['H'][2])
plt.plot(breakers_sigs['H'][1], breakers_sigs['H'][3])
plt.xlim([1.0, 1.7])
plt.title("Depth Profile Scalar")
plt.xlabel("Parameter Value")

plt.figure(4)
ax1 = plt.subplot(122)
plt.scatter(error_value_rels['the_sals']['error'], 
            error_value_rels['the_sals']['value'] )
hline = plt.ylim()
plt.plot(np.ones(2)*0.858, hline, c='r')
plt.title("Salinity Profile Scalar")
plt.xlabel("Nash-Sutcliffe Efficiency")
plt.ylim(hline)            
ax2 = plt.subplot(121)
plt.scatter(error_value_rels['H']['error'], 
            error_value_rels['H']['value'] )
hline = plt.ylim()
plt.plot(np.ones(2)*0.858, hline, c='r')
plt.title("Depth Profile Scalar")
plt.ylim(hline)
plt.xlabel("Nash-Sutcliffe Efficiency")
plt.ylabel("Parameter Value")
"""