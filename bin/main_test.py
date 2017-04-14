# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 21:02:22 2016

This is a script to autogenerate all the `glm.nml` files for all the cases
The first case will be sunny mild day with no variability


To run glm on linux, you need to get the file `l_fcompxe_2013.2.146_redist`
and later use the command when terminal starts:

    source /opt/intel/bin/compilervars.sh intel64

@author: Keith Arora-Williams
"""
from __future__ import division
import time
import os
import cPickle
import LakeModel
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import copy, sys
from collections import OrderedDict

plotting = False
plt.style.use('fivethirtyeight')
print "\nGLM_WRAPPER"
print time.ctime(time.time())
print ""

# these are code names for each trial
cases = ["testcase_5"]

# this is the directory you want to work in
mainFolder = os.path.dirname(os.getcwd())

# GLM starts from the command line and immediately searches in the current 
# directory for `glm2.nml`. So we need a directory for each case and a 
# directory to hold them all which will be:

superDir = os.path.join(mainFolder, 'glm_case_folders')
LakeModel.make_dir(superDir, verbose=False)

print "Adding all data file paths to memory"

met_fn = "BostonLoganAirportWeather.csv"
ceres_fn = "CERES_SSF_XTRK-MODIS_Edition3A_Subset_2010010102-2014010116.nc"
ceres_fn2 = "CERES_SYN1deg-Day_200508-201406.nc"

waterDataFiles = { 'in' : "Aberjona_Discharge.txt",
                   'fresh' : "Fresh_Discharge_WaterT_AirT_Cond_Precip.txt",
                   'Ho' : "Hobbs_Out_Discharge_WaterT_AirT_Cond_Precip.txt",
                   'Hi1' : "Hobbs_in1_discharge_waterT_cond.txt",
                   'Hi2' : "Hobbs_in2_discharge_waterT_cond.txt",
                   'Hi3' : "Hobbs_in3_discharge_waterT_cond.txt",
                   'Hi4' : "Hobbs_in4_discharge_waterT_cond.txt",
                   'ML' : "MysticLake_TimeDepthTemp.csv",
                   'out_d' : "MysticRiver_Discharge.txt",
                   'out2' : "MysticRiver2_Discharge.txt",
                   'out_t' : "MysticRiver_GH_WaterTemp_15min.txt",
                   'ground' : "Wilmington_groundwater.txt"}

met_path = os.path.join(mainFolder, 'weatherData', met_fn)
ceres_ssf = os.path.join(mainFolder, 'weatherData', ceres_fn)
ceres_SYN1 = os.path.join(mainFolder, 'weatherData', ceres_fn2)

waterDataPaths = {i:os.path.join(mainFolder, 'waterdata', waterDataFiles[i]) 
                  for i in waterDataFiles.keys()}

newDirs = [os.path.join(superDir, x) for x in cases]
_ = [LakeModel.make_dir(i, verbose=False) for i in newDirs]

print "Instantiating a lake"
test = LakeModel.Lake(cases[0], newDirs[0])

test.write_glm_config(verbose=False)

BOS_weather = LakeModel.GHCN_weather_data(met_path)
BOS_weather.clean_columns()

CERES_SSF = LakeModel.CERES_nc(ceres_ssf, "3")
CERES_data = LakeModel.CERES_nc(ceres_SYN1, "4")

inflow = LakeModel.USGS_water_data(waterDataPaths['in'])
mystic_d = LakeModel.USGS_water_data(waterDataPaths['out_d'])
mystic_t = LakeModel.USGS_water_data(waterDataPaths['out_t'])
Hobbs_o = LakeModel.USGS_water_data(waterDataPaths['Ho'])
Hobbs_i1 = LakeModel.USGS_water_data(waterDataPaths['Hi1'])
Hobbs_i2 = LakeModel.USGS_water_data(waterDataPaths['Hi2'])
Hobbs_i3 = LakeModel.USGS_water_data(waterDataPaths['Hi3'])
Hobbs_i4 = LakeModel.USGS_water_data(waterDataPaths['Hi4'])
GroundW = LakeModel.USGS_water_data(waterDataPaths['ground'])

waterObjects = [Hobbs_o, Hobbs_i1, Hobbs_i2, Hobbs_i3, Hobbs_i4, inflow, 
                mystic_d, mystic_t, GroundW]

for h in waterObjects:
    h.preprocess()
    h.read()
    h.removeQuals(verbose=False)

GroundW. df['GroundwaterLevel'] = (GroundW.df[GroundW.df.columns[0]]-96.6)*-1

Lake_temps = pd.read_csv(waterDataPaths['ML'], index_col=0, 
                         parse_dates=['datetime'])
                         
# Reindex columns with "elevations" and not "depths"
depth_rev = range(len(Lake_temps.columns))
depth_rev.reverse()
Lake_temps.columns = depth_rev

CERES_data.df1.cldarea_total_daily.fillna(method='pad', inplace=True)

SSF_D = CERES_SSF.ceres_df.resample('D').mean()
SSF_Filled = SSF_D.interpolate()
filledagain = SSF_Filled.fillna(method='bfill')

#Pack up met.csv variables
cloudFrac = CERES_data.df1['cldarea_total_daily']
cloudFrac.name = "Clouds"
longWaveIn = CERES_data.df1['sfc_comp_lw-down_all_daily']
longWaveIn.name = "LongWave"
shortWaveNet = CERES_data.df1['sfc_comp_sw-down_all_daily'] - CERES_data.df1['sfc_comp_sw-up_all_daily']
shortWaveNet.name = "ShortWave"
airTemp = (BOS_weather.df.t_max + BOS_weather.df.t_min) / 2.0
airTemp.name = "AirTemp"
relHumidity = filledagain.humidity
relHumidity.name = "RelHum"
windSpeed = BOS_weather.df.wind_speed
windSpeed.name = "WindSpeed"
precipR = BOS_weather.df.precip
precipR.name = "Rain"
precipS = BOS_weather.df.snow_fall
precipS.name = "Snow"

# Pack up variables for flow csvs
    
In_Discharge = LakeModel.CubicftPerS_to_MegalitersPerDay(inflow.df['Discharge, cubic feet per second (Mean)'])
In_Discharge.name = "INFLOW"

Out_Discharge = copy.deepcopy(In_Discharge)*0.8
Out_Discharge.name = "OUTFLOW"

In_Temp = Hobbs_o.df['Temperature, water, degrees Celsius (Mean)'].interpolate()
In_Temp.name = "TEMP"
In_Cond = Hobbs_o.df[Hobbs_o.df.columns[6]].interpolate() 
In_Salt = LakeModel.cond2sal(In_Cond)
In_Salt.name = "SALT"


data_pack = [shortWaveNet, longWaveIn, airTemp, relHumidity, windSpeed, 
             precipR, precipS, cloudFrac, Out_Discharge, In_Discharge, 
             In_Temp, In_Salt]
             
test.temporal_clipping(data_pack)

metPack = ['ShortWave','LongWave','AirTemp','RelHum','WindSpeed','Rain','Snow',
           'Clouds']
inPack = ['INFLOW','TEMP','SALT']
outPack = ['OUTFLOW']

test.create_metcsv(metPack)
test.create_flowcsvs(inPack, outPack)
test.write_aed_files(None, None, None, None)

## This version is the best version
test = LakeModel.run_model(test, verbose=False, forceRerun=True)
test = LakeModel.pull_output_nc(test, force=True)

# Pull Observations 
chem_dir = os.path.join(mainFolder, 'ChemData')
obs_pack = LakeModel.readInElisasAndPreheimData(chem_dir)
do_elisas, ph_elisas, do_preheim, no3_preheim, temp_preheim = obs_pack

# Assign Observations
test.assign_observation_matrix(temp_preheim, 'temp', 'celsius')
test.assign_observation_matrix(do_elisas, 'OXY_oxy', 'mg/L')
test.assign_observation_matrix(ph_elisas, 'CAR_pH', None)
test.assign_observation_matrix(do_preheim, 'OXY_oxy', 'mg/L')
test.assign_observation_matrix(no3_preheim, 'NIT_nit', 'mg/L')

# Score model against available observations
temp_NSE = test.score_variant(test, 'temp', lik_fxn='NSE')
print "Temp error score: {}".format(temp_NSE)
oxy_NSE = test.score_variant(test, 'OXY_oxy', lik_fxn='NSE')
print "Oxy error score: {}".format(oxy_NSE)
no3_NSE = test.score_variant(test, 'NIT_nit', lik_fxn='NSE')
print "NO3 error score: {}".format(no3_NSE)
#ph_NSE = test.score_variants(test, 'CAR_pH', lik_fxn='NSE')
#print "PH error score: {}".format(ph_NSE)

LakeModel.plotLakeandError(test, test.observations_sub, 'OXY_oxy', 9)
LakeModel.plotLakeProfiles(test, 10)

# Make skeleton dataframe with expected time/depth shape
# Pull out depth dfs, convert layers to depths, and prepare plot
time_depth_axis = pd.DataFrame(index=test.depth_dfs['NIT_amm'].index,
                               columns=Lake_temps.columns)

fig1_od, fig2_od = OrderedDict(), OrderedDict()
fig3_od, fig4_od = OrderedDict(), OrderedDict()
fig5_od, fig6_od = OrderedDict(), OrderedDict()
fig7_od, fig8_od = OrderedDict(), OrderedDict()

amm_df = LakeModel.extractDepthDf(test, 'NIT_amm', time_depth_axis)
fig1_od['Ammonia'] = amm_df

nit_df = LakeModel.extractDepthDf(test, 'NIT_nit', time_depth_axis)
fig1_od['Nitrate'] = nit_df

meth_df = LakeModel.extractDepthDf(test, 'CAR_ch4', time_depth_axis)
fig2_od['Methane'] = meth_df

meth_ox_df = LakeModel.extractDepthDf(test, 'CAR_ch4ox', time_depth_axis)
fig2_od['MethaneOxidation'] = meth_ox_df

nitrif_df = LakeModel.extractDepthDf(test,'NIT_nitrif', time_depth_axis)
fig3_od['Nitrification'] = nitrif_df

denit_df = LakeModel.extractDepthDf(test, 'NIT_denit', time_depth_axis)
fig3_od['Denitrification'] = denit_df

dic_df = LakeModel.extractDepthDf(test, 'CAR_dic', time_depth_axis)
fig4_od['Dissolved Inorganic Carbon'] = dic_df

oxy_df = LakeModel.extractDepthDf(test, 'OXY_oxy', time_depth_axis)
fig4_od['Oxygen'] = oxy_df

phs_df = LakeModel.extractDepthDf(test, 'PHS_frp', time_depth_axis)
fig5_od['Phosphorus'] = phs_df

phy_green_df = LakeModel.extractDepthDf(test, 'PHY_green', time_depth_axis)
fig5_od['Model Green Algae'] = phy_green_df

pH_df = LakeModel.extractDepthDf(test, 'CAR_pH', time_depth_axis)
fig6_od['pH'] = pH_df.round(decimals=3)

doc_df = LakeModel.extractDepthDf(test, 'OGM_doc', time_depth_axis)
fig6_od['Dissolved Organic Carbon'] = doc_df

aed_df_pack = [amm_df, nit_df, nitrif_df, denit_df]
aed_df_fns = ['AmmoniaVdepthVtime.csv',
              'NitrateVdepthVtime.csv',
              'NitrificationVdepthVtime.csv',
              'DenitrificationVdepthVtime.csv']
              
for aed_df, aed_fn in zip(aed_df_pack, aed_df_fns):
    aed_df.to_csv(aed_fn)


temp_df = LakeModel.extractDepthDf(test, 'temp', time_depth_axis)
sal_df = LakeModel.extractDepthDf(test, 'salt', time_depth_axis)
radiation_df = LakeModel.extractDepthDf(test, 'rad', time_depth_axis)
height_df = LakeModel.extractDepthDf(test, 'z', time_depth_axis)

fig7_od["Temperature"] = temp_df
fig7_od["Salinity"] = sal_df
fig8_od["Solar Radiation"] = radiation_df
fig8_od["Layer Height"] = height_df

def plot2profiles(df_lab_dict, fignum):
    fig = plt.figure(fignum, figsize=(18,9))
    labs, dfs = zip(*df_lab_dict.items())
    for num, lab, df_ in zip([1,2], labs, dfs):
        ax = fig.add_subplot(2,1,num)
        ax.set_title(lab, fontsize=14)
        cax = ax.imshow(df_.T, interpolation='nearest', aspect='auto')
        divider = make_axes_locatable(ax)
        cax1 = divider.append_axes("right", size="5%", pad="3%")
        plt.colorbar(cax, cax=cax1)
        
plot2profiles(fig1_od, 1)
plot2profiles(fig2_od, 2)
plot2profiles(fig3_od, 3)
plot2profiles(fig4_od, 4)
plot2profiles(fig5_od, 5)
plot2profiles(fig6_od, 6)
plot2profiles(fig7_od, 7)
plot2profiles(fig8_od, 8)


sys.exit()

#This is the potential model space 
test.read_variants('../test_files/optimize_these.txt', verbose=False)
#Here we want to create a new lake at a random starting point

do_edging_optimiztion = False
do_random_optimization = True

if do_edging_optimiztion:
    variant_lakes = test.test_variant_configs(Lake_temps, copycsvs=True, 
                                          verbose=False)
    variants_ = np.array(variant_lakes)
    y = variants_[:, 2]
    print "failed runs\n", len(variants_[np.isnan(y.astype(float))])
    if len(variants_[np.isnan(y.astype(float))]) > 0:
        variants_real = variants_[~np.isnan(y.astype(float))]
        y = variants_real[:, 2]
        print "best score: ", y.astype(float).max()
        print "parameter value: ", variants_[np.argmax(y.astype(float)), :]
    else:
        print "best score: ", y.astype(float).max()
        print "parameter value: ", variants_[np.argmax(y.astype(float)), :]

if do_random_optimization:
    # Find the best improvement
    start_time = time.time()
    rand_case_name = "randomize_5"
    total_trials = 25250-2500
    max_pickle_size = 500
    repeats = int(np.ceil(total_trials/ max_pickle_size))
    trial_counter = 0
    for mults in range(repeats):
        bad_lakes = []
        single_run = time.time()
        
        if (total_trials - trial_counter) < max_pickle_size:
            these_boots = total_trials - trial_counter
        else:
            these_boots = max_pickle_size
        errors = np.zeros(these_boots)
        bootstraps = np.arange(these_boots)
        this_err = np.nan
        random_lake = test.randomize(rand_case_name, True)
        bstps_results = {'config':{},
                         'error':{}}
        
        for rep in bootstraps:
            trial_counter+=1
            # retry if run fails or error is nan
            while random_lake.ran == False or np.isnan(this_err):
                random_lake = LakeModel.run_model(random_lake, verbose=True, 
                                                  forceRerun=True)        
                random_lake = LakeModel.pull_output_nc(random_lake, force=True)
                this_err = random_lake.score_variant(test, 'temp', lik_fxn='NSE')
                if random_lake.ran == False or np.isnan(this_err):
                    if len(bad_lakes) < 10000:
                        bad_lakes.append(copy.deepcopy(random_lake.glm_config))
                    del random_lake
                    random_lake = test.randomize(rand_case_name, True)
                    
            # Save the config and the error
            bstps_results['config'][rep] = copy.deepcopy(random_lake.glm_config)
            bstps_results['error'][rep] = this_err
            del random_lake
            random_lake = test.randomize(rand_case_name, True)
        
        
        print("--- %s, %s seconds ---" % ((time.time() - start_time),
                                          (time.time() - single_run)))
        pickleName = rand_case_name+"_res_{}.pickle".format(mults+1)
        to_be_pickled = [bstps_results, bad_lakes]
        results_pickle = os.path.join(random_lake.dir_path, pickleName)
        with open(results_pickle, 'wb') as f:
            cPickle.dump(to_be_pickled, f, -1)
            
    sys.exit()

# the lake with best likelihood is singled out

LakeModel.print_minimums(test)