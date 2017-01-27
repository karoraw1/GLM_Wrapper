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
import pandas as pd
import copy, sys

plotting = False
plt.style.use('fivethirtyeight')
print "\nGLM_WRAPPER"
print time.ctime(time.time())
print ""

# these are code names for each trial
cases = ["testcase_4"]

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
In_Salt = Hobbs_o.df[Hobbs_o.df.columns[6]].interpolate()
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

## This version is the best version
test = LakeModel.run_model(test, verbose=False, forceRerun=True)
test = LakeModel.pull_output_nc(test, force=True)
test.score_variant(Lake_temps, test, 'temp', lik_fxn='NSE')

LakeModel.plotLakeandError(test, Lake_temps)
#This is the potential model space 
test.read_variants('../test_files/optimize_these.txt', verbose=False)
#Here we want to create a new lake at a random starting point

variant_lakes = test.write_variant_configs(copycsvs=True, verbose=False)

#Run them all & score them all
variant_lakes_w_data = []
for i, l in enumerate(variant_lakes):
    statvfs = os.statvfs(os.getcwd())
    bytes_avail = statvfs.f_frsize * statvfs.f_bavail
    if bytes_avail < 104857600/2:
        sys.exit("Disk space less than 100 Mb, aborting")
    print "~{} Mb of disk space remaining".format(bytes_avail/1024/1024)    
    print "Running #", i+1, "out of", len(variant_lakes)
    print "%r = %r" % (l.variant_var, l.variant_val)
    ran_lake = LakeModel.run_model(l, verbose=True, forceRerun=True)
    data_lake = LakeModel.pull_output_nc(ran_lake, force=True)
    if data_lake.ran == True:
        variant_lakes_w_data.append(data_lake)
        try: 
            data_lake.score_variant(Lake_temps, test, 'temp', lik_fxn='NSE')
        except ValueError:
            sys.exit()
        data_lake.cleanup(verbose=True)
    else:
        print "Variant %s failed and was removed" % data_lake.glm_config['glm_setup']['sim_name']

# Find the best improvement 
vars_lik = {i.liklihood: (i.variant_var, i.variant_val) for i in variant_lakes_w_data if ~np.isnan(i.liklihood)}
max_lik = np.array(vars_lik.keys()).max()
print vars_lik[max_lik], max_lik

sys.exit()

start_time = time.time()

rand_case_name = "randomize_4"
total_trials = 25250
max_pickle_size = 500
repeats = np.ceil(total_trials/ max_pickle_size)
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
            this_err = random_lake.score_variant(Lake_temps, test, 'temp', 
                                                 lik_fxn='NSE')
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


"""
This is a plot of inflow, precip, lake, and cumulative inflow vols
plt.figure(figsize=(12,8))
plt.plot(test.csv_dict['lake.csv']['Volume'], label = "Lake Volume", alpha=0.7)
ax = plt.gca()
ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
plt.plot(test.csv_dict['lake.csv']['Tot Inflow Vol'], label = "Inflow Volume", alpha=0.7)
rain, _ = subsectbydate_2(BOS_weather.df.precip, test.csv_dict['lake.csv']['Volume'])
plt.plot(rain*10e7, label = 'Precipitation', alpha=0.7)
cumInVol = test.csv_dict['lake.csv']['Tot Inflow Vol'].cumsum()
ymin, ymax = plt.ylim()
plt.plot(cumInVol, label = "Cumulative Inflow Vol", alpha=0.7)
plt.ylim([ymin, ymax])
xmin, xmax = plt.xlim()
plt.xlim([xmin-10, xmax])
plt.legend(loc=1)
plt.ylabel('cubic meters, cubic meters/day or mm/day')
"""


"""
#plots the given inflow and the processed output inflow, needs modified alpha
plt.plot(test.csv_dict['lake.csv']['Tot Outflow Vol'], label="Outflow")
plt.plot(test.csv_dict['lake.csv']['Tot Inflow Vol'], label="Inflow")
plt.plot(test.csv_dict['lake.csv']['Overflow Vol'], label="Overflow")
plt.plot(100000*In_Discharge[In_Discharge.index[-1551:(-1551+731)]], label="A priori inflow")


if plot_now:
    plt.figure(1)
    plt.plot(test.csv_dict['lake']['LakeNumber'], c='g', label="Lake Number", alpha=0.5)
    plt.legend(loc='best')
    plt.figure(2)
    plt.plot(test.csv_dict['lake']['Tot Inflow Vol'], c='r', label="Inflow Vol", alpha=0.5)
    plt.plot(test.csv_dict['lake']['Tot Outflow Vol'], c='b', label="Outflow Vol", alpha=0.5)
    plt.legend(loc='best')
    plt.figure(3)
    plt.plot(test.csv_dict['lake']['Evaporation'], c='m', label="Evaporation", alpha=0.5)
    plt.legend(loc='best')
    plt.figure(4)
    plt.plot(test.csv_dict['lake']['Max Temp'], label='Max Temp')
    plt.plot(test.csv_dict['lake']['Min Temp'], label='Min Temp')
    plt.plot(test.csv_dict['lake']['Surface Temp'], label='Surface Temp')
    plt.legend()
    

base_temp = test.output_nc['temp'][:,:,0,0].data
variant_temps = [vLr.output_nc['temp'][:,:,0,0].data for vLr in variant_lakes_ran]
base_temp[base_temp > 1000] = 0
for vLr in variant_temps:
    vLr[vLr > 1000] = 0

error = [((vt - base_temp)**2).sum() for vt in variant_temps]
variant_names = [n.name for n in variant_lakes_ran]
error_df = pd.DataFrame(index = variant_names, columns = ['error'],
                        data = error)
error_df.sort_values('error', ascending=False, inplace=True)

1. atm_stab
2. max_layer_thick (2-21)
3. strm_hf_angle (10-16) (not in order)
4. min_layer_thick
5. strmbd_slope
6. lw_factor
7. wind_factor
8. rain_factor
9. seepage_rate
10. rh_factor

plt.figure()
for k in test.variants['likelihood'].keys():
    grid_size = test.variants['likelihood'][k].shape[0]
    if grid_size > 11:
        plt.plot(np.arange(grid_size), test.variants['likelihood'][k], label=k)
plt.legend(loc='best')

counter = 0
fig=plt.figure(figsize=(11,8))
for k in test.variants['likelihood'].keys():
    grid_size = test.variants['likelihood'][k].shape[0]
    if counter == 0:
        ax=fig.add_subplot(2,1,1)
        ax.set_ylabel('likelihood')
    if counter == 6:
        ax.legend(loc='best')
        ax=fig.add_subplot(2,1,2)
        ax.set_ylabel('likelihood')     
    if grid_size == 11:
        counter+=1
        ax.plot(np.arange(grid_size), test.variants['likelihood'][k], label=k)

ax.legend(loc='best')
"""


