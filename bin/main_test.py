# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 21:02:22 2016

This is a script to autogenerate all the `glm.nml` files for all the cases
The first case will be sunny mild day with no variability


To run glm, you need to get the file `l_fcompxe_2013.2.146_redist`
and later use the command when terminal starts:

    source /opt/intel/bin/compilervars.sh intel64

@author: Keith Arora-Williams
"""
from __future__ import division
import time
import os
import LakeModel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plotting = False
plt.style.use('fivethirtyeight')
print "\nGLM_WRAPPER"
print time.ctime(time.time())
print ""

# these are code names for each trial
cases = ["testcase"]

# this is the directory you want to work in
mainFolder = os.path.dirname(os.getcwd())

# GLM starts from the command line and immediately searches in the current 
# directory for `glm2.nml`. So we need a directory for each case and a 
# directory to hold them all which will be:

superDir = '/glm_case_folders'
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

newDirs = [mainFolder+superDir+"/"+x for x in cases]
LakeModel.make_dir(mainFolder+superDir)
test = LakeModel.Lake(cases[0], newDirs[0])

test.write_glm_config()

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

waterOjbects = [Hobbs_o, Hobbs_i1, Hobbs_i2, Hobbs_i3, Hobbs_i4, inflow, 
                mystic_d, mystic_t, GroundW]

for h in waterOjbects:
    h.preprocess()
    h.read()
    h.removeQuals()
    h.print_metadata()

GroundW. df['GroundwaterLevel'] = (GroundW.df[GroundW.df.columns[0]]-96.6)*-1

Lake_temps = pd.read_csv(waterDataPaths['ML'], index_col=0, 
                         parse_dates=['datetime'])

CERES_data.df1.cldarea_total_daily.fillna(method='pad', inplace=True)

SSF_D = CERES_SSF.ceres_df.resample('D').mean()
SSF_Filled = SSF_D.interpolate()
filledagain = SSF_Filled.fillna(method='bfill')

#Pack up met.csv variables
longWaveIn = CERES_data.df1['sfc_comp_lw-down_all_daily']
longWaveIn.name = "LongWave"
shortWaveNet = CERES_data.df1['sfc_comp_sw-down_all_daily'] - CERES_data.df1['sfc_comp_sw-up_all_daily']
shortWaveNet.name = "ShortWave"
airTemp = (BOS_weather.df.t_max + BOS_weather.df.t_min) / 2.0
airTemp.name = "AirTemp"
relHumidity = SSF_Filled.humidity
relHumidity.name = "RelHum"
windSpeed = BOS_weather.df.wind_speed
windSpeed.name = "WindSpeed"
precipR = BOS_weather.df.precip
precipR.name = "Rain"
precipS = BOS_weather.df.snow_fall
precipS.name = "Snow"

# Pack up variables for flow csvs 
Out_Discharge = inflow.df['Discharge, cubic feet per second (Mean)']*0
Out_Discharge.name = "OutDischarge"
In_Discharge = inflow.df['Discharge, cubic feet per second (Mean)']
In_Discharge.name = "InDischarge"
In_Temp = Hobbs_o.df['Temperature, water, degrees Celsius (Mean)'].interpolate()
In_Temp.name = "InflowTemp"
In_Salt = Hobbs_o.df[Hobbs_o.df.columns[6]].interpolate()
In_Salt.name = "InflowSalt"

data_pack = [shortWaveNet, longWaveIn, airTemp, relHumidity, windSpeed, 
             precipR, precipS, Out_Discharge, In_Discharge, In_Temp, In_Salt]
             
test.temporal_clipping(data_pack)

metPack = ['ShortWave','LongWave','AirTemp','RelHum','WindSpeed','Rain','Snow']
inPack = ['InDischarge','InflowTemp','InflowSalt']
outPack = ['OutDischarge']
test.create_metcsv(metPack)
test.create_flowcsvs(inPack, outPack)








