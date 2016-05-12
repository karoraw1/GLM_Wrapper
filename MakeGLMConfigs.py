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
import time
import os
import LakeModel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr

plotting = False

# Kling-Gupta Efficiency
def error_metrics(Obs, Sim):
    Err = Sim - Obs
    MSE = np.mean(Err**2)
    rMSE = np.sqrt(MSE)
    NSE = 1 - (MSE/np.var(Obs))
    r = pearsonr(Obs, Sim)
    a = np.std(Sim) / np.std(Obs)
    b = np.mean(Sim) / np.mean(Obs)
    KGE = 1- np.sqrt((r[0]-1)**2 + (a-1)**2 + (b-1)**2)
    return {'MSE':MSE, 'rMSE':rMSE, 'NSE':NSE, 'r':r, 'a':a, 'b':b, 'KGE':KGE}
    
plt.style.use('fivethirtyeight')

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
    
def make_dir(s):
    if os.path.exists(s):
        #sys.exit("I cannot build a house atop another house")
        print "%s exists\n" % os.path.basename(s)
    else:
        os.mkdir( s, 0760)

print "\nGLM_WRAPPER"
print time.ctime(time.time())
print ""

# these are code names for each trial
cases = ["sunnyDay"]

# this is the directory you want to work in
mainFolder = os.getcwd()

# GLM starts from the command line and immediately searches in the current 
# directory for `glm2.nml`. So we need a directory for each case and a 
# directory to hold them all which will be:

superDir = '/glm_case_folders'
make_dir(mainFolder+superDir)
met_fn = "BostonLoganAirportWeather.csv"
met_path = os.path.join(mainFolder, 'weatherData', met_fn)

ceres_fn = "CERES_SSF_XTRK-MODIS_Edition3A_Subset_2010010102-2014010116.nc"
ceres_ssf = os.path.join(mainFolder, 'weatherData', ceres_fn)
ceres_fn2 = "CERES_SYN1deg-Day_200508-201406.nc"
ceres_SYN1 = os.path.join(mainFolder, 'weatherData', ceres_fn2)

GEODISC_path = os.path.join(mainFolder, "Giovanni", "test_files")
adapt_data = os.path.join(GEODISC_path, 'test_data_adaptive.csv')
adapt_data3n = os.path.join(GEODISC_path, 'test_data_adaptive3n.csv')

context_pond = "Fresh_pond_air_and_water_temp.txt"
inflow_fn = "aberjona_15min_discharge_precip_streamU_gageheight.txt"
outflow_fn1 = "mystic_river_15min_water_temp_discharge2015_2016.txt"
outflow_fn2 = "alewifedailydischarge2005-2016.txt"
hobbs_fn = "HobbsBk_all.txt"
flow_fns = [inflow_fn, outflow_fn1, outflow_fn2, hobbs_fn, context_pond]
flow_ps = [os.path.join(os.getcwd(), 'waterdata', i) for i in flow_fns]

newDirs = [mainFolder+superDir+"/"+x for x in cases]


# the glm contains the following blocks:
# &glm_setup: General simulation info and mixing parameters

setup_vars = ['max_layers', 'min_layer_vol', 'min_layer_thick', 
              'max_layer_thick', 'Kw', 'coef_inf_entrain', 'coef_mix_conv', 
              'coef_wind_stir', 'coef_mix_shear', 'coef_mix_turb', 
              'coef_mix_KH', 'coef_mix_hyp']
              # kw is the extinction coefficient of PAR
wq_vars = ['wq_lib', 'ode_method', 'split_factor', 'bioshade_feedback', 
           'repair_state', 'multi_ben']           
morpho_vars = ["lake_name", "latitude", "longitude",
               "bsn_len", "bsn_wid", "bsn_vals"]
time_vars = ["timefmt", "start", "stop", "dt", "timezone"]
output_vars = ["out_dir", "out_fn", "nsave", "csv_lake_fname", "csv_point_nlevs",
            "csv_point_fname", "csv_point_at", "csv_point_nvars",
            "csv_point_vars", "csv_outlet_allinone", "csv_outlet_fname",
            "csv_outlet_nvars", "csv_outlet_vars", "csv_ovrflw_fname"]
init_vars = ["lake_depth", "num_depths", "the_depths", "the_temps", "the_sals",
             "num_wq_vars", "wq_names", "wq_init_vals" ]
met_vars = ["met_sw", "lw_type", "rain_sw", "snow_sw", "atm_stab", "catchrain",
            "rad_mode", "albedo_mode", "cloud_mode", "subdaily", "meteo_fl",
            "wind_factor", "sw_factor", "lw_factor", "at_factor", "rh_factor",
            "rain_factor", "ce", "ch", "cd", "rain_threshold", "runoff_coef"]
bird_vars = ["AP", "Oz", "WatVap", "AOD500", "AOD380", "Albedo"]
outflow_vars = ["num_outlet", "flt_off_sw", "outl_elvs", "bsn_len_outl", 
            "bsn_wid_outl", "outflow_fl", "outflow_factor"]
inflow_vars = ["num_inflows", "names_of_strms", "subm_flag", "strm_hf_angle", 
               "strmbd_slope", "strmbd_drag", "inflow_factor", "inflow_fl", 
               "inflow_varnum", "inflow_vars", "coef_inf_entrain"]

inflow_vals = [1, "'Aberjona'", ".false.", 65.0, 2.0, 0.0160, 1.0, 
               "'inflow.csv'", 4, ['FLOW', 'TEMP', 'SALT', 'OXY_oxy', 
               'SIL_rsi', 'NIT_amm', 'NIT_nit', 'PHS_frp', 'OGM_don', 
               'OGM_pon', 'OGM_dop', 'OGM_pop', 'OGM_doc', 'OGM_poc', 
               'PHY_green', 'PHY_crypto', 'PHY_diatom'], 0.]
outflow_vals = [1, ".false.", -215.5, 799, 399, "'outflow.csv'", 0.8]
bird_vals = [973, 0.279, 1.1, 0.033, 0.038, 0.2]
init_vals = [ 22.0, 5, [1, 5, 9, 13, 17, 21],
              [4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
              6, ["'OGM_don',", "'OGM_pon',",
                  "'OGM_dop',", "'OGM_pop',",
                  "'OGM_doc',", "'OGM_poc'"],
              np.array([[1.1, 1.2, 1.3, 1.2, 1.3],
                        [2.1, 2.2, 2.3, 1.2, 1.3],
                        [3.1, 3.2, 3.3, 1.2, 1.3],
                        [4.1, 4.2, 4.3, 1.2, 1.3],
                        [5.1, 5.2, 5.3, 1.2, 1.3],
                        [6.1, 6.2, 6.3, 1.2, 1.3]])]

setup_vals = [500, 0.025, 0.50, 1.500, 0.6, 0., 0.125, 
              0.23, 0.20, 0.51, 0.30, 0.5]
wq_vals = ["'aed2'", 1, 1, '.true.', '.true.', "'fabm.nml'", '.true.']             
morpho_vals = ["'UpperMysticLake'", 42.4317, -71.1483, 
               "1073.637,", "632.60,", 24]
time_vals = [2, "'2012-01-01 00:00:00'", "'2014-01-01 00:00:00'", 3600.0, 5.0]
output_vals = ["", "out", 12, "'lake'", 1, "'WQ_'", "17.", 2, 
            ["'temp',", "'salt',", "'OXY_oxy',"], ".false.", 'outlet_', 3,
            ["'flow',", "'temp',", "'salt',", "'OXY_oxy',"], "\"overflow\""]
met_vals = [".true.", "'LW_IN'", ".false.", ".false.", ".false.", ".false.", 
            1, 1, 4, ".false.", "'met_hourly.csv'", 1.0, 1.0, 1.0, 1.0, 1.0, 
            1.0, 0.0013, 0.0013, 0.0013, 0.01, 0.3]

value_dict = {var:val for var, val in zip(setup_vars, setup_vals)}
addtodict(wq_vals, wq_vars, value_dict)
addtodict(morpho_vals, morpho_vars, value_dict)
addtodict(time_vals, time_vars, value_dict)
addtodict(output_vals, output_vars, value_dict)
addtodict(init_vals, init_vars, value_dict)
addtodict(met_vals, met_vars, value_dict)
addtodict(bird_vals, bird_vars, value_dict)
addtodict(outflow_vals, outflow_vars, value_dict)
addtodict(inflow_vals, inflow_vars, value_dict)

test = LakeModel.Lake(cases[0], newDirs[0], value_dict)

BOS_weather = LakeModel.GHCN_weather_data(met_path)
BOS_weather.clean_columns()

CERES_SSF = LakeModel.CERES_nc(ceres_ssf, "3")
CERES_data = LakeModel.CERES_nc(ceres_SYN1, "4")
test.read_GEODISC_cloud_data('cloud_data_5pt.csv',  GEODISC_path)

inflow = LakeModel.USGS_water_data(flow_ps[0])
inflow.preprocess()
inflow.read()
inflow.print_metadata()

mystic = LakeModel.USGS_water_data(flow_ps[1])
mystic.preprocess()
mystic.read()
mystic.print_metadata()

alewife = LakeModel.USGS_water_data(flow_ps[2])
alewife.preprocess()
alewife.read()
alewife.print_metadata()

FreshPond = LakeModel.USGS_water_data(flow_ps[4])
FreshPond.preprocess()
FreshPond.read()

Hobbs = LakeModel.USGS_water_data(flow_ps[3])
Hobbs.preprocess()
Hobbs.read()
Hobbs.print_metadata()

if plotting == True:
    old_cols = Hobbs.df.columns.values
    old_cols[7] = 'Water Temperature'
    old_cols[15] = 'Air Temperature'
    old_cols[9] = 'Specific conductance'
    old_cols[17] = 'Precipitation'
    
    Hobbs.df.columns = old_cols
    to_drop = list(old_cols)
    to_drop.remove(old_cols[7])
    to_drop.remove(old_cols[15])
    to_drop.remove(old_cols[9])
    to_drop.remove(old_cols[17])
    temp_df = Hobbs.df.drop(to_drop, axis=1)
    temp_df=temp_df.fillna(method="pad")
    temp_ser3 = temp_df['Precipitation'].resample('D').sum()
    temp_df = temp_df.resample('D').mean()
    temp_df2 = subsectbydate_1(FreshPond.df, temp_df)
    temp_df2=temp_df2.fillna(method="pad")
    
    air = error_metrics(temp_df2['Temperature, air'].values,
                        temp_df['Air Temperature'].values)
    cond = error_metrics(temp_df2['Specific conductance'].values,
                         temp_df['Specific conductance'].values)
    water = error_metrics(temp_df2['Temperature, water'].values,
                          temp_df['Water Temperature'].values)
    precip = error_metrics(temp_df2['Precipitation, total'].values,
                           temp_ser3.values)


    plt.style.use('ggplot')
    f, ((ax3, ax2), (ax1, ax4)) = plt.subplots(2, 2, sharex='col', figsize=(24,12))
    ax4.plot(temp_ser3)
    ax4.plot(temp_df2['Precipitation, total'])
    ax4.set_title('Precipitation', fontsize=20)
    ax3.plot(temp_df['Air Temperature'])
    ax3.plot(temp_df2['Temperature, air'])
    ax2.plot(temp_df['Water Temperature'])
    ax2.plot(temp_df2['Temperature, water'])
    ax1.plot(temp_df['Specific conductance'])
    ax1.plot(temp_df2['Specific conductance'])
    ax1.set_title('Specific conductance', fontsize=20)
    ax2.set_title('Water Temperature', fontsize=20)
    ax3.set_title('Air Temperature', fontsize=20)
    ax4.set_alpha(0.5)
    
    
    ax1.text(733238, 1700, "KSE: "+"{:04.2f}".format(cond['KGE']), fontsize=18)
    ax1.text(733230.5, 1600, "NSE: "+"{:04.2f}".format(cond['NSE']), fontsize=18)
    ax1.set_ylabel('microSiemens/cm @ 25C', fontsize=20)
    
    ax2.text(732953, 2.5, "KSE: "+"{:04.2f}".format(water['KGE']), fontsize=18)
    ax2.text(732953, 0.5, "NSE: "+"{:04.2f}".format(water['NSE']), fontsize=18)
    ax2.set_ylabel('degrees C', fontsize=20)
    
    ax3.text(733238, -12, "KSE: "+"{:04.2f}".format(air['KGE']), fontsize=18)
    ax3.text(733238, -14.5, "NSE: "+"{:04.2f}".format(air['NSE']), fontsize=18)
    ax3.set_ylabel('degrees C', fontsize=20)
    
    ax4.text(732953, 2.8, "KSE: "+"{:04.2f}".format(precip['KGE']), fontsize=18)
    ax4.text(732953, 2.62, "NSE: "+"{:04.2f}".format(precip['NSE']), fontsize=18)
    ax4.set_ylabel('inches (sum)', fontsize=20)
    f.legend( ax1.get_lines(), ['Hobbs Brook', 'Fresh Pond'], loc = (0.405, .465), 
             fontsize=18, ncol=2)
    
    ## Saving is broken for some reason. Copying and pasting the above and below
    ## code into ipython sequentially works though
    #labs = ax4.get_xticklabels()
    #ax1.set_xticklabels(labs, fontsize=12, color='k')
    #ax4.set_xticklabels(labs, fontsize=12, color='k')
    #f.savefig("HobbsandFresh.png")

#test.mergedata()

if plotting == True:
    plt.style.use('fivethirtyeight')
    adaptive_n = pd.read_csv(adapt_data, index_col = 0, 
                             parse_dates = ['start', 'end'], 
                             infer_datetime_format = True)
    adaptive_3n = pd.read_csv(adapt_data3n, index_col = 0, 
                              parse_dates = ['start', 'end'], 
                              infer_datetime_format = True)
    
    adaptive = pd.concat([adaptive_3n, adaptive_n], ignore_index=True)
    flat_n = test.geodisc_clouds.sample(n=160)
    xlabs = ["Std. Deviation per Swath(n)", "Measurements per Swath (n)" ]
    col_names = ['cc_stds', 'cc_n']
    fignames = ['adaptiveVflat_stds.png', 'adaptiveVflat_n.png']
    colors1 = ['r', 'g']
    
    for xlab, col_name, fig_n, cs in zip(xlabs, col_names, fignames, colors1):
        s1 = flat_n[col_name].values
        s2 = adaptive[col_name].values
        plot_vals = pd.Series(np.concatenate((s1,s2)))
        labs1 = np.tile(np.array('0.3 Degree Threshold'), (160,))
        labs2 = np.tile(np.array('Non-Zero Minimum Threshold'), (160,))
        plot_grps = np.concatenate((labs1, labs2))
        ax1, ax2 = plot_vals.hist(by=plot_grps, figsize=(12,6), color=cs)
        ax1.set_xlabel(xlab)
        ax2.set_xlabel(xlab)
        ax1.set_ylabel("Swath Count")
        plt.savefig(fig_n, dpi=300)


def getStartTimes(df, col1_i, i):
    a_time  = df.iloc[i, col1_i]
    total_s = a_time.hour * 3600 + a_time.minute * 60 + a_time.second
    return total_s/3600.

def getPeriodicity(df, col1_i, time_d, i):
    date_format = '%Y-%m-%d %H:%M:%S'
    if i < time_d:
        return np.nan
    else:
        end_s = str(test.geodisc_clouds.iloc[i, col1_i])
        end_t = pd.datetime.strptime(end_s, date_format)
        start_s = str(test.geodisc_clouds.iloc[i-time_d, col1_i])
        start_t = pd.datetime.strptime(start_s, date_format)
    return (end_t - start_t).seconds/60./60.


if plotting == True:
    nondataCols = ['start', 'end', 'hdf5', 'i_out_2', 'checksum']
    cc_plt_df = test.geodisc_clouds.drop(nondataCols, axis=1)
    cc_plt_df.columns = ['Mean Cloud Cover', 'Cloud Cover STD',
                         'Cloud Cover Sample Size (n)', 'Latitude Std (deg)', 
                         'Longitude Std (deg)' ]
    
    rows, cols = cc_plt_df.shape
    cc_plt_df.insert(0, 'Time of Day (hrs)', np.zeros(rows))
    
    for n in range(rows):
        cc_plt_df.iloc[n, 0] = getStartTimes(test.geodisc_clouds, 1, n)
    
    cc_plt_df.hist(bins=30)

sorted_df = test.geodisc_clouds.sort_values('start', axis = 0, ascending = True)
Geodisc_clouds = sorted_df.copy()
Geodisc_clouds.set_index('start', drop=True, inplace=True, verify_integrity=True)


errors = {}
lon1df = CERES_data.df1.cldarea_total_daily.fillna(method='pad')
lon2df = CERES_data.df2.cldarea_total_daily.fillna(method='pad')
lon_dfs = [lon1df,lon2df ]
lon_keys = list(CERES_data.lons_n)
Geo_filled=Geodisc_clouds.interpolate()
for srs, key in zip(lon_dfs, lon_keys):
    cld_C, cld_G = subsectbydate_2(srs,Geo_filled.cc_means)
    cld_1 = cld_G.resample('D').mean()
    cld_1 = cld_1.interpolate()*100.
    errors[key] = error_metrics(cld_1, cld_C)

lon1df_t, _ = subsectbydate_2(lon1df,Geo_filled.cc_means)
lon2df_t, _ = subsectbydate_2(lon2df,Geo_filled.cc_means)
errors['other'] = error_metrics(lon1df_t, lon2df_t)


SSF_D = Geodisc_clouds.resample('D').mean()
Geo_D = CERES_SSF.ceres_df.resample('D').mean()
Geo_filled=SSF_D.interpolate()
SSF_Filled=Geo_D.interpolate()
cld_GE, cld_SSF = subsectbydate_2(Geo_filled,SSF_Filled.fillna(method='bfill'))
cld_GE, cld_SSF = subsectbydate_2(Geo_filled,SSF_Filled.fillna(method='bfill'))
errors['SSFvGEO'] = error_metrics(cld_SSF.cloud_frac.values, 
                                  (cld_GE.cc_means*cld_SSF.cloud_frac.max()).values)
lon_keys2 = ['df1', 'df2']
filledagain = SSF_Filled.fillna(method='bfill')
for srs, key in zip(lon_dfs, lon_keys2):
    cld_C, cld_G = subsectbydate_2(srs,filledagain.cloud_frac)
    cld_1 = cld_G.resample('D').mean()
    cld_1 = cld_1.interpolate()*cld_C.max()
    errors[key] = error_metrics(cld_1, cld_C)

if plotting == True:    
    f , ax5 = plt.subplots(1, 1, figsize=(12,12))
    ax5.plot(cld_SSF.cloud_frac, label="CERES SSF Data")
    ax5.plot((cld_GE.cc_means*cld_SSF.cloud_frac.max()), label="Goddard OMI Data")
    ax5.legend()
    ax5.set_ylim((-10, 100))
    ax5.set_xlim(734503+90,734503+185)
    f.savefig("CERES_SSFvOMMYCLD0.png")

if plotting == True:    
    sorted_df.insert(0, 'periodicity (hours)', np.zeros(rows))
    
    for m in range(rows):
        sorted_df.iloc[m, 0] = getPeriodicity(test.geodisc_clouds, 1, 1, m)
    
    sorted_df.set_index('start', drop=True, inplace=True, verify_integrity=True)
    sorted_df_pre_agg, new_idx_cols = LakeModel.TimeIdx(sorted_df)
    null_col = sorted_df_pre_agg.notnull().cc_means
    sorted_df_pre_agg['Null Values'] = null_col
    cc_aggs_mu = LakeModel.makeAggregations(sorted_df_pre_agg, new_idx_cols, np.mean)
    cc_aggs_std = LakeModel.makeAggregations(sorted_df_pre_agg, new_idx_cols, np.std )
    cc_aggs_sum = LakeModel.makeAggregations(sorted_df_pre_agg, new_idx_cols, np.sum)
    
    plt.figure(3)
    plt.clf()
    plt.scatter(flat_n.cc_stds , flat_n.cc_n, c = 'r')
    plt.scatter(adaptive.cc_stds , adaptive.cc_n, c = 'b')
    plt.ylabel('Measurements (n)')
    plt.xlabel('Standard Deviation (n)')



#
#import matplotlib.pyplot as plt
#plt.figure(1)
#plt.clf()
#im_path = "/home/login/GDrive/Documents/Landscape_Hydrology/Final Project/Presentation"
#im_name = os.path.join(im_path, "CERESBB.png")
#im = plt.imread(im_name)
#
#implot = plt.imshow(im)
## put a red dot, size 40, at 2 locations:
#
##
##plt.scatter(x=[10], y=[10], marker="x", s=10)
#xlim = 380/lat_c.mean()
#ylim = 480/lon_c.mean()
#plt.scatter(x=lat_c*xlim, y=lon_c*ylim, marker="x", s=10)
#plt.tight_layout()
#plt.show()

# TODO: plot aggregations and original data on each plot 
# TODO: use fft to determine variables within time series
# TODO: impute interpolated data based on autocorrelation
# TODO: read in inflow data
# TODO: read in outflow data
# TODO: maximize cross-correlation of identical data from different sources
# TODO: write out properly formatted input files based on input data
# TODO: write class method that runs glm by spawning independent processes
# TODO: test and debug execution
# TODO: run once and time it and measure memory footprint 
# TODO: write a function that generates & submits a batch file
# TODO: 
# 
# 
# each input file will be manually input as a file argument for vetting
# each will be accessible as a data frame object 


