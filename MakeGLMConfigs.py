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

plt.style.use('fivethirtyeight')


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
mainFolder = '/home/login/GDrive/Documents/Landscape_Hydrology/Final Project'

# GLM starts from the command line and immediately searches in the current 
# directory for `glm2.nml`. So we need a directory for each case and a 
# directory to hold them all which will be:

superDir = '/glm_case_folders'
make_dir(mainFolder+superDir)

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

met_path = os.path.join(mainFolder, "BostonLoganAirportWeather.csv")
ceres_fn = "CERES_SSF_XTRK-MODIS_Edition3A_Subset_2010010102-2014010116.nc"
ceres_path = os.path.join(mainFolder, "CERES", ceres_fn)
GEODISC_path = os.path.join(os.getcwd(), "Giovanni", "test_files")
adapt_data = os.path.join(GEODISC_path, 'test_data_adaptive.csv')
adapt_data3n = os.path.join(GEODISC_path, 'test_data_adaptive3n.csv')


test.read_GCHN(met_path)
test.read_CERES_nc(ceres_path)
test.read_GEODISC_cloud_data('cloud_data_5pt.csv',  GEODISC_path)

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


