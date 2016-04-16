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
import os
import LakeModel
import numpy as np

def addtodict(vals, vars, dict):
    for var, val in zip(vars, vals):
        dict[var] = val
    
def make_dir(s):
    if os.path.exists(s):
        #sys.exit("I cannot build a house atop another house")
        print "%s exists\n" % os.path.basename(s)
    else:
        os.mkdir( s, 0760)


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
met_path = mainFolder + "/" + "BostonLoganAirportWeather.csv"
ceres_path =mainFolder+"/CERES/"+"CERES_SSF_XTRK-MODIS_Edition3A_Subset_2010010102-2014010116.nc"

test.read_GCHN(met_path)
test.read_CERES_nc(ceres_path)
# each input file will be manually input as a file argument for vetting
# each will be accessible as a data frame object 

#input inflow.csv



