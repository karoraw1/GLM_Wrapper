# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:03:39 2016

@author: login
"""
import os
import numpy as np

def make_dir(s):
        if os.path.exists(s):
            #sys.exit("I cannot build a house atop another house")
            print "%s exists\n" % os.path.basename(s)
        else:
            os.mkdir( s, 0760)
            print "\n made dir %s" % os.path.basename(s)

def safe_dir(s):
    counter = 1
    # check if dir exists
    if os.path.exists(s):
        # check to see if it is empty
        if os.listdir(s):
            # if it is, add integers until a fresh name is found
            while os.path.exists(s) and counter < 1000:
                s = s + str(counter)
                counter +=1
            os.mkdir(s, 0760)
            print "\nNew folder name found & folder made"
        else:
            # if empty, skip ahead
            print "\nEmpty folder already exists, skipping creation"
    else:
        #if it doesn't exist, make it
        os.mkdir(s, 0760)
        print "\n made dir %s" % os.path.basename(s)
            
class Lake(object):
    
    # &glm_setup: General simulation info and mixing parameters
    setup_vars = ['sim_name', 'max_layers', 'min_layer_vol', 'min_layer_thick', 
                  'max_layer_thick', 'Kw', 'coef_inf_entrain', 'coef_mix_conv', 
                  'coef_wind_stir', 'coef_mix_shear', 'coef_mix_turb', 
                  'coef_mix_KH', 'coef_mix_hyp']
    # &wq_setup: Details about the coupling with the water quality model (eg. FABM or AED2)
    wq_vars = ['wq_lib', 'ode_method', 'split_factor', 'bioshade_feedback', 
               'repair_state', 'multi_ben']
    # &time: Time controls
    time_vars = ["timefmt", "start", "stop", "dt", "timezone"]
    # &morphometry: Lake morphometric information
    morpho_vars = ["lake_name", "latitude", "longitude", "bsn_len", "bsn_wid",
                   "bsn_vals"]
    # &output: Specification of output file details
    out_vars = ["out_dir", "out_fn", "nsave", "csv_lake_fname", 
                "csv_point_nlevs", "csv_point_fname", "csv_point_at", 
                "csv_point_nvars", "csv_point_vars", "csv_outlet_allinone", 
                "csv_outlet_fname", "csv_outlet_nvars", "csv_outlet_vars", 
                "csv_ovrflw_fname"]
    # &init_profiles: Setting initial conditions (depth profiles) of GLM and WQ variables
    init_vars = ["lake_depth", "num_depths", "the_depths", "the_temps", "the_sals",
             "num_wq_vars", "wq_names", "wq_init_vals" ]
    # &meteorology: Information about surface forcing and meteorology data
    met_vars = ["met_sw", "lw_type", "rain_sw", "snow_sw", "atm_stab",
                "catchrain", "rad_mode", "albedo_mode", "cloud_mode",
                "subdaily", "meteo_fl", "wind_factor", "sw_factor",
                "lw_factor", "at_factor", "rh_factor", "rain_factor", "ce",
                "ch", "cd", "rain_threshold", "runoff_coef"]
    # &inflows: Information about inflows

    # &outflows: Information about outflows
    out_vars = ["num_outlet", "flt_off_sw", "outl_elvs", "bsn_len_outl", 
                "bsn_wid_outl", "outflow_fl", "outflow_factor"]
    # &bird: Optional block to input parameters for the Bird solar radiation model
    bird_vars = ["AP", "Oz", "WatVap", "AOD500", "AOD380", "Albedo"]
    
    def __init__(self, name, dir_path, value_dict):
        
        self.area = np.array([  77373.8, 99499.81098205, 124404.590076,
                              150649.64403947, 168609.73032884, 187580.97034902,
                              207020.52054899, 225715.75772409, 245219.05086584,
                              267726.76155245, 295060.85089366, 323723.31628528,
                              348395.92365881, 368505.54150352, 389179.36005101,
                              410522.51421128, 432500.14410228, 455050.82895077,
                              477078.71675217, 499253.70320025, 521932.40482378,
                              535409.21116882, 547660.91761877, 560051.22])
        self.elevation =  np.array([-23.384, -22.32382609, -21.26365217, 
                                    -20.20347826, -19.14330435, -18.08313043, 
                                    -17.02295652, -15.96278261, -14.9026087, 
                                    -13.84243478, -12.78226087, -11.72208696,
                                    -10.66191304, -9.60173913, -8.54156522,
                                    -7.4813913, -6.42121739, -5.36104348, 
                                    -4.30086957, -3.24069565, -2.18052174, 
                                    -1.12034783, -0.06017391, 1.])

        self.name = name
        self.dir_path = dir_path
        make_dir(dir_path)
        self.glm_path = self.dir_path+"/glm2.nml"
        self.value_dict = value_dict
        self.value_dict['sim_name'] = "'"+name+"'"
        self.value_dict["out_fn"] = value_dict["out_fn"] + "_"+name
        #make output dir
        self.value_dict["out_dir"] = dir_path+"/output"
        safe_dir(value_dict["out_dir"])
        # start writing config file
        glm_handle = open(self.glm_path, 'w+')    
        
        glm_handle.write("&glm_setup\n")
        for var in self.setup_vars:
            glm_handle.write("\t{0} = {1}\n".format(var, value_dict[var]))

        glm_handle.write("/\n&wq_setup\n")
        for var in self.wq_vars:
            glm_handle.write("\t{0} = {1}\n".format(var, value_dict[var]))

        glm_handle.write("/\n&morphometry\n")
        for var in self.morpho_vars:
            glm_handle.write("\t{0} = {1}\n".format(var, value_dict[var]))

        # Input Bathy data
        glm_handle.write("\tH = ")
        layers = range(len(self.elevation))
        for idx in layers:
            if ((idx%4) == 0) and idx !=0:
                glm_handle.write("\n\t    ")            
            if idx !=layers[-1]:
                glm_handle.write("{0}, ".format(self.elevation[idx]))
            else:
                glm_handle.write("{0}\n".format(self.elevation[idx]))
        glm_handle.write("\tA = ")

        for idx in layers:
            if ((idx%4) == 0) and idx !=0:
                glm_handle.write("\n\t    ")
            if idx !=layers[-1]:
                glm_handle.write("{0}, ".format(self.area[idx]))
            else:
                glm_handle.write("{0}\n/\n".format(self.area[idx]))
        glm_handle.write("&time\n")
        for var in self.time_vars:
            glm_handle.write("\t{0} = {1}\n".format(var, value_dict[var]))
        glm_handle.write("/\n&output\n")
        
        for var in self.out_vars:
            val = value_dict[var]            
            if type(val) == list:
                glm_handle.write("\t{0} = {1}\n".format(var, val[0]))
                for idx in np.array(range(len(val)-1))+1:
                    glm_handle.write("\t\t\t\t     {}\n".format(val[idx]))
            elif var == 'out_dir' or var == 'out_fn':
                glm_handle.write("\t{0} = '{1}'\n".format(var, val))
            else:
                glm_handle.write("\t{0} = {1}\n".format(var, val))
        
        glm_handle.write("/\n&init_profiles\n")
        for var in self.init_vars:
            val = value_dict[var]
            if type(val) == float or type(val) == int:
                glm_handle.write("\t{0} = {1}\n".format(var, val))
            elif type(val) == list:
                if type(val[0]) == str:
                    glm_handle.write("\t{0} =     {1}\n".format(var, val[0]))
                    for idx in np.array(range(len(val)-1))+1:
                        glm_handle.write("\t\t\t       {}\n".format(val[idx]))
                else:
                    glm_handle.write("\t{0} = {1},".format(var, val[0]))
                    for idx in np.array(range(len(val)-2))+1:
                        glm_handle.write("{},".format(val[idx]))
                    glm_handle.write("{}\n".format(val[-1]))
            else:
                glm_handle.write("\t{0} =  ".format(var))
                depths, wqs = val.shape
                for idx in range(wqs):
                    glm_handle.write("{0}, ".format(val[0][idx]))
                for row in np.array(range(depths-1))+1:
                    glm_handle.write("\n                    ")
                    for elmt in range(wqs):
                        glm_handle.write("{0}, ".format(val[row][elmt]))
       
        glm_handle.write("\n/\n&meteorology\n")
        for var in self.met_vars:
            val = value_dict[var]
            glm_handle.write("\t{0} = {1}\n".format(var, val))
        glm_handle.write("/\n&bird_model\n")
        for var in self.bird_vars:
            val = value_dict[var]
            glm_handle.write("\t{0} = {1}\n".format(var, val))
        glm_handle.write("/\n&outflows\n")
        for var in self.out_vars:
            val = value_dict[var]
            glm_handle.write("\t{0} = {1}\n".format(var, val))
        glm_handle.close()

        
        
        
        
        
        
        