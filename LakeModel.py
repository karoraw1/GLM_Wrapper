# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:03:39 2016
ncdump() was created by this dude:
http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html#code
JD_converter was created by this dude: (and modified by me)
https://asimpleweblog.wordpress.com/2010/06/20/julian-date-calculator/

@author: login
"""
import peakutils
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import JD_converter as jd

def printNaNWarning(df, label, fill_method=None):
    print "\nNaN Values read into {0} ({1})".format(label, str(type(df)))
    print df.isnull().sum()
    return None
    
def pandasfillwrap(seriesORdf, fill_method):
    if fill_method=="ZERO":
        print "\tSetting them equal to zero"
        seriesORdf.fillna(value=0, inplace=True)
    elif fill_method in ['backfill', 'bfill', 'pad', 'ffill']:
        seriesORdf.fillna(method=fill_method, inplace=True)
        print "\t using Pandas {} method".format(fill_method)
    else:
        print "\t FILL METHOD NOT RECOGNIZED"
        print "\tLeaving them alone"
    return None

# TODO: def standardscaling():
# TODO: def makeAggregations():
# TODO: def plotAggregations():
# TODO: redo this to apply to aggregations 

def printAutocorr(df, threshold=None):
    print "Autocorrelation peaks"    
    data_cols = {}
    ac_df = pd.DataFrame(index=df.index, columns=df.columns)
    ignorable = ['lat', 'lon', 'time_J', 'altitude']
    for col in df.columns:
        if df[col].dtype != '<M8[ns]' or col not in ignorable:
            print "\t{}:".format(col)
            autocorrs = np.array([df[col].autocorr(i) for i in range(1, len(df[col]))])
            data_cols[col] = autocorrs
            peaks = peakutils.indexes(autocorrs, thres=0.5, min_dist=100)
            #peakvals = {autocorrs[j]:j for j in peaks}
            print "\t\t{} peaks detected".format(len(peaks))
            print "\t\tTallest peak ({0}) @ lag {1}".format(autocorrs.max(),
                                                            autocorrs.argmax())
        else:
            ac_df = ac_df.drop(col, 1)
            
    for key in data_cols.keys():
            ac_df[key] = data_cols[key]
            
    return ac_df


# TODO: def plotfrequencydomain():



def insertTimeColumns(df, dates=None, insertAt=0):
    """this inserts day of the year, month of the year, and season of the year
    columns into a dataframe for use by the groupby function"""
    dates = pd.date_range(dates[0], periods=len(dates))
    month_i = dates.month
    day_i = dates.dayofyear
    year_i = dates.year
    fall, winter, spring, summer = [10,11,12], [1,2,3], [4,5,6], [7,8,9]
    season_i = np.zeros(len(month_i))
    for idx in range(len(month_i)):
        if month_i[idx] in fall:
            season_i[idx] = 1
        elif month_i[idx] in winter:
            season_i[idx] = 2
        elif month_i[idx] in spring:
            season_i[idx] = 3
        elif month_i[idx] in summer:
            season_i[idx] = 4
        else:
            print "unexpected error"

    name = ["year_i", "month_i", "season_i", "day_i"]
    new_cols = [year_i, month_i, season_i, day_i]
    for ind, name, col in zip(range(insertAt,4), name, new_cols):
        df.insert(ind, name, col)
    return df

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print "NetCDF Global Attributes:"
        for nc_attr in nc_attrs:
            print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print "NetCDF dimension information:"
        for dim in nc_dims:
            print "\tName:", dim 
            print "\t\tsize:", len(nc_fid.dimensions[dim])
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print "NetCDF variable information:"
        for var in nc_vars:
           if var not in nc_dims:
                print '\tName:', var
                print "\t\tdimensions:", nc_fid.variables[var].dimensions
                print "\t\tsize:", nc_fid.variables[var].size
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars

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
    output_vars = ["out_dir", "out_fn", "nsave", "csv_lake_fname", 
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
    inflow_vars = ["num_inflows", "names_of_strms", "subm_flag", 
                   "strm_hf_angle", "strmbd_slope", "strmbd_drag", 
                   "inflow_factor", "inflow_fl", "inflow_varnum", 
                   "inflow_vars", "coef_inf_entrain"]
    # &outflows: Information about outflows
    outflow_vars = ["num_outlet", "flt_off_sw", "outl_elvs", "bsn_len_outl", 
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
        
        for var in self.output_vars:
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
        for var in self.outflow_vars:
            val = value_dict[var]
            glm_handle.write("\t{0} = {1}\n".format(var, val))
        glm_handle.write("/\n&inflows\n")        
        for var in self.inflow_vars:
            val = value_dict[var]
            if type(val) != list:
                glm_handle.write("\t{0} = {1}\n".format(var, val))
            else:
                glm_handle.write("\t{0} = '{1}',\n".format(var, val[0]))
                for idx in np.array(range(len(val)-1))+1:
                    glm_handle.write("                  ")
                    glm_handle.write("'{0}',\n".format(val[idx]))    
        glm_handle.write("/")
        glm_handle.close()
        self.met_path = None
        self.raw_met = None
    
    def read_GCHN(self, met_path):
        self.met_path = met_path
        self.raw_met = pd.read_csv(met_path, parse_dates=[2], 
                                   infer_datetime_format=True, 
                                   na_values=[-9999])
        
        # TODO: play with dates in the way that we did in 
        # WaterBalance.py to get different resamples & groupbys        
        dates = self.raw_met.DATE
        
        #PRCP - Precipitation (tenths of mm)        
        self.precip = self.raw_met.PRCP/10000
        #AWND - Average daily wind speed (tenths of meters per second)        
        self.wind_speed = self.raw_met.AWND/10
        #TMAX - Maximum temperature (tenths of degrees C)
        self.t_max = self.raw_met.TMAX/10
        self.t_min = self.raw_met.TMIN/10
        #SNOW - Snowfall (mm) -> meters
        self.snow_fall = self.raw_met.SNOW/1000
        #SNWD - Snow depth (mm) -> meters
        self.snow_depth = self.raw_met.SNWD/1000
        self.clean_cols = [self.precip, self.wind_speed, self.t_max, self.t_min, 
                      self.snow_fall, self.snow_depth]
        col_labels = ["precip", "wind_speed", "t_max", "t_min", "snow_fall",
                      "snow_depth"]
        self.net_GCHN = pd.DataFrame(index=self.raw_met.DATE, columns=col_labels)
        for idx in range(len(self.clean_cols)):
            self.net_GCHN.iloc[:, idx] = self.clean_cols[idx].values
        
        printNaNWarning(self.net_GCHN, "GCHN data", None)
        self.GCHN_ac = printAutocorr(self.net_GCHN)        
        
    def read_CERES_nc(self, ceres_path):
        rootgrp = Dataset(ceres_path, "r", format="NETCDF3_CLASSIC")
        nc_attrs, nc_dims, nc_vars = ncdump(rootgrp, False)
        # Temps in Kelvin        
        # Fluxes in u'Watts per square meter'
        # wind vectors in meters per second
        # humididity & cloud cover in ?
        toExtract = ["lon", "lat", "Time_of_observation", 
                     "Surface_skin_temperature", "Surface_wind___U_vector",
                     "Surface_wind___V_vector", 
                     "Column_averaged_relative_humidity", 
                     "Altitude_of_surface_above_sea_level",
                     "PSF_wtd_MOD04_cloud_fraction_land",
                     "CERES_net_LW_surface_flux___Model_B",
                     "CERES_net_SW_surface_flux___Model_B"]
        betterNames = ["lat", "lon", "time_J", "temp", "windU", "windV",
                       "humidity", "altitude", "cloud_frac", "LW_rad", "SW_rad"]
        self.ceres_ = {}
        for var, nam in zip(toExtract, betterNames):
            self.ceres_[nam] = pd.Series(rootgrp.variables[var][:])
        gDays = [jd.caldate(j) for j in self.ceres_["time_J"]]
        self.ceres_['date'] = pd.to_datetime(np.array(gDays))
        
        self.ceres_df = pd.DataFrame(index=self.ceres_['date'], 
                                     columns=betterNames)
        for key in self.ceres_.keys():
            if key == 'temp':
                self.ceres_df[key] = self.ceres_[key].values-273.15
            else:
                self.ceres_df[key] = self.ceres_[key].values

        printNaNWarning(self.ceres_df, "CERES data", None)
        self.ceres_ac = printAutocorr(self.ceres_df)
        
        ceres_daily = self.ceres_df.resample("D", how='mean')
        #pass dataframe to insert new index columns
        self.ceres_daily = insertTimeColumns(ceres_daily, ceres_daily.index)
        
#        annualsum = data.groupby(data['year']).agg(np.sum)
#        monthlysum = data.groupby(data['month']).agg(np.sum)
#        dailysum = data.groupby(data['day']).agg(np.sum)
#        seasonalsum = data.groupby(data['season']).agg(np.sum)
#        #ceres_ needs to be converted into a dataframe        
        rootgrp.close()
        

                                   

        
        
        
        
        
        
        