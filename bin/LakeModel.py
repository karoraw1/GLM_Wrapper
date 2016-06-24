# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:03:39 2016
ncdump() was created by this dude:
http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html#code
JD_converter was created by this dude: (and modified by me)
https://asimpleweblog.wordpress.com/2010/06/20/julian-date-calculator/

@author: login
"""
from scipy.stats import pearsonr
import re
import mmap
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import JD_converter as jd

def import_data(return_type):
    
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
                
    expectations = { "setup": setup_vals, "wq": wq_vals, 'morpho': morpho_vals, 
                    'time': time_vals, 'output': output_vals, 
                    'init': init_vals, 'met': met_vals, 'bird': bird_vals, 
                    'outflow': outflow_vals, 'inflow': inflow_vals}
                    
    parameters =  { "setup": setup_vars, "wq": wq_vars, 'morpho': morpho_vars,
                   'time': time_vars, 'output': output_vars, 'init': init_vars,
                   'met': met_vars, 'bird': bird_vars, 'outflow': outflow_vars,
                   'inflow': inflow_vars }
                   
    if return_type == 'Mystic_Vals':
        return expectations
    elif return_type == "Parameters":
        return parameters
    elif return_type == "Custom":
        print "Custom Return Type Selected"

def error_metrics(Obs, Sim):
    Err = Sim - Obs
    MSE = np.mean(Err**2)
    rMSE = np.sqrt(MSE)
    NSE = 1 - (MSE/np.var(Obs))
    r = pearsonr(Obs, Sim)
    a = np.std(Sim) / np.std(Obs)
    b = np.mean(Sim) / np.mean(Obs)
    # Kling-Gupta Efficiency
    KGE = 1- np.sqrt((r[0]-1)**2 + (a-1)**2 + (b-1)**2)
    return {'MSE':MSE, 'rMSE':rMSE, 'NSE':NSE, 'r':r, 'a':a, 'b':b, 'KGE':KGE}

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

def printNaNWarning(df, label, fill_method=None):
    print "\nNaN Values read into {0} ({1})".format(label, str(type(df)))
    print df.isnull().sum()
    print "Total measurements per variable: {}".format(df.shape[0])
    return None

def plotTemporalDist(df, fignum, clear=True, bins=None):
    days = df.index.dayofyear
    if bins == None:
        buckets = len(np.unique(days))
    else:
        buckets = bins
        
    plt.figure(fignum)
    
    if clear == True:
        plt.clf()
        
    n, bins, patches = plt.hist(days, buckets, facecolor='green', alpha=0.75)
    
    return n, bins, patches

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

def z_score(df):
    z_df = pd.DataFrame(index=df.index, columns=df.columns)    
    for col in df.columns:
        z_df[col] = (df[col] - np.mean(df[col])) / np.std(df[col], ddof=1)
    return z_df
    
def makeAggregations(df, indices, fxn):
    aggs = {}    
    for idx in indices:
        new_df = df.copy().groupby(df[idx]).agg(fxn)
        aggs[idx] = new_df.drop(indices, axis=1)
    
    return aggs
    
def show_agg_resolution(aggs, col_name):
    day_frame = pd.DataFrame(index=aggs['day_i'].index, columns=aggs.keys())
    for col in day_frame.columns:
        ser = aggs[col][col_name]
        if col == 'day_i':
            new_idx = ser.index
        elif col == 'week_i':
            new_idx = ser.index*7+2
        elif col == 'month_i':
            new_idx = ser.index*30+6
        elif col == 'season_i':
            new_idx = ser.index*(366/4)+2
        elif col == 'year_i':
            new_idx = pd.Index([366], dtype='Int64')
            
        if col == 'year_i':
            mod_ser = pd.Series(data=ser.values.mean(), index=new_idx)
            day_frame[col] = mod_ser.reindex(day_frame.index)
        else:
            mod_ser = pd.Series(data=ser.values, index=new_idx)
            day_frame[col] = mod_ser.reindex(day_frame.index)
        
    return day_frame


# TODO: def plotAggregations():
# TODO: redo this to apply to aggregations 

def printAutocorr(aggs, threshold=None):
    print "Autocorrelation peaks"    
#    data_cols = {}
#    ac_df = pd.DataFrame(index=df.index, columns=df.columns)
#    
#    for col in df.columns:
#        if df[col].dtype != '<M8[ns]':
#            print "\t{}:".format(col)
#            autocorrs = np.array([df[col].autocorr(i) for i in range(1, len(df[col]))])
#            data_cols[col] = autocorrs
#            peaks = peakutils.indexes(autocorrs, thres=0.5, min_dist=100)
#            #peakvals = {autocorrs[j]:j for j in peaks}
#            print "\t\t{} peaks detected".format(len(peaks))
#            print "\t\tTallest peak ({0}) @ lag {1}".format(autocorrs.max(),
#                                                            autocorrs.argmax())
#            
#    for key in data_cols.keys():
#            ac_df[key] = data_cols[key]
    for key in aggs.keys():
        df = aggs[key].dropna()
        print key
        for col in df.columns:
            if col[-1] != "i":
                autocorrs = np.array([df[col].autocorr(i) for i in range(1, len(df[col])/2)])
                print "\t{0} max autocorr @ {1} = {2:.2f}".format(col, np.argmax(autocorrs)+1, autocorrs.mean())
            
    return None


# TODO: def plotfrequencydomain():



def TimeIdx(df, dates=None, insertAt=0):
    """this inserts day of the year, month of the year, and season of the year
    columns into a dataframe for use by the groupby function"""

    if dates == None:
        dates = df.index
    else:
        pass
    
    dates = pd.date_range(dates[0], periods=len(dates))
    month_i = dates.month
    week_i = dates.week
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

    i_names = ["year_i", "month_i", "season_i", 'week_i', "day_i"]
    new_cols = [year_i, month_i, season_i, week_i, day_i]
    for ind, name, col in zip(range(insertAt,len(new_cols)), i_names, new_cols):
        df.insert(ind, name, col)
    return df, i_names

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
    
    def __init__(self, name, dir_path):
        
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
        self.parameters = import_data("Parameters")

    def fill_parameters_set(self, val_set = "Mystic_Vals"):
        
        self.expectations = import_data(val_set)
        self.expectations['output'][1] = self.expectations['output'][1] + "_"+self.name
        self.expectations['output'][0] = self.dir_path+"/output"        

        safe_dir(self.expectations['output'][0])
        # start writing config file
        glm_handle = open(self.glm_path, 'w+')
        
        glm_handle.write("&glm_setup\n")        
        for var,val in zip(self.parameters["setup"],self.expectations['setup']):
            glm_handle.write("\t{0} = {1}\n".format(var, val))

        glm_handle.write("/\n&wq_setup\n")
        for var,val in zip(self.parameters["wq"],self.expectations['wq']):
            glm_handle.write("\t{0} = {1}\n".format(var, val))

        glm_handle.write("/\n&morphometry\n")
        for var, val in zip(self.parameters['morpho'],self.expectations['morpho']):
            glm_handle.write("\t{0} = {1}\n".format(var, val))

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
        for var,val in zip(self.parameters["time"],self.expectations['time']):
            glm_handle.write("\t{0} = {1}\n".format(var, val))
            
        glm_handle.write("/\n&output\n")           
        for var, val in zip(self.parameters["output"], self.expectations["output"]):
            if type(val) == list:
                glm_handle.write("\t{0} = {1}\n".format(var, val[0]))
                for idx in np.array(range(len(val)-1))+1:
                    glm_handle.write("\t\t\t\t     {}\n".format(val[idx]))
            elif var == 'out_dir' or var == 'out_fn':
                glm_handle.write("\t{0} = '{1}'\n".format(var, val))
            else:
                glm_handle.write("\t{0} = {1}\n".format(var, val))
        
        glm_handle.write("/\n&init_profiles\n")
        for var, val in zip(self.parameters["init"],self.expectations["init"]):
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
        for var, val in zip(self.parameters["met"],self.expectations["met"]):
            glm_handle.write("\t{0} = {1}\n".format(var, val))

        glm_handle.write("/\n&bird_model\n")
        for var, val in zip(self.parameters["bird"],self.expectations["bird"]):
            glm_handle.write("\t{0} = {1}\n".format(var, val))
        
        glm_handle.write("/\n&outflows\n")
        for var, val in zip(self.parameters["outflow"],self.expectations["outflow"]):
            glm_handle.write("\t{0} = {1}\n".format(var, val))
            
        glm_handle.write("/\n&inflows\n")        
        for var, val in zip(self.parameters["inflow"],self.expectations["inflow"]):
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
        
    def read_GEODISC_cloud_data(self, fname=None, data_dir=None):
        if fname == None:
            self.geodisc_f = 'cloud_data.csv'
        else:
            self.geodisc_f = fname
        if data_dir == None:
            self.geodisc_d = os.path.join(os.getcwd(), "GEODESC")
        else:
            self.geodisc_d = data_dir

        self.geodisc_p = os.path.join(self.geodisc_d, self.geodisc_f)        
        
        try:
            self.geodisc_clouds = pd.read_csv(self.geodisc_p, index_col = 0,
                                              parse_dates = ['start', 'end'],
                                              infer_datetime_format = True)
        except IOError:
            print "no file detected"
            raise IOError
        
class USGS_water_data(object):

    def __init__(self, fn):
        self.fn = fn
    
    def preprocess(self):
        f = open(self.fn)    
        s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        self.alt1 = False  
        
        try:
            header = s.find("#    DD parameter   Description")
            s.seek(header)
        except ValueError:
            header = s.find("#    DD parameter statistic")
            s.seek(header)
            self.alt1=True            

        s.seek(header)
        s.readline()
        text = s.readline()
        words = [splits for splits in text.split(" ") if splits is not ""]
        dds, params, descrips  = [], [], []
        
        while words[0] != "#\n":
            dds.append(words[1])
            
            if self.alt1 == False:
                params.append(words[2])
                descrips.append(" ".join(words[3:-1])[:-1])
            else:
                params.append(words[2]+"_"+words[3])
                descrips.append(" ".join(words[4:-1])[:-1])
                
            text = s.readline()
            words = [splits for splits in text.split(" ") if splits is not ""]
        s.seek(0)        
        line = s.readline()
        com_lines = 0
        while line[0] == '#':
            line = s.readline()
            com_lines +=1
            
        self.cols = line
        re.split(r'\t+', self.cols)
        self.skipRows = []
        self.com_lines, counter2 = com_lines, com_lines
        
        while line[:4]!='USGS':
            line = s.readline()
            counter2+=1
            self.skipRows.append(counter2)
            
        s.close()
        f.close()
        
        col_keys = [i+"_"+j for i, j in zip(dds, params)]            
        self.col_dict = {i:j for i, j in zip(col_keys, descrips)}
        self.qual_cols = [i+"_cd" for i in col_keys]
        
    def read(self):
        self.df = pd.read_csv(self.fn, sep="\t", header=self.com_lines, 
                                  skiprows=self.skipRows, index_col=[2], 
                                  parse_dates=[2], low_memory=False,
                                  infer_datetime_format=True)
        
        new_cols = []
        for vals in self.df.columns:
            if self.col_dict.has_key(vals):
                new_cols.append(self.col_dict[vals])
            else:
                new_cols.append(vals)
        self.df.columns = new_cols

    def print_metadata(self):
        print "\n Metadata for %r" %os.path.basename(self.fn)
        
        for i in self.df.columns:
            if i[-3:] == '_cd':
               qual_col = self.df[i]
               print i
               print "\t%r" % qual_col.value_counts().index
               print "\t%r" % qual_col.value_counts().values
        print ""

class GHCN_weather_data(object):
    
    def __init__(self, met_path):
        self.met_path = met_path
        self.raw_met = pd.read_csv(met_path, parse_dates=[2], 
                                   infer_datetime_format=True, 
                                   na_values=[-9999])

    def clean_columns(self):
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
        self.df = pd.DataFrame(index=self.raw_met.DATE, columns=col_labels)

        for idx in range(len(self.clean_cols)):
            self.df.iloc[:, idx] = self.clean_cols[idx].values
        
        self.gchn_zs = z_score(self.df)
        printNaNWarning(self.df, "GCHN data", None)

                                           
class CERES_nc(object):
    
    def __init__(self, ceres_path, type):
        if type == "3":
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
    
            ignorable = ['lat', 'lon', 'time_J', 'altitude']
    
            for key in self.ceres_.keys():
                if key == 'temp':
                    self.ceres_df[key] = self.ceres_[key].values-273.15
                elif key in ignorable:
                    self.ceres_df = self.ceres_df.drop(key, 1)
                elif key not in self.ceres_df.columns.values:
                    pass
                else:
                    self.ceres_df[key] = self.ceres_[key].values
            
            self.ceres_zs = z_score(self.ceres_df)
            self.ceres_d_zs = self.ceres_zs.resample("D", how='mean')
            self.ceres_d_zs = self.ceres_d_zs.ix[:-1]        
            printNaNWarning(self.ceres_df, "CERES data", None)
            self.ceres_reindex_zs, self.i_names = TimeIdx(self.ceres_d_zs)
            self.ceres_mean_aggs = makeAggregations(self.ceres_reindex_zs, 
                                                    self.i_names, np.mean)
                                                    
            self.humidity_scales = show_agg_resolution(self.ceres_mean_aggs, 
                                                       "humidity")
            self.cloud_scales = show_agg_resolution(self.ceres_mean_aggs, 
                                                       "cloud_frac")
            printAutocorr(self.ceres_mean_aggs)
            #pass dataframe to insert new index columns
            rootgrp.close()
        elif type == "4":
            import xarray as xr
            data = xr.open_dataset(ceres_path)
            vars = data.keys()
            self.lats_n = data['lat'].values
            self.lons_n = data['lon'].values
            self.time_n = data['time'].values
            other_vars = vars[5:]
            time_ns, lats_ns, lons_ns = data[vars[3]].values.shape
            
            row1 = data[vars[3]].values[:,0,0]
            row2 = data[vars[4]].values[:,0,0]
            stack = np.vstack((row1, row2))
            for var in other_vars:
                row3 = data[var].values[:,0,0]    
                stack = np.vstack((stack, row3))
            
            self.df1 = pd.DataFrame(index=self.time_n, columns=vars[3:], data = stack.T.copy())
                
            
            row1 = data[vars[3]].values[:,0,1]
            row2 = data[vars[4]].values[:,0,1]
            stack = np.vstack((row1, row2))
            for var in other_vars:
                row3 = data[var].values[:,0,1]
                stack = np.vstack((stack, row3))
            
            self.df2 = pd.DataFrame(index=self.time_n, columns=vars[3:], data = stack.T.copy())
