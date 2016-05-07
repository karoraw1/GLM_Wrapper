# -*- coding: utf-8 -*-
"""
Created on Wed May  4 14:36:30 2016

@author: Keith
"""
import subprocess
import numpy as np
import h5py
import os
#import time

def dl_parse_hdf(data_urls, data_vals, i, build = False):
    d_u = data_urls[i]
    if (d_u[0:3] == 'ftp') and (d_u[:-4:-1] == '5eh'):
        #start = time.time()
        wget = subprocess.Popen("wget -nv " + d_u, shell = True, executable = '/bin/bash')
        wget.communicate()
        hdf_f = os.path.basename(d_u)
        #print "Time to download: ", (time.time()-start)
        try:
            f = h5py.File(hdf_f)
        except IOError:
            hdf_f = hdf_f + ".1"
            f = h5py.File(hdf_f)
        
        try: 
            group = f['HDFEOS']['SWATHS']['CloudProduct']
        except KeyError:
            if hdf_f[-1] == '1':
                hdf_f = hdf_f[:-2]
            else:
                hdf_f = hdf_f + ".1"
            f = h5py.File(hdf_f)
            group = f['HDFEOS']['SWATHS']['CloudProduct']
            
        geoLoc = group['Geolocation Fields']
        dFields = group['Data Fields']
        latIndex = geoLoc['Latitude'].value
        lonIndex = geoLoc['Longitude'].value
        cloudFrac = dFields['CloudFractionforO2'].value
        naVal = dFields['CloudFractionforO2'].fillvalue
        lake = np.array([42.4317, -71.1483])
        filter_Nan = cloudFrac != naVal

        deg_lim = np.arange(0, 1, 0.1)
        hits = np.zeros(len(deg_lim))
        for idx, m in enumerate(deg_lim):
            lat_fil = (abs(latIndex - lake[0]) < m) 
            lon_fil = (abs(lonIndex - lake[1]) < m)
            hits[idx] = ((lat_fil & lon_fil & filter_Nan).sum())
        
        if (len(hits.nonzero()[0]) != 0):
            min_lim = deg_lim[hits.nonzero()[0][0]]
        else:
            min_lim = 1.0
        
        filter_pos = (abs(latIndex - lake[0]) < min_lim) & (abs(lonIndex - lake[1]) < min_lim)

        filter_Nan = cloudFrac != naVal
        full_mask = filter_pos & filter_Nan
        
        good_recs = cloudFrac[full_mask].copy()
        good_lats = latIndex[full_mask].copy()
        good_lons = lonIndex[full_mask].copy()
        
        if build == False:
            subprocess.call('rm '+hdf_f, shell = True, executable = '/bin/bash')
        
        if (good_recs.shape[0]!=0):
            return [i, good_recs.mean(), good_recs.std(), good_recs.shape[0], 
                    good_lats.std(), good_lons.std(), hdf_f[-4:], min_lim]
        else:
            return [i, np.nan, np.nan, 0, np.nan, np.nan, hdf_f[-4:], min_lim]
    else:
        print "not a hdf file"
        return [i, np.nan, np.nan, np.nan, np.nan, np.nan, hdf_f[-4:], np.nan]
        
