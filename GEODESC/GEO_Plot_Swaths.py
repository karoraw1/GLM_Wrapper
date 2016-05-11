# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:41:42 2016

@author: login
"""
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import pandas as pd
import numpy as np
import h5py
from netCDF4 import Dataset, num2date
import time, calendar, datetime
from mpl_toolkits.basemap import Basemap, cm


test_f = "test.he5"
test_time = "2014-01-01 16:51:19"

f = h5py.File(test_f)
group = f['HDFEOS']['SWATHS']['CloudProduct']
timeIndex = group['nTimes'].value
geoLoc = group['Geolocation Fields']
dFields = group['Data Fields']
time_2 = geoLoc['Time'].value

sliced = False
if sliced == True:
    slicer_c = abs(np.arange(0, 65, 5)-1)
    slicer_r = abs(np.arange(0, 1644, 5)-1)
    latIndex = geoLoc['Latitude'].value[:, slicer_c]
    lonIndex = geoLoc['Longitude'].value[:, slicer_c]
    cloudFrac = dFields['CloudFractionforO2'].value[:, slicer_c]
    latIndex = latIndex[slicer_r, :]
    lonIndex = lonIndex[slicer_r, :]
    cloudFrac = cloudFrac[slicer_r, :]
else:
    latIndex = geoLoc['Latitude'].value
    lonIndex = geoLoc['Longitude'].value
    cloudFrac = dFields['CloudFractionforO2'].value
    
naVal = dFields['CloudFractionforO2'].fillvalue
mask = cloudFrac != naVal

#
masking = False
if masking == True:
    true_CF = cloudFrac[mask]
    true_lat = latIndex[mask]
    true_lon = lonIndex[mask]
else:
    true_CF = cloudFrac
    true_lat = latIndex
    true_lon = lonIndex
    
lake = np.array([42.4317, -71.1483])

filter_Nan = cloudFrac != naVal

deg_lim = np.arange(0, 2, 0.1)
hits = np.zeros(len(deg_lim))
for idx, m in enumerate(deg_lim):
    lat_fil = (abs(latIndex - lake[0]) < m) 
    lon_fil = (abs(lonIndex - lake[1]) < m)
    hits[idx] = ((lat_fil & lon_fil & filter_Nan).sum())

plot_filter = False

if plot_filter == True:
    plt.scatter(deg_lim, hits, s=60, label="Adjacent Measurements")
    plt.scatter(deg_lim[2], hits[2], s=100, c='r', label= "Non-Zero Minimum")
    plt.xlim([0,2])
    plt.ylim([-10, 500])
    plt.xlabel("Degrees distance (Lat & Long)")
    plt.ylabel("Observations (n)")
    plb.rc('font', family='serif', size=18)
    
    if (len(hits.nonzero()[0]) != 0):
        min_lim = deg_lim[hits.nonzero()[0][0]]
    else:
        min_lim = 1.0
else:
    pass


lake_lat = np.array([lake[0]-0.5, lake[0]+0.5])
lake_lon = np.array([lake[1]-0.5, lake[1]+0.5])
latcorners = np.array([lake[0]-1, lake[0]+1])
loncorners = np.array([lake[1]-1, lake[1]+1])
center_lat = latcorners.mean()
center_lon = loncorners.mean()

# create polar stereographic Basemap instance.

fig = plt.figure(figsize=(24,12))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
map = Basemap(llcrnrlon=-100, llcrnrlat=-80,
            urcrnrlon=80, urcrnrlat=80, \
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc',\
            lat_0=center_lat, lon_0=center_lon)
map.drawcoastlines(color='w')
#map.fillcontinents()
# draw parallels
map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,1])
map.drawparallels(np.arange(-90,90,30),labels=[1,1,0,1])
#map.drawparallels(np.array([42., 43.,]),labels=[1,1,0,1])
# draw meridians
#m.drawmeridians(np.array([-72.,-71.]),labels=[1,1,0,1])
ax.set_title('OMI Cloud Fraction Swath')
x, y = map(true_lat,true_lon)
map.drawmapboundary(fill_color='#99ffff')
map.scatter(x,y,marker='o',color='k', s=true_CF*4)
map.fillcontinents(color='#cc9966',lake_color='#99ffff', alpha=0.5)
plt.show()
