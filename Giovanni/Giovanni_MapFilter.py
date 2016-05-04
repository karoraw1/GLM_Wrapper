# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:41:42 2016

@author: login
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import h5py
from netCDF4 import Dataset, num2date
import time, calendar, datetime
from mpl_toolkits.basemap import Basemap, cm


test_f = "OMI-Aura_L2-OMMYDCLD_2014m0101t1651-o50349_v003-2015m0228t154805.he5"
test_time = "2014-01-01 16:51:19"
path2obj = "/HDFEOS/SWATHS/CloudProduct/Data Fields/CloudFractionforO2"

f = h5py.File(test_f)
group = f['HDFEOS']['SWATHS']['CloudProduct']
timeIndex = group['nTimes'].value

geoLoc = group['Geolocation Fields']
dFields = group['Data Fields']
time_2 = geoLoc['Time'].value
latIndex = geoLoc['Latitude'].value
lonIndex = geoLoc['Longitude'].value
cloudFrac = dFields['CloudFractionforO2']
naVal = cloudFrac.fillvalue
mask = cloudFrac.value != naVal

true_CF = cloudFrac[mask]
true_lat = latIndex[mask]
true_lon = lonIndex[mask]

lake = np.array([42.4317, -71.1483])
lake_lat = np.array([lake[0]-0.5, lake[0]+0.5])
lake_lon = np.array([lake[1]-0.5, lake[1]+0.5])
latcorners = np.array([lake[0]-1, lake[0]+1])
loncorners = np.array([lake[1]-1, lake[1]+1])
center_lat = latcorners.mean()
center_lon = loncorners.mean()

# create polar stereographic Basemap instance.

fig = plt.figure(figsize=(24,12))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
m = Basemap(llcrnrlon=loncorners[0], llcrnrlat=latcorners[0],
            urcrnrlon=loncorners[1], urcrnrlat=latcorners[1], \
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc',\
            lat_0=center_lat, lon_0=center_lon)
m.drawcoastlines()
#m.fillcontinents()
# draw parallels
m.drawparallels(np.array([42., 43.,]),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.array([-72.,-71.]),labels=[1,1,0,1])
ax.set_title('Great Circle from New York to London')
x, y = m(true_lat,true_lon)
m.drawmapboundary(fill_color='#99ffff')
#m.fillcontinents(color='#cc9966',lake_color='#99ffff')
m.scatter(x,y,3,marker='o',color='k')
plt.show()
