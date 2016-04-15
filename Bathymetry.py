# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 18:40:15 2016
Bathymetry measurements were done using ImageJ. The original image was edited
to remove text and make contiguous contours. Additional images were created 
for each contour with all inner contours removed, so as to eliminate the borders
from being included in the area measurements. These were called Layers1-9
and are in the BathymetryLayers folder. An additional final layer was made of
the entire lake to compare with recorded values of surface area. For our
purposes we are only going to simulate the main body of the lake, however
the beach and inlet areas would also be included in recorded measurements of
surface area. 

The scaling was 0.5848 pixels/meter.



@author: login
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D


area = np.array([77373.80, 148475.73, 202472.97, 257818.95, 338552.69, 
                 397077.50, 460778.04, 524802.66, 560051.22])
                 
elevation =  np.array([-23.384, -20.336, -17.288, -14.24, -11.192, 
                       -8.144, -5.096, -2.048, 1])
surf = elevation.max()
floor = elevation.min()
intervals = np.round(surf-floor, decimals=0)
elevation_int = np.linspace(floor, surf, intervals)
radius = np.sqrt(area/np.pi)
radii_int = np.interp(elevation_int, elevation, radius)

t = np.linspace(0,2*np.pi, num=len(radii_int));
fig = plt.figure()
ax = fig.gca(projection='3d')   
for r, z in zip(radii_int, elevation_int):
    ax.plot(r*np.cos(t), r*np.sin(t), z)
plt.show()

