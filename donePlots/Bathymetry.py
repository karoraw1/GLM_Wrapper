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

area = np.array([ 77373.80, 148475.73, 202472.97, 
                 257818.95, 338552.69, 397077.50, 
                 460778.04, 524802.66, 560051.22])
             
elevation =  np.array([-23.384, -20.336, -17.288, 
                       -14.24, -11.192, -8.144, 
                       -5.096, -2.048, 1.00])       


volume = np.zeros(len(area))
layers = range(0,len(area))
alphas, betas = np.zeros(len(area)), np.zeros(len(elevation))

volume[0] = area[0] * (elevation[1]-elevation[0])/3
for i in range(1,len(area)):
    upper_plus_average = area[i-1] + 0.5*(area[i]-area[i-1])
    thickness = (elevation[i]-elevation[i-1])
    volume[i] = volume[i-1] + upper_plus_average*thickness

alphas[0] = (np.log10(volume[1]/volume[0]) / np.log10(elevation[1]/elevation[0])) 
betas[0] = (np.log10(area[1]/area[0]) / np.log10(elevation[1]/elevation[0]))

for i in range(1,len(area)):
    denominator = np.log10(elevation[i-1]/elevation[i])
    alphas[i] = (np.log10(volume[i-1]/volume[i]) / denominator ) 
    betas[i] = (np.log10(area[i-1]/area[i]) / denominator )


n_layers = range(0,len(area))
meters = np.round(elevation[-1]-elevation[0], decimals=0)
refined_depths = np.linspace(elevation[0], elevation[-1], meters*10)
refined_vols = []
refined_areas = []

for a, b, h_b, v_b, i in zip(alphas, betas, elevation, volume, range(len(betas))):
    for h in refined_depths:
        if h > h_b and h < 
            refined_vols.append()
    



"""
radii_int = np.interp(elevation_int, elevation, radius)

t = np.linspace(0,2*np.pi, num=len(radii_int));
fig = plt.figure()
ax = fig.gca(projection='3d')   
for r, z in zip(radius, elevation):
    ax.plot(r*np.cos(t), r*np.sin(t), z)
plt.show()

plt.xlabel("meters")
"""
