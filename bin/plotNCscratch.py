# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 13:31:30 2016

@author: login
"""

import os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

pwd = os.path.dirname(os.getcwd())
out_p = os.path.join(pwd,'GLM_Executables','examples_2.2','warmlake','fabm')
out_f = os.path.join(out_p, 'output.nc')
rootgrp = Dataset(out_f, "r")
temp = rootgrp['temp'][:,:50,0,0].data
temp[temp > 1000] = -1
plt.imshow(np.flipud(temp[:200,:].T))
plt.colorbar()

