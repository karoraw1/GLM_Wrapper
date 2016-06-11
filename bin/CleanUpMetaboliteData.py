# -*- coding: utf-8 -*-
"""
Created on Mon May 16 15:47:12 2016



@author: login
"""

import pandas as pd
import numpy as np
import os, sys

def strip(text):
    try:
        return text.strip()
    except AttributeError:
        return text


data_fn = "MysticLakeIons_Clean.csv"
data2_fn = "MysticLakeIons_Part2_Clean.csv"
data_path = os.path.join("..","waterData", data_fn)
data2_path = os.path.join("..","waterData", data2_fn)

ions_df = pd.read_csv(data_path, parse_dates=["Collection Date"],
                      converters = {'Sample Name' : strip})
ions2_df = pd.read_csv(data2_path, parse_dates=["Collection_Date"],
                      converters = {'Sample Name' : strip})

matchDepthandDate = {}
counter = 0
for i in xrange(ions2_df.shape[0]):
    new_meas = ions2_df.iloc[i,:]
    test1=False
    test2=False
    for j in xrange(ions_df.shape[0]):    
        old_meas = ions_df.iloc[j, :]
        test1=np.equal(new_meas.Collection_Date, old_meas["Collection Date"])
        try:
            o_m = int(old_meas["Depth"])
            n_m = int(new_meas["Depth"])
            test2=np.equal(n_m, o_m)
            
            if (test1 and test2) == True:
                matchDepthandDate[counter] = (new_meas.Collection_Date,
                                              new_meas["Depth"], 
                                              old_meas["Collection Date"],
                                              old_meas["Depth"])
                counter+=1

        except ValueError:
            o_m = str(old_meas["Depth"])
            n_m = str(new_meas["Depth"])
            test2= n_m == o_m
            
            if (test1 and test2) == True:
                matchDepthandDate[counter] = (new_meas.Collection_Date,
                                              new_meas["Depth"], 
                                              old_meas["Collection Date"],
                                              old_meas["Depth"])
                counter+=1

