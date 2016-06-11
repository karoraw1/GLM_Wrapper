# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:41:42 2016

This downloads xml metadata files from a text list 
It reads the start, end, and corresponding data file name.
It writes this data out to a CSV
If passed a command line argument, a test will be performed on the first
16 files

@author: Keith
"""

import collections, math, random, time, zipfile
from six.moves import urllib
from six.moves import xrange  # pylint: disable=redefined-builtin
import subprocess
import pandas as pd
import numpy as np
import os, sys
import mmap


    ## Turn list of lists into individual array variables
    i_out, cc_means, cc_stds, cc_n, lat_stds, lons_stds, checksum, min_lim = np.array(zip(*data_vals))
    new_cols = [i_out, cc_means, cc_stds, cc_n, lat_stds, lons_stds, checksum, min_lim]
    new_labels = ["i_out_2", "cc_means", "cc_stds", "cc_n", "lat_stds", 
                  "lons_stds", "checksum", "min_lim"]

    for colu, labe in zip(new_cols, new_labels):
        clean_clouds[labe] = colu

    cleanup = clean_clouds.checksum != '.he5'
    
    if (cleanup.sum() != 0) and build == False:
        cleanuplist = list(clean_clouds.hdf5[cleanup].values)
        for leftover in cleanuplist:
            os.remove(leftover)
    
    if clean_clouds.columns[0][:7] == 'Unnamed':
        clean_clouds.drop(clean_clouds.columns[0], axis=1)
        
    
    if test == False:
        d_csv_p = os.path.join(giovanni_path, "cloud_data_first_build.csv")
        print "cloud fetch complete"
    else:
        d_csv_p = os.path.join(giovanni_path, "build_test_data.csv")
        print "testing complete"
    
    clean_clouds.to_csv(d_csv_p)
    
    
    
