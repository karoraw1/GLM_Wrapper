#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 11:38:23 2017

@author: login
"""

import pandas as pd
import numpy as np

def abund_vec_to_list(sp_ser):
    sp_list = []
    for idx in sp_ser.index:
        sp_count = sp_ser.ix[idx]
        for cnt in range(int(sp_count)):
            sp_list.append(idx)
    return np.array(sp_list)

def rarefy_table(otu_table):
    """
    Rarefies and otu table (pandas dataframe) to minimum lib size
    """
    lib_size = int(otu_table.sum(axis=1).min())
    print "Minimum library size is {}".format(lib_size)
    data_ = np.zeros(otu_table.shape)
    rarefiedTable = pd.DataFrame(data = data_, columns=otu_table.columns,
                                 index=otu_table.index)
    for row in otu_table.index:
        def add_cnt_to_col(c, n):
            rarefiedTable.ix[row, c] = n
            return None
        sp_ser = otu_table.ix[row, :]
        sp_pool = abund_vec_to_list(sp_ser)
        shrunk_pool = np.random.choice(sp_pool, size=(1, lib_size), replace=False)
        cols, cnts = np.unique(shrunk_pool, return_counts=True)
        map(add_cnt_to_col, cols, cnts)

    return rarefiedTable


        

        

    