#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 20:13:04 2017

@author: login
"""
import os
import pandas as pd 
import numpy as np
from collections import OrderedDict

inflow_fn = 'inflow_base.csv'
inflow_dn = '/Users/login/Documents/GLM_Wrapper/glm_case_folders/aed_manual'
inflow_path = os.path.join(inflow_dn, inflow_fn)
# Units 
aed_inflow_params = OrderedDict()
aed_inflow_params['OXY_oxy'] = 225
aed_inflow_params['SIL_rsi'] = 12.5
aed_inflow_params['NIT_amm'] = 12.5
aed_inflow_params['NIT_nit'] = 27.6
aed_inflow_params['PHS_frp'] = 0.25
aed_inflow_params['OGM_don'] = 70
aed_inflow_params['OGM_pon'] = 5
aed_inflow_params['OGM_dop'] = 10
aed_inflow_params['OGM_pop'] = 10
aed_inflow_params['OGM_doc'] = 80
aed_inflow_params['OGM_poc'] = 700
aed_inflow_params['PHY_green'] = 40
aed_inflow_params['PHY_crypto'] = 40
aed_inflow_params['PHY_diatom'] = 40 

inflow_orig = pd.read_csv(inflow_path,index_col=0)
inflow_mod = inflow_orig.copy()

for name, value in aed_inflow_params.items():
    print name
    inflow_mod[name] = np.ones((inflow_orig.shape[0],))*value

inflow_mod_fn = 'inflow.csv'
inflow_mod_path = os.path.join(inflow_dn, inflow_mod_fn)
inflow_mod.to_csv(inflow_mod_path, index_label='time')



















    

