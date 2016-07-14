# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 12:32:26 2016

@author: login
"""
import os
p1 = '/Users/login/Documents/GLM_Wrapper/glm_case_folders/testcase'
p2 = '/Users/login/Documents/GLM_Wrapper/GLM_Executables/examples_2.2/coldlake/fabm'

handles = [open(os.path.join(i,'glm2.nml'), 'r') for i in [p1, p2]]

for meL, youL in zip(handles[0], handles[1]):
    print meL, youL

closed = [handle.close() for handle in handles]    