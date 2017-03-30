#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 14:50:14 2017

@author: login
"""

import os, sys
from collections import OrderedDict

def load_aed_var_types(var_sub_type):
    """ 
    Loads the names of the variables for the files in which they are unnamed
    """
    if var_sub_type == 'zoo':
        zoo_list = []
        return zoo_list
    elif var_sub_type == 'phyto':
        phyto_data = OrderedDict()

        phyto_basic = ['p_name', 'p_initial', 'p0', 'w_p', 'Ycc']
        phyto_data['Basic'] = phyto_basic

        phyto_grow = ['R_growth', 'fT_Method', 'theta_growth','T_std', 'T_opt', 
                      'T_max']
        phyto_data['Growth'] = phyto_grow

        phyto_light = ['lightModel', 'I_K', 'I_St', 'KePHY']
        phyto_data['Light'] = phyto_light

        phyto_resp = ['f_pr', 'R_resp', 'theta_resp', 'K_fres', 'K_fdom']
        phyto_data['Respiration'] = phyto_resp

        phyto_data['Salinity'] = ['salTol','S_Bep', 'maxSP', 'Sop']

        phyto_nitro = ['simDINUptake', 'simDONUptake', 'simNFixation', 
                       'simINDynamics', 'N_o', 'K_N', 'INcon', 'INmin', 'INmax', 
                       'UNmax', 'gthRedNFix', 'NFixationRate']
        phyto_data['Nitrogen'] = phyto_nitro

        pytyo_phospho = ['simDIPUptake', 'simIPDynamics', 'Po', 'KP', 'IPcon', 
                         'IPmin', 'IPmax', 'UPmax']
        phyto_data['Phosphorus'] = pytyo_phospho
        phyto_data['Silica'] = ['simSiUptake', 'Sio', 'KSi', 'Sicon']
        phyto_list = []
        for group in phyto_data.keys():
            phyto_list += phyto_data[group]
        return phyto_list
    elif var_sub_type == 'pathogen':
        pathogen_vars = ["p_name", "coef_grwth_uMAX", "coef_grwth_Tmin",
                         "coef_grwth_Tmax","coef_grwth_T1","coef_grwth_T2",
                         "coef_grwth_Kdoc","coef_grwth_ic","coef_mort_kd20",
                         "coef_mort_theta","coef_mort_c_SM","coef_mort_alpha",
                         "coef_mort_beta","coef_mort_c_PHM","coef_mort_K_PHM",
                         "coef_mort_delta_M","coef_mort_fdoc","coef_light_kb_vis",
                         "coef_light_kb_uva","coef_light_kb_uvb",
                         "coef_light_cSb_vis","coef_light_cSb_uva",
                         "coef_light_cSb_uvb","coef_light_kDOb_vis",
                         "coef_light_kDOb_uva","coef_light_kDOb_uvb",
                         "coef_light_cpHb_vis","coef_light_cpHb_uva",
                         "coef_light_cpHb_uvb","coef_light_KpHb_vis",
                         "coef_light_KpHb_uva","coef_light_KpHb_uvb",
                         "coef_light_delb_vis","coef_light_delb_uva",
                         "coef_light_delb_uvb","coef_pred_kp20",
                         "coef_pred_theta_P","coef_sett_fa","coef_sett_w_path"]
        return pathogen_vars
    else:
        sys.exit("Illegal variable requested")


def exclamationNotFirst(string):
    if string.strip()[0] != '!':
        return True
    else:
        return False

def ampersandFirst(string):
    if string[0] == '&':
        return True
    else:
        return False


def read_aed_nml(file_name=None):
    if not file_name:
        test_aed_fn = 'sample_aed2.nml'
        test_path = os.path.join(os.path.dirname(os.getcwd()), 'test_files')
        file_name = os.path.join(test_path, test_aed_fn)        
    # open file and read into list of lines
    with open(file_name, 'r') as f:
        content = f.read().split("\n")
    # remove all empty lines
    c_list_full = filter(None, content)
    # remove all comment lines
    params_only = filter(exclamationNotFirst, c_list_full)
    # preemtively get module names into a list
    block_names = filter(ampersandFirst, params_only)
    # set up table to store each module's variable names 
    param_dict = OrderedDict()
    # set up table, within larger table, to make values accesible by var name
    for block in block_names:
        param_dict[block[1:]] = OrderedDict()
    
    # Data assimilation list parse
    for l in params_only:
        l = l.strip()
        # Which submodule am i in? 
        if l[0] == '&':
            curr_block = l[1:]
        # While variable am I dealing with? 
        elif '=' in l:
            # remove inline comments
            if '!' in l:
                l = l.split('!')[0]
            # Assign variable name and value to different objects
            split_line = l.split('=')
            curr_var = split_line[0].strip()
            curr_val = split_line[1].strip()
            # differentiate between single values and value arrays
            if curr_val[-1] == ',':
                curr_val = curr_val[:-1]
                param_dict[curr_block][curr_var] = [curr_val]
            else:
                param_dict[curr_block][curr_var] = curr_val
        elif l == '/':
            # block end signifiers
            pass
        else:
            # remove inline comments
            if '!' in l:
                l = l.split('!')[0]
            # most value arrays are split across multiple lines
            curr_val = l
            if curr_val[-1] == ',':
                curr_val = curr_val[:-1]
            param_dict[curr_block][curr_var].append(curr_val)
    
    # final pass through dict to detect and split up single-line value arrays  
    for block in param_dict.keys():
        for sub_block in param_dict[block].keys():
            vals = param_dict[block][sub_block]
            if type(vals) != list:
                if ',' in vals:
                    new_vals = vals.split(',')
                    param_dict[block][sub_block] = new_vals

    return param_dict

def read_cell_params(var_type, file_name=None):

    if not file_name:
        if var_type == 'phyto':
            test_fn = 'sample_aed2_phyto_pars.nml'
        elif var_type == 'pathogen':
            test_fn = 'sample_aed2_pathogen_pars.nml'
        else:
            sys.exit("Illegal var type specified")

        test_path = os.path.join(os.path.dirname(os.getcwd()), 'test_files')
        file_name = os.path.join(test_path, test_fn)
    
            
    with open(file_name, 'r') as f:
        content = f.read().split("\n")
    content = filter(None, content)
    content = filter(exclamationNotFirst, content)
    one_block_name = content[0][1:]
    content = content[1:]
    one_sub_block = 'pd'
    param_dict = {one_block_name: {one_sub_block:OrderedDict()}}
    for l in content:
        if l == '/' or l[0] == '&':
            pass
        elif '=' in l:
            split_line = l.split('=')
            one_sub_block = split_line[0].strip()
            param_line = split_line[1].split(",")
            param_list = [i.strip() for i in param_line if i != '']
            print len(param_list)
            param_dict[one_block_name][one_sub_block][param_list[0]] = param_list
        else:
            param_line = l.split(",")
            param_list = [i.strip() for i in param_line if i != '']       
            param_dict[one_block_name][one_sub_block][param_list[0]] = param_list
            print len(param_list)
    return param_dict
    
def read_zoo_params(file_name=None):
    if not file_name:
        test_fn = 'sample_aed2_zoop_pars.nml'
        test_path = os.path.join(os.path.dirname(os.getcwd()), 'test_files')
        file_name = os.path.join(test_path, test_fn)
    
    with open(file_name, 'r') as f:
        content = f.read().split("\n")
    content = filter(None, content)
    content = filter(exclamationNotFirst, content)
    one_block_name = content[0][1:]
    content = content[1:]
    param_dict = {one_block_name: {'pd':OrderedDict()}}
    for l in content:
        if l == '/' or l[0] == '&':
            pass
        else:
            param_line = l.split(",")
            param_list = [i.strip() for i in param_line if i != '']       
            param_dict[one_block_name]['pd'][param_list[0]] = param_list
            print len(param_list)
    return param_dict
    
def write_aed2_eco_configs(config_dict, out_path, config_type):
    """Write zoo / pathogen / phyto configs"""
    block_name = '&'+config_dict.keys()[0]+"\n"
    first_level = config_dict.keys()[0]
    second_level = config_dict[first_level].keys()[0]
    real_keys = config_dict[first_level][second_level].keys()
    
    if config_type != 'zoo':
        first_line_prefix = 'pd = '
    else:
        first_line_prefix = ''
    
    with open(out_path, 'w') as out_handle:
        out_handle.write(block_name)
        out_handle.write(first_line_prefix)
        for idx, rk in enumerate(real_keys):
            real_value = config_dict[first_level][second_level][rk]
            real_string = ",\t".join(real_value) + "\n"
            if idx != 0 and config_type != 'zoo':
                out_handle.write("\t"+real_string)
            else:
                out_handle.write(real_string)
        out_handle.write("/")
    return 0

def write_aed2_config(param_dict, out_path):
    block_prints = ['&'+i for i in param_dict.keys()]
    block_names = param_dict.keys()
    all_lines = []
    for bp, bn in zip(block_prints, block_names):
        all_lines.append(bp)
        sub_blocks = param_dict[bn].keys()
        for sb in sub_blocks:
            this_line = "\t"+sb+' = '
            this_val = param_dict[bn][sb]
            if type(this_val) == list:
                val_string = ", ".join(this_val)
            elif type(this_val) != str and type(this_val) != list:
                val_string = str(this_val)
            else:
                val_string = this_val
            this_line += val_string
            all_lines.append(this_line)
        all_lines.append("/")
    
    print_lines = [i+"\n" for i in all_lines]
    with open(out_path, 'w') as f:
        f.writelines(print_lines)
    return 0


param_dict = read_aed_nml()
phyto_dict = read_cell_params('phyto')
patho_dict = read_cell_params('pathogen')
zoo_dict = read_zoo_params()

test_file_dir = os.path.join(os.path.dirname(os.getcwd()), 'test_files')
phyto_out_file = os.path.join(test_file_dir, 'auto_written_phyto_pars.nml')
_ = write_aed2_eco_configs(phyto_dict, phyto_out_file, 'phyto')

patho_out_file = os.path.join(test_file_dir, 'auto_written_pathogen_pars.nml')
_ = write_aed2_eco_configs(patho_dict, patho_out_file, 'pathogen')

zoo_out_file = os.path.join(test_file_dir, 'auto_written_zoo_pars.nml')
_ = write_aed2_eco_configs(zoo_dict, zoo_out_file, 'zoo')

aed_out_file = os.path.join(test_file_dir, 'auto_aed2.nml')
_ = write_aed2_config(param_dict, aed_out_file)

