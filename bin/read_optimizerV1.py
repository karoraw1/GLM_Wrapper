# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 21:17:27 2016

@author: login
"""
import numpy as np

opt_f = '../test_files/optimize_these.txt'

def optimize_lake(Lake, path):
    meta_config = Lake.glm_config
    block_key_option = []
    with open(opt_f, 'r') as f:
        for line in f:
            if line[0] == '#':
                pass
            elif line[0] == '&':
                new_block = line[1:].strip()
                assert new_block in meta_config.keys()
                print "\n", new_block
            elif line[0] == '-':
                splitted = line[2:].split('=')
                spit_cleaned = [i.strip() for i in splitted]
                sub_block = spit_cleaned[0]
                print "\t%r" % sub_block
                [val_space, dType] = spit_cleaned[1].split("],[")
                val_space, dType = val_space[1:], dType[:-1]
                val_split = [l.strip() for l in val_space.split(",")]
                
                if val_split[-1] == 'None':
                    val_split.pop()
                    
                if dType == 'float':
                    val_split = [float(l) for l in val_split]
                elif dType == 'int':
                    val_split = [int(l) for l in val_split]
                        
                if dType == 'bool':
                    for i in val_split:
                        print "\t\t%r" % i
                        block_key_option.append((new_block, sub_block, i))
                elif len(val_split) != 3 and dType != 'bool':
                    for i in val_split:
                        block_key_option.append((new_block, sub_block, i))
                        print "\t\t%r" % i
                elif len(val_split) == 3:
                    arange = np.arange(val_split[0], val_split[1]+val_split[2],
                                       val_split[2])
                    for i in arange:
                        j = int(np.floor(i*1000))
                        j = j/1000.
                        print "\t\t%r" % j
                        block_key_option.append((new_block, sub_block, j))
                else:
                    print "\n\nunparsed option\n\n"
                
                        
            else:
                print "Unparsable line"
                print line
    print "%i alternate cases observed" % len(block_key_option)
    
    return block_key_option
            