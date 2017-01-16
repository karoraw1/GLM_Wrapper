# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 18:18:56 2016

@author: login
"""
import os

base_name = 'glm_'
base_path = '../MARCCTEST'
shared_node_range = [1, 2, 4, 8, 16, 24]
parallel_node_range = list(range(1,7))
execution_range = [2, 4, 6, 8, 10, 20, 40, 80]
partitions = ["shared", "parallel"]
trios = []
for j in execution_range:
    for i, l in zip(shared_node_range, parallel_node_range):
        for k in partitions:
            trios.append((i, j, k, l))
        
edits = ["file_handle", "job_name", "nodes", "partition"]

replacement_pack = {i:[] for i in edits}
command_steps = ["module load netcdf/intel/4.3.3.1 glm",
                 "cd /home-3/karoraw1@jhu.edu/work/kaw/GLM_Wrapper/GLM_Executables/examples_2.2/coldlake/fabm",
                 "glm"]
master_name = base_path+"/master_batch.sh"
if os.path.exists(master_name):
    os.remove(master_name)
    print "Removed"+master_name
                 
handle2 = open(base_path+"/master_batch.sh", "a")

for i in trios:
    if i[2] == 'shared':
        nodes = i[0]        
    else:
        nodes = i[3]        
    fn = base_name+str(nodes)+"_nodes_"+str(i[1])+"_calls_"+i[2]+".sh"
    if not os.path.exists(os.path.join(base_path, fn)):
        handle = open(os.path.join(base_path, fn), "w+")
        replacement_pack['file_handle'].append(handle)
        replacement_pack['nodes'].append(str(nodes))
        replacement_pack['partition'].append(str(i[2]))
        handle.write("#!/bin/bash -l\n\n#SBATCH\n#SBATCH --job-name=")
        handle.write(str(fn))
        handle.write("\n#SBATCH --time=03:00:00\n#SBATCH --mail-type=begin,end\n")
        handle.write("#SBATCH --mail-user=karoraw1@jhu.edu\n""#SBATCH --nodes=")
        handle.write(str(nodes))
        handle.write("\n#SBATCH --partition=")
        handle.write(str(i[2]))
        handle.write("\n\n")
        for execs in range(i[1]):
            for line in range(len(command_steps)):
                handle.write(command_steps[line]+"\n")
            handle.write("\n")
        
        handle.close()
    handle2.write("sbatch "+fn+"\n")
    
handle2.close()

