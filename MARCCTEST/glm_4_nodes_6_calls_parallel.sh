#!/bin/bash -l

#SBATCH
#SBATCH --job-name=glm_4_nodes_6_calls_parallel.sh
#SBATCH --time=03:00:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --nodes=4
#SBATCH --partition=parallel

