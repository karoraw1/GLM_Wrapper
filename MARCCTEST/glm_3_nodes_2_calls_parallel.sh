#!/bin/bash -l

#SBATCH
#SBATCH --job-name=glm_3_nodes_2_calls_parallel.sh
#SBATCH --time=03:00:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --nodes=3
#SBATCH --partition=parallel

