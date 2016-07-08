#!/bin/bash -l

#SBATCH
#SBATCH --job-name=glm_2_nodes_10_calls_parallel.sh
#SBATCH --time=03:00:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --nodes=2
#SBATCH --partition=parallel

