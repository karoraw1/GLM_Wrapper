#!/bin/bash -l

#SBATCH
#SBATCH --job-name=glm_16_nodes_20_calls_shared.sh
#SBATCH --time=03:00:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --nodes=16
#SBATCH --partition=shared

