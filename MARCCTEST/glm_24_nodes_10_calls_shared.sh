#!/bin/bash -l

#SBATCH
#SBATCH --job-name=glm_24_nodes_10_calls_shared.sh
#SBATCH --time=03:00:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --nodes=24
#SBATCH --partition=shared

