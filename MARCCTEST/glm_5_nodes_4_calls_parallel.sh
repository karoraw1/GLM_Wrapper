#!/bin/bash -l

#SBATCH
#SBATCH --job-name=glm_5_nodes_4_calls_parallel.sh
#SBATCH --time=03:00:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --nodes=5
#SBATCH --partition=parallel

module load netcdf/intel/4.3.3.1 glm
cd /home-3/karoraw1@jhu.edu/work/kaw/GLM_Wrapper/GLM_Executables/examples_2.2/coldlake/fabm
glm

module load netcdf/intel/4.3.3.1 glm
cd /home-3/karoraw1@jhu.edu/work/kaw/GLM_Wrapper/GLM_Executables/examples_2.2/coldlake/fabm
glm

module load netcdf/intel/4.3.3.1 glm
cd /home-3/karoraw1@jhu.edu/work/kaw/GLM_Wrapper/GLM_Executables/examples_2.2/coldlake/fabm
glm

module load netcdf/intel/4.3.3.1 glm
cd /home-3/karoraw1@jhu.edu/work/kaw/GLM_Wrapper/GLM_Executables/examples_2.2/coldlake/fabm
glm

