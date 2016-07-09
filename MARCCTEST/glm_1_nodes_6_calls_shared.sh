#!/bin/bash -l

#SBATCH
#SBATCH --job-name=glm_1_nodes_6_calls_shared.sh
#SBATCH --time=03:00:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --nodes=1
#SBATCH --partition=shared

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

module load netcdf/intel/4.3.3.1 glm
cd /home-3/karoraw1@jhu.edu/work/kaw/GLM_Wrapper/GLM_Executables/examples_2.2/coldlake/fabm
glm

module load netcdf/intel/4.3.3.1 glm
cd /home-3/karoraw1@jhu.edu/work/kaw/GLM_Wrapper/GLM_Executables/examples_2.2/coldlake/fabm
glm

