#!/bin/bash -l

#SBATCH

#SBATCH --job-name=kaw_build
#SBATCH --time=24:00:00
#SBATCH --nodes=24
#SBATCH --mail-type=end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=build_db.err
#SBATCH --output=build_db.out

module load python/2.7.10

MPATH="/home-3/karoraw1@jhu.edu/work/kaw/GLM_Wrapper/Giovanni"

python $MPATH/GiovanniMetaData.py 2>&1
