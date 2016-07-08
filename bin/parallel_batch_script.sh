#!/bin/bash -l

#SBATCH

#SBATCH --job-name=no24_4_
#SBATCH --time=00:30:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=parallel_glm24qNodes1.err
#SBATCH --output=parallel_glm24qNodes1.out
#SBATCH --nodes=1
#SBATCH --partition=parallel

./run01.sh 
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh

./run01.sh
./run02.sh
