#!/bin/bash

# set array to max K
#SBATCH --array=1-20

# allocate more memory (solving memory error)
#SBATCH --ntasks=2

#Set job options 
#SBATCH -D /home/cjfiscus/projects/structure/array_run # job dir
#SBATCH -o /home/cjfiscus/projects/structure/array_run/slurm-log/struct10k-stdout-%A_%a.txt # stdout 
#SBATCH -e /home/cjfiscus/projects/structure/array_run/slurm-log/struct10k-stderr-%A_%a.txt # stderr
#SBATCH -J strut10
set -e
set -u

# load the module to use
module load structure-console

# run the module using K = array task id, set separate output file for each run to prevent overwriting
structure -K $SLURM_ARRAY_TASK_ID -o outfile_$SLURM_ARRAY_TASK_ID_f
# above options overwrite options set in mainparams 