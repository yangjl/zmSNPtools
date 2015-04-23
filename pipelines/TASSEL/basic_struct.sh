#!/bin/bash
#SBATCH -D /home/cjfiscus/projects/structure/ # running directory
#SBATCH -o /home/cjfiscus/projects/structure/slurm-log/10krandSNPs-stdout-%j.txt
#SBATCH -e /home/cjfiscus/projects/structure/slurm-log/10krandSNPs-stderr-%j.txt
#SBATCH -J randSNPs
set -e
set -u

module load structure-console

structure 

# runs with options configured in mainparams and extraparams in running directory
# ^^ these files need to be present in running directory for this script to work 
#
# additional options can be added after "structure" above to change parameters in mainparams and extraparams without editing files