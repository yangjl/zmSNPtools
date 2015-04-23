#!/bin/bash

#SBATCH -D /home/cjfiscus/projects/tassel/ # running directory
#SBATCH -o /home/cjfiscus/projects/tassel/slurm-log/glm-stdout-%j.txt
#SBATCH -e /home/cjfiscus/projects/tassel/slurm-log/glm-stderr-%j.txt
#SBATCH --ntasks=2 # allocate 16 GB ram
#SBATCH -J glm
#SBATCH--mail-user=cjfiscus@ucdavis.edu
#SBATCH--mail-type=END # emails when done
#SBATCH--mail-type=FAIL # emails if fails
set -e
set -u

# This script uses TASSEL to run a GWAS using GLM
# Hapmap file must be in format FILENAME.hmp.txt for this to work
# SNPs in Hapmap file must be in order
# Must include phenotype (.txt) and Q matrix (.txt) 
# Runs on bigmem

module load jdk
module load tassel/5

run_pipeline.pl -Xmx16g -fork1 -h ./SNP_all_lines.hmp.txt -fork2 -importGuess ./phenotype.txt -fork3 -importGuess ./qmatrix_3.txt -combine4 -input1 -input2 -input3 -intersect -glm -glmOutputFile ./glm_3 -runfork1 -runfork2 -runfork3


# Runs TASSEL pipeline with a minimum of 512 MB RAM and a maximum of 16 GB RAM (2 compute node on bigmem)
# imports file SNP… as hapmap (-h)
# outputs results to the same directory under the name “glm_3”