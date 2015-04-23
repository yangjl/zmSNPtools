#!/bin/bash

#SBATCH -D /home/cjfiscus/projects/tassel/ # running directory
#SBATCH -o /home/cjfiscus/projects/tassel/slurm-log/mlm-stdout-%j.txt
#SBATCH -e /home/cjfiscus/projects/tassel/slurm-log/mlm-stderr-%j.txt
#SBATCH -J mlmsansq
#SBATCH --ntasks=2 # allocate 16 GB ram/ 2 compute nodes
#SBATCH--mail-user=cjfiscus@ucdavis.edu
#SBATCH--mail-type=END # emails when done
#SBATCH--mail-type=FAIL # emails if fails
set -e
set -u

# This script uses TASSEL to run an MLM excluding Q
# Hapmap file must be in format FILENAME.hmp.txt for this to work
# SNPs in Hapmap file must be in order
# Runs on bigmem

module load jdk
module load tassel/5

run_pipeline.pl -Xmx16g -fork1 -k kinship.txt -fork2 -importGuess phenotype.txt -fork3 -h SNP_all_lines.hmp.txt -combine4 -input3 -input2 -intersect -combine5 -input4 -input1 -mlm -export mlm_output_noq -runfork1 -runfork2 -runfork3 -runfork4

# Runs TASSEL pipeline with a minimum of 512 MB RAM and a maximum of 16 GB RAM (2 compute nodes on bigmem)
# imports file SNP… as hapmap (-h)
# outputs results to the same directory under the name “mlm_output_noq*”