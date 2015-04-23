#!/bin/bash
#SBATCH -D /home/cjfiscus/projects/gs/ # running directory
#SBATCH -o /home/cjfiscus/projects/gs/convert-stdout-%j.txt # std out dir
#SBATCH -e /home/cjfiscus/projects/gs/convert-stderr-%j.txt #std error dir
#SBATCH -J fasta # jobname
set -e
set -u

# load module to be used
module load samtools

# Convert .bam file to .fasta file using samtools
# Then copy this file to project directory
samtools view path_to_file.bam | awk '{OFS="\t"; print ">"$1"\n"$10}' - > /home/cjfiscus/projects/gs/W64A.fasta