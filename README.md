zmSNPtools
==========

A collection of python, perl and R codes designed to process a range of basic, \
large-scale analyses for SNPs of maize diversity panel.

Usage:
=========
Copy the piece of python code to the directory in your $PATH:
for example:
```
$: cp impute4diallel/impute4diallel.py ~/bin/impute4diallel

```
Then, you will be ready to run it
```
# to find the help:
$: impute4diallel -help
```

Description
=========
## fpSNP
find a set of population specific fingerprint SNPs

#### update of the program
10/3/2014: v0.1 python code tested with a small set of simulated data

#### usage:

```
git clone git@github.com:yangjl/zmSNPtools.git
cd to/your/path
cp zmSNPtools/packages/fpSNP/fpSNP_v0.1.py zmSNPtools/bin/fpSNP
chmod +x zmSNPtools/bin/fpSNP
export PATH=$PATH:full/path/of/zmSNPtools/bin/
fpSNP -h

```

## gerpIBD
Compute the GERP score in a given IBM regions across different genotypes.

#### update of the program
12/17/2014: beta version for testing run

#### impute4diallel
SNP Imputation for maize diallel population from parental SNP panel to F1 hybrids.

## dsf2GWAS
transform density SNP format (DSF) to GWAS (gensel4.2 supported) and other formats

## top2ref.pl
translate illumina top/bot SNP call to ref/alt based on the reference genome.
#### usage:
```
##cloned from: https://github.com/VertebrateResequencing/vr-codebase
##The overall codebase developed and used by the Vertebrate Resequencing group at the Sanger Institute
# added perl module
export PERL5LIB=$PERL5LIB:~/Documents/Github/zmSNPtools/modules
cp packages/top2ref/top2ref.pl bin/top2ref
# use samtools to index reference genome in fasta format
sed -i -- 's/chromosome:AGPv2:.*chromosome /chr/g' ZmB73_RefGen_v2.fasta
samtools faidx ZmB73_RefGen_v2.fasta

cat snps.top | topbot-to-fwd-strand -r ref.fa > snps.map 
```







Update
=========
### General
1. changed gitignore
2. unwatch the MyPack and test folders


#### packages:
1. updated a SNP merge package **snp3merge**: merge hapmap1, hapmap2 and RNA-seq data
2. updated **impute4diallel**
3. updated **snpfrq**
4. updated **merge4ames**

#### MyPack
1. updated and tested **snp_sampling** program for XBSA project;


