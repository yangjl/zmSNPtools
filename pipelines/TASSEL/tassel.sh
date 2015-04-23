#!/bin/bash

cd /Users/cjfiscus/tassel4-standalone # change working directory

# run GLM using TASSEL tutorial data
#
# importGuess is used but when analyzing data using TASSEL the filetype should be imported using the specific command.  
# For instance, without using -h for hapmap file TASSEL will not import all of the data, just something like first 30 SNPs

./run_pipeline.pl -fork1 -importGuess ./Tut_Data/d8_sequence.phy -fork2 -importGuess ./Tut_Data/mdp_population_structure.txt -fork3 -importGuess ./Tut_Data/mdp_traits.txt -combine4 -input1 -input2 -input3 -intersect -glm -glmMaxP 1e-3 -glmOutputFile ./test -runfork1 -runfork2 -runfork3

# run_pipeline.pl >> running the actual command line script for TASSEL
# forking is important so that you can do multiple things at a time 
# d8_sequence is hapmap seq., mdp_population_structure is q matrix, mdp_traits is phenotype file  
