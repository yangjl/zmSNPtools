#!/bin/bash

# sort a hapmap file by chromosome and position 

head -1 SNP_rand_10k.hapmap.txt > newhapmap.txt; for i in {1..10}; do grep "chr$i\t"  SNP_rand_10k.hapmap.txt | sort -n -k 4 >> newhapmap.txt; done