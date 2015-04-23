scripts
=======
10k_structarray_10.sh

Slurm script to use structure to infer population structure from SNP data using pre-installed module on cluster.  Runs as an array job for interval 1- max K.  Dependent on mainparams and extraparams files included with structure console.  K value in mainparams will be overwritted by $SLURM_ARRAY_TASK_ID. 
-----

basic_struct.sh

Basic slurm script to use structure to infer population structure from SNP data using pre-installed module on cluster.  Runs as a single job using parameters supplied by mainparams and extraparams files.  
-----

combine.sh

Experimental script used to attempt to develop a way of editing K value in mainparams file for structure input on cluster.  Requires that mainparams file is split before line to be processed.  Script is successful but does not work well on cluster due to the speed at which mainparams is accessed for each task in array.  Problem solved by 10k_structarray_10.sh (above). 
-----

con_to_fasta.sh

Slurm script to use samtools to convert a .bam file to a .fasta file in a different directory.  
-----

dup_snp_for_struct.py

Python 3 script written duplicate SNP data after it has been transposed.  This is to fulfill data format requirements for structure.  Use after file has been processed with hapmap_to_structuretxt
-----

export_selected_regions.py

Python 3 script to take a .txt file with SNPs (either a single chromosome or absolute position in genome) and export SNP data for SNPs in a designated genomic region. 
-----

export_snp_names_for_multiple_regions(_v2).py

Python 3 script to export SNP names from MULTIPLE given positions in the genome using specially formatted input file.  Similar to export_selected_regions.py. See script for details.

VERSION 2 adds ability to specify boundaries to add to given bounds (ie +-5000 bp)
-----

export_SNP_names.py

Functionally equivalent to export_selected_regions.py but only exports the SNP name (for highlighting in R). 
-----

find_sig_snps.py

Python 3 script to export SNP data from GWAS results based on a supplied P value threshold.
-----  

format_tassel_to_r_graph.py

Python 3 script that takes GWAS (GLM or MLM) output from TASSEL and formats it for graphing with R.  Output is file with columns SNP, CHR, BP, P. 
-----

hapmap_to_structuretxt

R script written by rossibarra to transpose SNP hapmap table.  Use before dup_snp_for_struct.py.
-----

pdf_to_jpeg.app

Automator application to convert pdf files to jpeg image files.  Used to process graphs created by R and saved as large pdfs. Must be run on Macintosh OS.   
-----

q_matrix.sh 

Bash script that formats a Q matrix from structure for use in TASSEL.  Prompts users for number of populations (K) and adds the appropriate heading to the file.  
-----

slurm_glm.sh

Slurm script to run GWAS using General Linear Model (GLM) in TASSEL.  
-----

slurm_mlm_sansq.sh

Slurm script to run GWAS using MLM in TASSEL but excluding the Q matrix.  Useful for comparisons between Q + K and K models.  
-----

slurm_mlm.sh

Slurm script to run GWAS using Mixed Linear Model (MLM) in TASSEL. 
-----

sort_hapmap_chr_pos.sh

Bash script written by rossibarra to sort a hapmap file with randomly selected scripts by both chromosome and position.  This sorting is required to use the hapmap file in TASSEL.    
-----

split_hapmap_by_chr.py

Python 3 script written to extract SNP data by chromosome to separate files for sorting.  
-----

tassel_load_results.sh

Bash script to import table results into TASSEL.  Useful for graphing results of previous analyzes.   
-----

tassel.sh

Bash script to run a GWAS using GLM in TASSEL pipeline using the tutorial data.    
-----