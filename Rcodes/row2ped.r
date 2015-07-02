# Purpose: match up the barcoded data with pedigree
# author: Jinliang Yang
# update: 11/18/2011



###################################################################
# pwd: your pathway of input and output files
# inputpheno: a phenotype file with column named "Row"
# inputped: part of the FB with column "Row", and all "" been replaced by characters
# output: your output file name
###################################################################

row2ped <- function(pwd= "/Users/yangjl/Documents/workingSpace/Diallel_NAMF1/", 
inputpheno="Off-PVP-2011.csv",
inputped="FB_2011_offpvp.csv"){
	
# setting your working directory:
	setwd(pwd);
# read in the phenotype file:
	pheno <- read.csv(inputpheno, header=TRUE);
# read in the pedigree file:
	ped <- read.csv(inputped, header=TRUE);
# merge the file together
	mfile <- merge(pheno, ped, by="Row", all.x=TRUE);
# report the results:
	print(paste("your input row # is ", nrow(pheno), sep=""));
	print(paste(sum(is.na(mfile$Pedigree)), " rows could not find matched pedigree", sep=""));
# return the results:
	return(mfile);
	
}





