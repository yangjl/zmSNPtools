# load marker data
NEGSc2.GBS80=read.delim("NEGSc2.GBS80.MAF5.txt", header=T, na.string="NA", sep="\t")
dim(NEGSc2.GBS80)
NEGSc2.GBS80[1:5,1:10]
NEGSc2.GBS80=NEGSc2.GBS80[,3:382]
genotype=as.matrix(t(NEGSc2.GBS80))
dim(genotype)

#imputation of GBS genotype data with randomForest method
## load package "missForest"
library(itertools, lib="~/bin/Rlib")
library(randomForest, lib="~/bin/Rlib")
library(missForest, lib="~/bin/Rlib")
Markers_impute=missForest(genotype, maxiter = 10, ntree = 100, variablewise = FALSE, decreasing = FALSE, verbose = FALSE, replace = TRUE, classwt = NULL, cutoff = NULL, strata = NULL, sampsize = NULL, nodesize = NULL, maxnodes = NULL)
Markers_impute=Markers_impute$ximp
write.table(Markers_impute, file="NEGSc2_RFtr100.F01miss80.txt", sep=",")
