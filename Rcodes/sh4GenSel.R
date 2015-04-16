# Jinliang Yang
# Purpose: control GenSel running
# date: July.10.2014
# location: server.9


sh4GenSel <- function(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
                      sh="CL_test.sh", pi=0.995, findsale ="no",
                      geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr", 
                      pheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
                      map="/Users/yangjl/Documents/linkage.map",
                      chainLength=1000, burnin=100, varGenotypic=1.4, varResidual=2
                      ){
  #####
  setwd(pwd)
  inp <- gsub("sh", "inp", sh);
  log <- gsub("sh", "log", sh)
  
  GenSel_inp(inp=inp, geno=geno, pheno=pheno, chainLength=chainLength, map=map, pi=pi,findsale ="no",
             burnin=burnin, varGenotypic=varGenotypic, varResidual=varResidual)
  
  cat(paste("#Note:", Sys.time(), sep=" "),
      paste("GenSel4R", inp, ">", log),
      paste("python ~/bin/send_email.py -s ", "'", sh, "'", sep=""),
      file=sh, sep="\n", append=FALSE);
  
  message(paste("In this path: ", pwd, sep=""), "\n",
          paste("RUN: sh ", sh))
  
}


GenSel_inp <- function(inp="CL_test.inp", pi=0.995, findsale ="no",
                       geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr", 
                       pheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
                       map="/Users/yangjl/Documents/linkage.map",
                       chainLength=1000, burnin=100, varGenotypic=1.4, varResidual=2){
  
    cat(paste("// gensel input file written", Sys.time(), sep=" "), 
        
        "analysisType Bayes",
        "bayesType BayesC",
        paste("chainLength", chainLength, sep=" "),
        paste("burnin", burnin=burnin, sep=" "),
        paste("probFixed", pi, sep=" "),
        
        paste("varGenotypic",  varGenotypic, sep=" "),
        paste("varResidual",  varResidual, sep=" "),
        "nuRes 10",
        "degreesFreedomEffectVar 4",
        "outputFreq 100",
        "seed 1234",
        "mcmcSamples yes",
        "plotPosteriors no",
        paste("FindScale", findsale),
        "modelSequence no",
        "isCategorical no",
        "linkageMap AGPv2",
        "addMapInfoToMarkers yes",
        "windowBV yes",
        "",
        "// markerFileName",
        paste("markerFileName", geno, sep=" "), 
        "",
        "// phenotypeFileName",
        paste("phenotypeFileName", pheno, sep=" "),
        "",
        "// mapOrderFileName",
        paste("mapOrderFileName", map, sep=" "),
        
        file=inp, sep="\n"
    )	
}
