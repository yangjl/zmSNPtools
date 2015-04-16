

##########
#dsnp <- read.csv("~/MyProjects/KRN_GWAS_v3/GWAS3_proj/data/S.table10_elite_geno.csv")
###
#col3 <- dsnp[, 1:3]
#names(col3) <- c("snpid", "chr", "pos")
#source("~/Documents/Rcodes/snpid_chrpos.R")
#col3 <- snpid_chrpos(df=col3, which_dir="snpid2pos")
#elitegeno <- merge(col3, dsnp, by.x="snpid", by.y="SNPID")

#dsnp2GenABEL(dsf3=elitegeno, geno_cols=4:ncol(elitegeno), output="test2.raw")

dsnp2GenABEL <- function(dsf3=elitegeno, geno_cols=4:211, output="test.raw"){
  #DSF3: snpid, chr, pos
  #genotype must be coded as 0 (missing), 1 (for AA), 2 (for AB) and 3 (for BB)
  #ref:http://www.genabel.org/GenABEL/convert.snp.text.html
  if(sum(names(dsf3) %in% c("snpid", "chr", "pos")) == 3){
    
    #1st line: The first line of this file contains IDs of all study subjects.
    cat(names(dsf3)[geno_cols],file=output, append=FALSE, sep="\t")
    cat("\n",file=output, append=TRUE, sep="")
    cat(dsf3$snpid,file=output, append=TRUE, sep="\t")
    cat("\n",file=output, append=TRUE, sep="")
    cat(dsf3$chr,file=output, append=TRUE, sep="\t")
    cat("\n",file=output, append=TRUE, sep="")
    cat(dsf3$pos,file=output, append=TRUE, sep="\t")
    cat("\n",file=output, append=TRUE, sep="")
    
    geno <- t(dsf3[, geno_cols])
    for(i in 1:ncol(geno)){
      cat(
        #The 5th line contains the data for SNP, which is listed first on the second line. 
        #The first column of this line specifies the genotype for the person, 
        #who is listed first on the line 1; the second column gives the genotype 
        #for the second person, so on. The genotypes are coded as 0 (missing), 1 (for AA), 
        #2 (for AB) and 3 (for BB). Here is a small example:
        
        geno[,i],
        file=output, append=TRUE, sep="\t"
      )
      cat("\n", file=output, append=TRUE, sep="")
    }  
  }else{
    stop("You should input DSF3 format!")
  }
  
}

