## Jinliang Yang
## 10.8th, 2013

#Indel coding: A="-"deletion and T="+" insertion: 
chr10 <- read.delim("hmp12rna_chr10", header=TRUE, colClasses = classes)
dim(chr10)
#[1] 2220183      33
table(chr10$source)

hmp1 <- read.delim("../Data/hmp1_dsnp/maizeHapMapV1_B73RefGenV2_20110309_chr10.hmp.txt_filtered",header=TRUE)
hmp2 <- read.delim("../Data/hmp2_dsnp/maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt_filtered", header=TRUE)
rnaseq <- read.delim("../Data/rnaseq/rnaseq_snp+indel_chr10.txt", header=TRUE)
hmp1$source <- "hmp1"
hmp2$source <- "hmp2"
rnaseq$source <- "rna"

merge_test(randidx=2110, sour="hmp12")

########################################################
merge_test <- function(randidx=8497, sour ="hmp12", als=NULL){
  
  
  if(!is.null(als)){
    myid <- subset(chr10, source==sour & alleles %in% als)$id[randidx]
  } else{
    myid <- subset(chr10, source==sour)$id[randidx]
  }
  if(sour == 'hmp12'){
    tem <- rbind(subset(hmp1, id == myid),subset(hmp2, id == myid), subset(chr10, id == myid) )
  }
  
  if(sour == 'hmp1rna'){
    tem <- rbind(subset(hmp1, id == myid), subset(rnaseq, id == myid), subset(chr10, id == myid) )
  }
  
  if(sour == 'hmp2rna'){
    tem <- rbind(subset(hmp2, id == myid), subset(rnaseq, id == myid), subset(chr10, id == myid) )
  }
  if(sour == 'hmp12rna'){
    tem <- rbind(subset(hmp1, id == myid), subset(hmp2, id == myid),
                 subset(rnaseq, id == myid), subset(chr10, id == myid) )
  }
  return(tem)
}

