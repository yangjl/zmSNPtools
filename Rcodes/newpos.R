
newpos <- function (dataframe, GAP = 5e+06, version = "v3") 
{
  d <- dataframe
  if (!("chr" %in% names(d) & "pos" %in% names(d))){
    stop("Make sure your data frame contains columns chr and pos")
  }
    
  if(version == "v2"){
      cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v2.csv")
  }
  if(version == "v3"){
      cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v3.csv")
  }
  
  cl$accumpos <- cl$BP
  cl <- cl[order(cl$CHR), ]
  d$newpos <- d$pos;
  for (i in 2:10) {
    cl[cl$CHR == i, ]$accumpos <- cl[cl$CHR == (i - 1), ]$accumpos + cl[cl$CHR == i, ]$accumpos + GAP
    d[d$chr == i, ]$newpos <- d[d$chr == i, ]$pos + cl[cl$CHR == (i - 1), ]$accumpos + GAP
  }
  return(d)
}

chrline_tick <- function(GAP=5e+06, version = "v3"){
  #xscale:
  if(version == "v2"){
      cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v2.csv")
  }
  if(version == "v3"){
      cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v3.csv")
  }
  
  names(cl) <- c("chr", "snp", "pos")
  cl <- newpos(cl, GAP=GAP)
    
  cl$ticks <- cl$pos[1]/2
  cl$chrlines <- cl$pos[1]+GAP/2
  for(i in 2:10){
    cl$ticks[i] <- cl$newpos[i-1] + (cl$newpos[i]-cl$newpos[i-1])/2;
    cl$chrlines[i] <- cl$newpos[i]+ GAP/2;  
  }
  return(cl)
}

#cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v3.csv")
#cl$SNP <- "ATCG"
#cl <- cl[, c(1,3,2)]
#write.table(cl, "~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v3.csv", sep=",", row.names=FALSE, quote=FALSE)
