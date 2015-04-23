# Jinliang Yang
# Purpose: control GenSel running
# date: April 2015


collect_gsout <- function(dir = "slurm-scripts/wholeset", fileptn ="out"){
  files <- list.files(path = dir, pattern = fileptn)
  
  out <- data.frame(file=files, genvar=-9, resvar=-9, totvar=-9, h2=-9)
  for(i in 1:nrow(out)){
    text <- readLines(paste(dir, out$file[i], sep="/"))
    val1 <- grep("Residual Variance", text, value=TRUE)
    val1 <- gsub(".* = ", "", val1)
    
    val2 <- grep("Genetic  Variance", text, value=TRUE)
    val2 <- gsub(".* = ", "", val2)
    
    val3 <- grep("Total Variance", text, value=TRUE)
    val3 <- gsub(".* = ", "", val3)
    
    val4 <- grep("Proportion of Variance a/c Markers", text, value=TRUE)
    val4 <- gsub(".* = ", "", val4)
    
    out$resvar[i] <- val1
    out$genvar[i] <- val2
    out$totvar[i] <- val3
    out$h2[i] <- val4
  }
  return(out)
}
