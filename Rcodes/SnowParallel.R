## Jinliang Yang
## July 25th, 2014


##################################################################

SnowParallel <- function(range=1:1000, cpus=8, myfun=mrk2RNA, ...){
  ### need further test
  library(snow)
  library(snowfall)
  start <- proc.time();
  
  sfInit(parallel=TRUE, cpus=cpus)
  sfExportAll();
  sfLibrary(glmnet);
  output <- sfLapply(range, myfun, ...);
  sfRemoveAll()
  sfStop()
  
  return(output)
  cat("Total Run Time:", proc.time() - start, "\n");		
  
}

