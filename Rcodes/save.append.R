## some save functions

save.append <- function(list=list, file="test.RData", delete=FALSE,
  description="first load"){
  
  mydate <- date()
  if(delete == TRUE){
    if(file.exists(file)){
      file.remove(file)
    }
  }

  if(file.exists(file)){
    ob <- load(file);
    if(exists("log")){
      log[[mydate]] <- data.frame(object=list, des=description)
      save(list=c(list, ob), file=file)
    }else{
      log <- list()
      log[[mydate]] <- data.frame(object=list, des=description)
      save(list=c(list, ob, "log"), file=file)
    }
    message("list appended to objects!")
    
  }else{
    log <- list()
    log[[mydate]] <- data.frame(object=list, des=description)
    
    message("new RData created with log!")
    save(list=c(list, "log"), file=file)  
  }
}

