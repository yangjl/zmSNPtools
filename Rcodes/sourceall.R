sourceall <- function(){
    file <- list.files(pattern="[.]R$", path="lib", full.names=TRUE)
    for(i in 1:length(file)){
        try(source(file[i]))
    }
}
