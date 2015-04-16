snpid_chrpos <- function(df=col3, which_dir="snpid2pos"){
  if(which_dir == "snpid2pos"){
    if(sum(names(df) %in% "snpid") ==1){
      df$chr <- as.numeric(as.character(gsub("_.*", "", df$snpid)))
      df$pos <- as.numeric(as.character(gsub(".*_", "", df$snpid)))
    }else{
      stop("df should have col: 'snpid'! ")
    }
  }
  if(which_dir == "pos2snpid"){
    if(sum(names(df) %in% c("chr", "pos")) == 2){
      df$snpid <- paste(df$chr, df$pos, sep="_")
    }else{
      stop("df should have cols: 'chr' and 'pos'! ")
    }
  }
  return(df)
}
