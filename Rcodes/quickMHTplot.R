# Jinliang Yang
# Purpose: quick plot of GWAS results
# start: 2.11.2012
# updated: 7/14/2014
# add the RNA-seq LM regression data

# location: 129.186.85.7
quickMHTplot <- function(res=res, cex=.9, pch=16, col=rep(c("slateblue", "cyan4"), 5), 
                         GAP=5e+06, yaxis=NULL,
                         col2plot="ModelFreq", ... ){
  
  source("~/Documents/Rcodes/newpos.R")
  res <- newpos(res, GAP = GAP)
  chrtick <- chrline_tick(GAP = GAP)
  
  #### setup the cavon
  if(is.null(yaxis)){
    plot(x=-1000, y=-1000,  type="p", xaxt="n", xlab="", 
         xlim=c(0, max(chrtick$chrlines)), ylim=c(0, max(res[, col2plot], na.rm=TRUE)*1.3 ),
         ...)
  }else{
    plot(x=-1000, y=-1000,  type="p", xaxt="n", yaxt="n", xlab="",
         xlim=c(0, max(chrtick$chrlines)),
         ...)
    axis(side=2, at=yaxis, labels=yaxis)
  }
  axis(side=1, at=chrtick$ticks, labels=c("chr1", "chr2", "chr3", "chr4", "chr5", 
                                          "chr6", "chr7", "chr8", "chr9", "chr10"))
  abline(v=chrtick$chrlines, col="grey")
  
  for(i in 1:10){
    points(x=subset(res, chr==i)$newpos, y=res[res$chr==i, col2plot],
         pch = pch, col=col[i], cex=cex);
  }  
}
