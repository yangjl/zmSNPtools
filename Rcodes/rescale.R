#http://cran.r-project.org/doc/contrib/Lemon-kickstart/rescale.R
# linearly transforms a numeric object to fit a different range.
# could use a few more sanity checks

rescale<-function(x,newrange) {
  if(nargs() > 1 && is.numeric(x) && is.numeric(newrange)) {
    # if newrange has max first, reverse it
    if(newrange[1] > newrange[2]) {
      newmin<-newrange[2]
      newrange[2]<-newrange[1]
      newrange[1]<-newmin
    }
    xrange<-range(x)
    if(xrange[1] == xrange[2]) stop("can't rescale a constant vector!")
    mfac<-(newrange[2]-newrange[1])/(xrange[2]-xrange[1])
    invisible(newrange[1]+(x-xrange[1])*mfac)
  }
  else {
    cat("Usage: rescale(x,newrange)\n")
    cat("\twhere x is a numeric object and newrange is the min and max of the new range\n")
  }
}