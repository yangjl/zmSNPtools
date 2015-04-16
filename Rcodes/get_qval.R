#

get_qval <- function(pval=pval1, pcol="P1df", method="fdr"){
  qval <- p.adjust(pval[,pcol], method=method)
  return(qval)
}

