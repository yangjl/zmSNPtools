# This code was stolen from Turner
# last modified: 8/15, 2011

# To source this file into an environment to avoid cluttering the global workspace, put this in Rprofile:
# my.env <- new.env(); sys.source("/Users/yangjl/Documents/workingSpace/MyRCodes/Rprofile_general.r", my.env); attach(my.env)

#-----------------------------------------------------------------------
#     Load packages, set options and cwd, set up database connection
#-----------------------------------------------------------------------

## Load packages
# require(ggplot2)  #plotting
# require(plyr)		#data manipulation
# require(reshape)	#data manipulation
# require(sqldf)	#manipulate data frame with SQL
# require(Hmisc)	#frank harrell's miscellaneous functions

## Sets the working directory to C:/R
setwd("/Users/yangjl/Documents/workingSpace/")

## Don't show those silly significanct stars
#options(show.signif.stars=FALSE)

## Do you want to automatically convert strings to factor variables in a data.frame?
#options(stringsAsFactors=TRUE)

## Hard code the US repository for CRAN so it doesn't ask me every time.
#r <- getOption("repos")             
#r["CRAN"] <- "http://cran.us.r-project.org"
#options(repos = r)
#rm(r)

## Some SQLite stuff I don't use any more because I switched to MySQL
# require(RSQLite)
# channel <- dbConnect(SQLite(), "C:/cygwin/home/sturner/dbs/sdt.sqlite")
# query <- function(...) dbGetQuery(channel,...)

## Set up ODBC connection for MySQL localhost, and make it easy to query a database with query() function.
#require(RODBC) # The rprofile script will fail here if you don't have RODBC installed.
#channel <- odbcConnect("localhost")
#query <- function(...) sqlQuery(channel, ...)


#-----------------------------------------------------------------------
#                             Functions
#-----------------------------------------------------------------------

## Convert selected columns of a data frame to factor variables
factorcols <- function(d, ...) lapply(d, function(x) factor(x, ...)) 

## Returns a logical vector TRUE for elements of X not in Y
"%nin%" <- function(x, y) !(x %in% y) 

## Returns names(df) in single column, numbered matrix format.
n <- function(df) matrix(names(df)) 

## Single character shortcuts for summary() and head().
s <- base::summary
h <- utils::head

## ht==headtail, i.e., show the first and last 10 items of an object
ht <- function(d) rbind(head(d,10),tail(d,10))

## Show the first 5 rows and first 5 columns of a data frame or matrix
hh <- function(d) d[1:5,1:5]

## Read data on clipboard.
read.cb <- function(...) read.table(file="clipboard", ...) 

## Paste without a separator.
glue <- function(...) paste(...,sep="")

## name("test.png") results in "C:/R/2010-04-20-test.png" if running this in C:/R on April 20 2010.
name <- function(filename="filename") glue(getwd(),"/",Sys.Date(),"-",filename)

## Takes a dataframe and a column name, and moves that column to the front of the DF.
moveColFront <- function(d=dataframe, colname="colname") {
	index <- match(colname, names(d))
	cbind(d[index],d[-index])
}

## Permutes a column in a data.frame, sets seed optionally
permute <- function (dataframe, columnToPermute="column", seed=NULL) {
	if (!is.null(seed)) set.seed(seed)
	colindex <- which(names(dataframe)==columnToPermute)
	permutedcol <- dataframe[ ,colindex][sample(1:nrow(dataframe))]
	dataframe[colindex] <- permutedcol
	return(dataframe)
}

## Summarize missing data in a data frame. Return a list (lpropmiss) or data frame (propmiss)
lpropmiss <- function(dataframe) lapply(dataframe,function(x) data.frame(nmiss=sum(is.na(x)), n=length(x), propmiss=sum(is.na(x))/length(x)))
propmiss <- function(dataframe) {
	m <- sapply(dataframe, function(x) {
		data.frame(
			nmiss=sum(is.na(x)), 
			n=length(x), 
			propmiss=sum(is.na(x))/length(x)
		)
	})
	d <- data.frame(t(m))
	d <- sapply(d, unlist)
	d <- as.data.frame(d)
	d$variable <- row.names(d)
	row.names(d) <- NULL
	d <- cbind(d[ncol(d)],d[-ncol(d)])
	return(d[order(d$propmiss), ])
}

## Make a pretty QQ plot of p-values
qq = function(pvector, ...) {
    if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
	pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
    o = -log10(sort(pvector,decreasing=F))
	#e = -log10( 1:length(o)/length(o) )
	e = -log10( ppoints(length(pvector) ))
	plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
	abline(0,1,col="red")
}

## Get the proportion variation explained. See this website for more details: http://goo.gl/jte8X
rsq <- function(predicted, actual) 1-sum((actual-predicted)^2)/sum((actual-mean(actual))^2)

## Correlation matrix with p-values. See http://goo.gl/nahmV for documentation of this function
cor.prob <- function(X, dfr = nrow(X) - 2) {
     R <- cor(X)
	 above <- row(R) < col(R)
	 r2 <- R[above]^2
	 Fstat <- r2 * dfr / (1 - r2)
	 R[above] <- 1 - pf(Fstat, 1, dfr)
	 R[row(R)==col(R)]<-NA
     R
}

## This function accepts a GLM object and does a LR chi-square test on the fit.
lrt <- function (modelobject) {
	lrtest.chi2 <- model$null.deviance - model$deviance # Difference in deviance between model with intercept only and full model.  This is the likelihood ratio test statistic (-2(log(L))).
	lrtest.df   <- model$df.null - model$df.residual # Difference in DF.  Make sure this equals the number of predictors in the model!
	fitpval     <- 1-pchisq(lrtest.chi2,lrtest.df)
	cat("Likelihood ratio test on model fit:\n\n")
	data.frame(lrtest.chi2=lrtest.chi2,lrtest.df=lrtest.df,fitpval=fitpval) #Output gives you the chisquare, df, and p-value.
}

## This function does the same thing as lrtest in the Design package, but doesn't do as much checking.
## Remember, the lrt has to test the same model (model fit on same observations)
## Also the drop1(fullmodel,test="Chisq") does something similar.
lrt2 <- function (full,reduced) {
	if (reduced$deviance<=full$deviance) stop ("Reduced model not worse than full.")
	if (reduced$df.residual<=full$df.residual) stop ("Reduced model doesn't have more degrees of freedom.")
	lrtest.chi2 <- reduced$deviance-full$deviance
	lrtest.df   <- reduced$df.residual - full$df.residual
	fitpval        <- 1-pchisq(lrtest.chi2,lrtest.df)
	cat("Likelihood ratio test on two models:\n\n")
	data.frame(lrtest.chi2=lrtest.chi2,lrtest.df=lrtest.df,fitpval=fitpval)
}	

## This gets the overall anova p-value out of a linear model object
lmp <- function (modelobject) {
	if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
	f <- summary(modelobject)$fstatistic
	p <- pf(f[1],f[2],f[3],lower.tail=F)
	attributes(p) <- NULL
	return(p)
}

## Function for arranging ggplots. use png(); arrange(p1, p2, ncol=1); dev.off() to save.
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
	dots <- list(...)
	n <- length(dots)
	if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
	if(is.null(nrow)) { nrow = ceiling(n/ncol)}
	if(is.null(ncol)) { ncol = ceiling(n/nrow)}
        ## NOTE see n2mfrow in grDevices for possible alternative
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
	ii.p <- 1
	for(ii.row in seq(1, nrow)){
	ii.table.row <- ii.row	
	if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
		for(ii.col in seq(1, ncol)){
			ii.table <- ii.p
			if(ii.p > n) break
			print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
			ii.p <- ii.p + 1
		}
	}
}

## Imputes the median value of a vector, matrix, or data frame.
## Stolen from na.roughfix function in the randomForest package.
medianImpute <- function (object=NULL, ...) {
	if (class(object) == "data.frame") {
		isfac <- sapply(object, is.factor)
		isnum <- sapply(object, is.numeric)
		if (any(!(isfac | isnum))) 
			stop("dfMedianImpute only works for numeric or factor")
		roughfix <- function(x) {
			if (any(is.na(x))) {
				if (is.factor(x)) {
					freq <- table(x)
					x[is.na(x)] <- names(freq)[which.max(freq)]
				}
				else {
					x[is.na(x)] <- median(x, na.rm = TRUE)
				}
			}
			x
		}
		object[] <- lapply(object, roughfix)
		return(object)
	}
	else if(is.atomic(object)) {
		d <- dim(object)
		if (length(d) > 2) 
			stop("vectorMedianImpute can't handle objects with more than two dimensions")
		if (all(!is.na(object))) 
			return(object)
		if (!is.numeric(object)) 
			stop("vectorMedianImpute can only deal with numeric data.")
		if (length(d) == 2) {
			hasNA <- which(apply(object, 2, function(x) any(is.na(x))))
			for (j in hasNA) object[is.na(object[, j]), j] <- median(object[, 
				j], na.rm = TRUE)
		}
		else {
			object[is.na(object)] <- median(object, na.rm = TRUE)
		}
		return(object)
	}
	else stop("Object is not a data frame or atomic vector")
}

## Cast a string of letters to an array of integers/10. 
## Useful for converting genotypes in string format to array of imputed dosages.
## E.g. If homozygous major allele coded 0, heterozygote 1, homozygote minor 2, missing NA
##      then you can save a ton of space by recoding genotypes as a text string with a single
##      character (A through U) representing an imputed genotype dosage between 0 and 2. 
##      Missing values can be coded zero in the text string, and are recoded to NA here.
## Here's an example:
## > cast("ABCD0K0U")
## [1] 0.0 0.1 0.2 0.3  NA 1.0  NA 2.0
cast <- function(letters) {
    as.numeric(sapply(strsplit(letters,NULL), function(x) {
        x<-sub("0",NA,x)
        x<-sub("A",0.0,x)
        x<-sub("B",0.1,x)
        x<-sub("C",0.2,x)
        x<-sub("D",0.3,x)
        x<-sub("E",0.4,x)
        x<-sub("F",0.5,x)
        x<-sub("G",0.6,x)
        x<-sub("H",0.7,x)
        x<-sub("I",0.8,x)
        x<-sub("J",0.9,x)
        x<-sub("K",1.0,x)
        x<-sub("L",1.1,x)
        x<-sub("M",1.2,x)
        x<-sub("N",1.3,x)
        x<-sub("O",1.4,x)
        x<-sub("P",1.5,x)
        x<-sub("Q",1.6,x)
        x<-sub("R",1.7,x)
        x<-sub("S",1.8,x)
        x<-sub("T",1.9,x)
        x<-sub("U",2.0,x)
        x
    }
    ))
}

# Did you make it this far?
message("\n******************************\nSuccessfully loaded Rprofile_general.r\n******************************")