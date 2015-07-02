
#================= Functions used ========================================

# Function asking the user to selec an integer among given intergers
userQueryInteger<-function (msg, allowed) {
	repeat {
		allowMsg <- paste("[", paste(allowed, collapse = "/"), "] ", sep = "")
		outMsg <- paste(msg, allowMsg)
		cat(outMsg) 
		ans <- as.numeric(readLines(n = 1))
		if (is.na(ans)) {cat("The answer should be numeric, try again.\n")}
		else {if (ans %in% allowed) 
			{break}
			else {cat(paste(ans, "is not a valid response, try again.\n"))}}
	}
	return(ans)
}


# Function asking the user to choose a value in an interval
userQueryNB<-function (msg, allowed.min, allowed.max) {
	repeat {
		allowMsg <- paste("[min=", allowed.min, "; max=", allowed.max, "] ", sep = "")
		outMsg <- paste(msg, allowMsg)
		cat(outMsg) 
		ans <- as.numeric(readLines(n = 1))
		if (is.na(ans)) {cat("The answer should be numeric, try again.\n")}
		else {if (ans>=allowed.min & ans <=allowed.max) 
			{break}
			else {cat(paste(ans, "is not a valid response, try again.\n"))}}
	}
	return(ans)
}


# Function used to perform one iteration of the EM algo
EM.algo <- function(x,m,v,p)
{
 	n.dsn <- length(m)
 	n.x <- length(x)
 	w <- array(0, dim=c(n.x,n.dsn))
	minv=1e-100 
 	for (i in 1:n.dsn){
 		w[,i] <- p[i]*dnorm(x,m[i],sqrt(v[i]))
 	}
 	wt <- apply(w, 1, sum)
 	wt <- array(wt, dim=c(n.x,n.dsn))
 	w <- w/wt
 	Q <- 0
 	for (i in 1:n.dsn){
 		Ind <- (p[i]*dnorm(x,m[i],sqrt(v[i]))!=0)
 		Q <- Q +
		sum(w[Ind,i]*log(p[i]*dnorm(x[Ind],m[i],sqrt(v[i]))))
 	}
 	
 	m0 <- (x%*%w)/(apply(w, 2,sum))
 	for (i in 1:n.dsn){
 		v0[i]<- ((x-m0[i])^2%*%w[,i])/sum(w[,i])
 	}
 	minv0<-min(v0)
	
	if (minv0==0) {
		a=which.min(v0)
		v0[a]=minv
		minv0=minv
	}
	
	p0 <- apply(w,2,mean)
 	para <- list(m0=m0,v0=v0,p0=p0,Q=Q,w=w)
	
 	return(para)
}	


# Function initiating the EM algo and testing at which iteration the algo should stop
EM <- function(x,m0,v0,p0,eps, n)
{	
 	n0 <- 0
 	n.dsn <- length(m0)
 	n.x <- length(x)
 	f0 <- 0
 	for (i in 1:n.dsn){
 		Ind <- (p0[i]*dnorm(x,m0[i],sqrt(v0[i]))!=0)
 		f0 <- f0 + sum(log(p0[i]*dnorm(x[Ind],m0[i],sqrt(v0[i]))))
 	}
 	w0 <- NULL
 	para0 <- list(m0=m0,v0=v0,p0=p0)
 	repeat{
 		n0 <- n0+1
 		iter <- EM.algo(x,m=para0$m0,v=para0$v0,p=para0$p0) 
 		para <- list(m0=iter$m0, v0=iter$v0, p0=iter$p0)
 		f <- iter$Q;
 		w0 <- iter$w
		
 		if(n0>n | abs(f-f0)<eps | f<f0){break;}
 		f0 <- f
 		para0 <- para
 	}
 	cat("Iter: ",n0,"\n m0:",para0$m0,"\n v0:",para0$v0,"\n p0:",para0$p0,"\n f:",f0, "\n", fill=T)
	
 	return(list(para0=para0,f0=f0,w0=w0))
}


# Function plotting the the histogram and the density of your data and superimpose the components and the sum of the components
# Used to check if the EM algo identified the best components
col=c("orange", "black", "blue", "pink", "lightblue", "brown", "purple", "darkgreen", "violet")
multi.compo.plot<-function(y, x, x.seq) {
	hist(x, freq=FALSE, breaks=100, main="", xlab="", xlim=range(x), ylim=c(0,1))
	lines(density(x), col="red", xlim=range(x), ylim=c(0,1))
	par(new=TRUE)
	plot(x.seq,apply(y,1,sum),type="l",col="green", xlim=range(x), ylim=c(0,1), xlab="", ylab="")
	
	for (j in 1:compo) {
		par(new=TRUE)
		plot(x.seq, y[,j],type="l",lwd=2,col=col[j], xlim=range(x), ylim=c(0,1), xlab="", ylab="")
	}
	
	par(new=FALSE)
}


#======================== How I run the functions ==============================================
# Let call your.data the variable for which you want to identify the components
# First you have to initiate the EM algo by giving the number of components you want,
# their mean, variance and proportion.
# I am using the distribution of the data (histogram) to estimate these numbers as followed.

compo=as.vector(NULL)
m.ini=as.vector(NULL)
v.ini=as.vector(NULL)
p.ini=as.vector(NULL)

for (i in 1) {
hist(your.data, freq=FALSE, breaks=50)
lines(density(your.data), col="red")
compo=c(compo, userQueryInteger(msg="Please evaluate the number of components from the graph", allowed=seq(1,9, by=1)))
	
for (j in 1:compo[i]) {
	m.ini=c(m.ini, userQueryNB(msg=paste("Please evaluate the logratio at the peak of component no ", j," from the graph",sep=""), allowed.min=-3, allowed.max=3))
	v.ini=c(v.ini, userQueryNB(msg=paste("Please evaluate the variance of component no ", j," from the graph",sep=""), allowed.min=0, allowed.max=1))
	p.ini=c(p.ini, userQueryNB(msg=paste("Please evaluate the proportion of component no ", j," from the graph (sum of the proportion of the ", compo[i], " components have to be exactly 1)",sep=""), allowed.min=0, allowed.max=1))
}
repeat{
	if (sum(p.ini)!=1) {
		som=sum(p.ini)
		p.ini=0
		cat(paste("Sum of the proportions is ", som,". It is not exactly equal to 1, try again.\n", sep=""))
		
		for (j in 1:compo[i]) {
			p.ini=c(p.ini, userQueryNB(msg=paste("Please evaluate the proportion of component no ", j," from the graph (sum of the proportion of the ", compo[i], " components have to be exactly 1)",sep=""), allowed.min=0, allowed.max=1))
		}
		
		p.ini=p.ini[-1]
	} else {break}
}
}

#==================================================================================
# After initiating the EM algo, you can run it and check if the components found are correct.   

v0=v.ini
your.data.EM  <- EM(x=your.data, m0=m.ini, v0=v.ini, p0=p.ini, eps=.1e-5, n=1000) # n = max number of iterations; eps = cutoff to stop the iterations (I have never changed it)
m1=c(your.data.EM$para0$m0)
v1=c(your.data.EM$para0$v0)
p1=c(your.data.EM$para0$p0)
your.data.seq=seq(min(your.data),max(your.data),length=500)
y=rep(0, times=500)
	
for (i in 1:compo) {
	y=data.frame(y, p1[i]*dnorm(your.data.seq,mean=m1[i], sd=sqrt(v1[i])))
}
y=y[,-1]
	
a=multi.compo.plot(y, your.data, your.data.seq)

your.data=data.frame(your.data)
your.data[,c(paste("proba.compo", seq(1,compo, by=1),sep=""))]=your.data.EM$w0[,c(1:compo)]
head(your.data)

	