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
