# TODO: Add comment
# 
# Author: Sancar
###############################################################################

setwd(paste(Sys.getenv("R_HOME") ,"./../projects/",collapse=""))

setwd("./DataFusion_Priebe/")
   source("./src/simulation_math_util_fn.R")
   require(MASS)
   require(plotrix)

	windows()
    r<-1.5
	p<-2
	q<-2
    n<- 8
	c<- 0
	box.dim<-2
	sigma<- r^2*diag(p)
	alpha.mc <- mvrnorm(n, rep(0,p),sigma)
	print(dim(alpha.mc))
	## n pairs of matched points
	xlist <- matched_rnorm(n, p, q, c,r, alpha=alpha.mc[1:n, ],sigma.alpha=sigma)
	X1 <- xlist$X1
	X2 <- xlist$X2
	par(mfcol=c(1,2))
	par(mfg=c(1,1,1,2))
	plot(0,0,xlim=c(-box.dim,box.dim),ylim=c(-box.dim,box.dim),asp=1)
	mapply(draw.circle,alpha.mc[,1],alpha.mc[,2],
			MoreArgs=list(radius=2/r,nv=100,border=NULL,col=NA,lty=3,lwd=1))
	points(X1,type="p",col="red",cex=0.5)
	
	points(alpha.mc,type="p",pch=20,col="black")
	text(15,15,bquote(Xi[1]))
	
	
	par(mfg=c(1,2,1,2))
	plot(0,0,xlim=c(-box.dim,7),ylim=c(-box.dim,7),asp=1)
	mapply(draw.circle,alpha.mc[,1],alpha.mc[,2],
			MoreArgs=list(radius=r,nv=100,border=NULL,col=NA,lty=3,lwd=1))
	points(X2,type="p",col="blue",,cex=0.5)
	points(alpha.mc,type="p",pch=20,col="black")
	text(15,15,bquote(Xi[2]))
