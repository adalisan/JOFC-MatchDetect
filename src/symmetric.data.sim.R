# TODO: Add comment
# 
# Author: Sancar
###############################################################################
require(MASS)
require(MCMCpack)
require(shapes)
require(vegan)
require(klaR)
setwd("./DataFusion")
source("./lib/simulation_math_util_fn.R")

r<-10000
simplex.dim<-4
d<-simplex.dim-1
t<-100
Wchoice <-"avg"
oos <- T
w.vals<-c(0.1,0.4,0.5,0.8,0.95,0.99)


psuedo.counts<- rep(t,simplex.dim)
normalize.rows<-function(X){
	for (i in 1:nrow(X)){
		X[i,]<- X[i,]/sum(X[i,])
	}	
	return(X)
}
displace<-diag(simplex.dim)*3*t/4
coord.matrix <-matrix(rep(psuedo.counts,simplex.dim),simplex.dim,simplex.dim)+displace

coord.matrix<- normalize.rows(coord.matrix)

#print(coord.matrix)


sym.orig <- rep(1/simplex.dim,simplex.dim)
nmc<- 100

agg.test.stats.null <- c()
agg.test.stats.alt <- c()
agg.test.stats.null.true.proj<-c()
agg.test.stats.null.iden.proj<-c()
agg.test.stats.alt.iden.proj<- c()
agg.test.stats.null.hub.pt <- c()
agg.test.stats.null.fringe.pt <- c()
agg.test.stats.alt.hub.pt <- c()
agg.test.stats.null.JOFC <- c()
agg.test.stats.alt.JOFC <- c()

nonid.Proc<-0
id.Proc <-TRUE

Bary.to.Euc <- matrix(
		c(1,1,1,
				-1,-1,1,
				-1,1,-1,
				1,-1,-1),simplex.dim,simplex.dim,byrow=FALSE)
use.Euc.coords <-FALSE
corr.Bary.coords<- TRUE

run.jofc.analysis<-FALSE


windows()
discarded.dims.vec.1<-rep(0,nmc)

discarded.dims.vec.2<-rep(0,nmc)
init.proj.plot<-dev.cur()
windows()
second.proj.plot<-dev.cur()
for (mc in 1:nmc){
	
	dirich.data<- matched_rdirichlet(n=simplex.dim, p=simplex.dim-1, r=r, q=3, c=0, alpha=coord.matrix)
	X1<- dirich.data$X1
	X2<- dirich.data$X2
	dev.set(init.proj.plot)
#	s3d.1 <-quadplot(e=X1)
#	quadlines(centerlines(4), sp = s3d.1, lty = "dashed")
#	s3d.2 <-quadplot(e=X2)
#	quadlines(centerlines(4), sp = s3d.2, lty = "dashed")
	
	
	# Let PCA find the exact prin. components for now.
	col.vars<-apply(X1,2,var)
	discarded.dim.1<-which.min(col.vars)
	discarded.dims.vec.1[mc] <- discarded.dim.1
	col.vars<-apply(X2,2,var)
	discarded.dim.2<-which.min(col.vars)
	discarded.dims.vec.2[mc] <- discarded.dim.2

	
	#pr.1 <-prcomp(X1)
	#pr.2 <-prcomp(X2)
	if (use.Euc.coords){
	X1.b <- X1
	X1.b[,discarded.dim.1]<-0
	X2.b <- X2
	X2.b[,discarded.dim.2]<-0
	X1t<- X1.b%*%t(Bary.to.Euc)
	X2t<- X2.b%*%t(Bary.to.Euc)
    } else{
		#X1t<- X1[,-discarded.dim.1]
		#X2t<- X2[,-discarded.dim.2]
		X1t<- X1
		X1t[,discarded.dim.1]<-X1t[,simplex.dim]
		X1t<-X1t[,1:(simplex.dim-1)]
		X2t<- X2
		X2t[,discarded.dim.2]<-X2t[,simplex.dim]
		X2t<-X2t[,1:(simplex.dim-1)]
		
		
		X1t <- normalize.rows(X1t)
		X2t <- normalize.rows(X2t)
		
	}
	if (mc==1) new.plot<-triplot(X1t)
	s2d.1<-tripoints(X1t,col="red")
	s2d.2<-tripoints(X2t,col="blue")
	
	for (k in 1:simplex.dim){
		trilines(rbind(X1t[k,],X2t[k,]),col="black")
	}
	#tripoints(s2d.1)
	#tripoints(s2d.2)
	#X1t.cent<- X1t-matrix(rep(colMeans(X1t),nrow(X1t)),
	#		nrow=nrow(X1t),ncol=ncol(X1t),byrow=TRUE)
	#X2t.cent<- X2t-matrix(rep(colMeans(X2t),nrow(X2t)),
	#		nrow=nrow(X2t),ncol=ncol(X2t),byrow=TRUE)
	X1t.cent<- X1t - matrix(rep(1/d,nrow(X1t)*ncol(X1t)),
			nrow=nrow(X1t),ncol=ncol(X1t),byrow=TRUE)
X2t.cent<- X2t - matrix(rep(1/d,nrow(X2t)*ncol(X2t)),
			nrow=nrow(X2t),ncol=ncol(X2t),byrow=TRUE)
		

pr.dim<-dim(X1t)[2]-1
	
	
	Proc.T.2<-list(R=diag(pr.dim),s=1)
	# Match 1st condition to
	if (!id.Proc){
	if (use.Euc.coords){
		Proc.T.2<-procOPA(X1t,X2t,scale=FALSE,reflect=FALSE)	
	} else if (corr.Bary.coords){
		Proc.T.2<-procOPA(X1t[,1:pr.dim],X2t[,1:pr.dim],scale=FALSE,reflect=FALSE)
	} else{
		Proc.T.2<-procOPA(X1t,X2t,scale=FALSE,reflect=FALSE)	

	}
} else {
	
		Proc.T.2$R<-diag(dim(Proc.T.2$R)[1])
		Proc.T.2$s <-1
	}
	
	print (det(Proc.T.2$R))
#	Proc.T.2$R<-diag(pr.dim)
#	Proc.T.2$s<-1
	#Proc.T<-procrustes(X2t[,],X1t[,],dilation=FALSE)
	
#Now the test statistic is computed for both matched points
# and the unmatched points
# i-th point in first condition(X) vs i-th point in second condition(Y)
# i-th point in first condition(X) vs all other points in second condition(Y)
dev.set(second.proj.plot)	

for (i in 1:simplex.dim){
		
		train.data.X <- X1t.cent[-i,]
		train.data.Y <- X2t.cent[-i,]
		test.data.X <-  X1t.cent[i,]
		test.data.Y <-  X2t.cent[i,]
#		alt.pt.ind<-sample(1:(simplex.dim-1),1)
#		alt.pt.ind.i <- alt.pt.ind
#		if (i<=alt.pt.ind)
#			alt.pt.ind.i<- alt.pt.ind.i+1
#		alt.pt<-train.data.Y [alt.pt.ind,]
		alt.pt<-train.data.Y
#		
#		
		#Proc.T<-procrustes(train.data.Y,train.data.X,dilation=TRUE)
		#Cheat

			print(mc)
#			print(Proc.T$R)
#			print(Proc.T$s)
#			print(Proc.T$tt)
#			
			#Proc.T$R<-diag(simplex.dim-1)
		#Proc.T$s<-1
		if (use.Euc.coords){		
		test.stats.null <-sum((test.data.X -test.data.Y%*%Proc.T.2$R*Proc.T.2$s)^2) 
		test.stats.alt  <-rowSums((matrix(rep(test.data.X,dim(alt.pt)[1]),
									dim(alt.pt)[1],dim(alt.pt)[2],byrow=TRUE) -alt.pt%*%Proc.T.2$R*Proc.T.2$s)^2)
	    } else if (corr.Bary.coords) {
			
		
#		
#	
    	test.data.Y.R    <- test.data.Y[1:(pr.dim)]%*%Proc.T.2$R*Proc.T.2$s
		test.data.Y.R    <- test.data.Y.R +rep(1/d,pr.dim)
		test.data.Y.bary <- c(test.data.Y.R,1-sum(test.data.Y.R))
		
		#Proc.T.2$R<- diag(pr.dim)
		alt.test.data.Y.R <- alt.pt[,1:(pr.dim)]%*%Proc.T.2$R*Proc.T.2$s
		alt.test.data.Y.R <- alt.test.data.Y.R +matrix(1/d,simplex.dim-1,pr.dim)
		alt.test.data.Y.bary <- cbind(alt.test.data.Y.R ,1-rowSums(alt.test.data.Y.R ))
		
		test.data.X.bary <- test.data.X +rep(1/d,d)
		test.stats.null <-sum((test.data.X.bary - test.data.Y.bary)^2)
		if (mc==1) new.plot<-triplot(test.data.X.bary)
		tripoints(test.data.X.bary,col="red",pch=3)
		tripoints(test.data.Y.bary,col="blue",pch=3)
		line.btw.matched <- rbind(test.data.X.bary,test.data.Y.bary)
		print(dim(line.btw.matched))
		trilines(line.btw.matched,col="black")
		test.stats.alt  <-rowSums((test.data.X.bary- alt.test.data.Y.bary)^2)
	} else{
	}
		
#		test.stats.null <-sum((test.data.X -test.data.Y)^2) 
#		test.stats.alt  <-sum((test.data.X -alt.pt)^2)
		
		agg.test.stats.null <- c(agg.test.stats.null,test.stats.null)
		agg.test.stats.alt <- c(agg.test.stats.alt,test.stats.alt)
		if (discarded.dim.1==discarded.dim.2){
		  agg.test.stats.null.iden.proj<-c(agg.test.stats.null.iden.proj,test.stats.null)
		  agg.test.stats.alt.iden.proj<-c(agg.test.stats.alt.iden.proj,test.stats.alt)
	     }
		 if ((i==discarded.dim.1)|(i==discarded.dim.2))
		   agg.test.stats.null.hub.pt<-c(agg.test.stats.null.hub.pt,test.stats.null)
	      if ((i!=discarded.dim.1)&(i!=discarded.dim.2))
		   agg.test.stats.null.fringe.pt<-c(agg.test.stats.null.fringe.pt,test.stats.null)
	   
#		 if ((i==discarded.dim.1)|(alt.pt.ind.i==discarded.dim.2)){
#			 
#			 agg.test.stats.alt.hub.pt<-c(agg.test.stats.alt.hub.pt,test.stats.alt)
#		 }
			 
		if (discarded.dim.1!=discarded.dim.2){
			if ((i!=discarded.dim.1)&&(i!=discarded.dim.2))
			agg.test.stats.null.true.proj <- c(agg.test.stats.null.true.proj,test.stats.null)
		}
		
	}	
	
	if (run.jofc.analysis){

num.test.pts <- 100
train.reps<-50
	for (i in 1:simplex.dim){
	train.coord.matrix <-coord.matrix[-i,]
      X1<-matrix()
	X2<-matrix()
	for (j in 1:train.reps){
	train.dirich.data<- matched_rdirichlet(n=(simplex.dim-1), p=simplex.dim-1,
		 r=r, q=3, c=0, alpha=train.coord.matrix)
	X1<- rbind(X1,train.dirich.data$X1)
	X2<- rbind(X1,train.dirich.data$X2)
      }

#	alpha.matched<-rdirichlet(rep(1,4))
      alpha.matched <- coord.matrix[i,]
	test.dirich.data<- matched_rdirichlet(n=num.test.pts, p=simplex.dim-1,
		 r=r, q=3, c=0, alpha=alpha.matched)
	alpha.unmatched<-rdirichlet(num.test.pts,rep(1,4))
	test.dirich.data.alt<- matched_rdirichlet(n=num.test.pts, p=simplex.dim-1,
		 r=r, q=3, c=0, alpha=alpha.unmatched)
	
		Y1 <- test.dirich.data$X1
		Y20 <- test.dirich.data$X2
		alt.pt.index <- sample(1:(simplex.dim-1),1)

		Y2A <- test.dirich.data.alt$X2
		if (oos){
			D1<- dist(X1)
			D2<- dist(X2)
			D10<-dist(rbind(X1,Y1))
			D20<-dist(rbind(X2,Y20))
			D2A<-dist(rbind(X2,Y2A))
			D.oos.1<-dist(Y1)
			D.oos.2.null<-dist(Y20)
			D.oos.2.alt<-dist(Y2A)
		}
		else{
			D1<- dist(X1)
			D2<- dist(X2)
			
			
		}
		if (Wchoice == "avg") {
			W <- (D1 + D2)/2
		} else if (Wchoice == "sqrt") {
			W <- sqrt((D1^2 + D2^2)/2)
		}
		M <- omnibusM(D1, D2, W)
		
		ideal.omnibus.0  <- as.matrix(dist(rbind(X1,X2,Y1,Y20)))
		ideal.omnibus.A  <- as.matrix(dist(rbind(X1,X2,Y1,Y2A)))
		
		
		X <- cmdscale(M, d-1)
		
		
		## ==== jofc ====
		if (oos){
		run.jofc(D1,D2,D10,D20,D2A,D.oos.1,D.oos.2.null,D.oos.2.alt, ideal.omnibus.0,ideal.omnibus.A,
				n=train.reps*(simplex.dim-1),m=num.test.pts,
				d=(simplex.dim-1),c.val=0,
				model="dirichlet",oos=TRUE,Wchoice=Wchoice,
				separability.entries.w=FALSE,
				wt.equalize=FALSE,assume.matched.for.oos=TRUE,
				oos.use.imputed=FALSE,
				compare.pom.cca = FALSE,
				w.vals=w.vals,
				verbose=TRUE)
		
				
				)
		
		Y1t  <- oosMDS(dist(rbind(X1, Y1)), X1t)
		Y20t <- oosMDS(dist(rbind(X2, Y20)), X2t)
		Y2At <- oosMDS(dist(rbind(X2, Y2A)), X2t)
	   }
	   else{
		   Y1t<- X[i,]
		   
		   Y20t<- X[simplex.dim+i,]
		   train.X <- X[-c(i,simplex.dim+i),]
		   Y2At <- train.X[alt.pt.index+simplex.dim-1,]
#		   alt.pt.index.X2 <- simplex.dim+alt.pt.index
#		   if (i<=alt.pt.index)
#				alt.pt.index.X2<- alt.pt.index.X2 +1
#		   Y2At<- X[alt.pt.index.X2,]
		   #print(rbind(Y1t,Y20t,Y2At))
	   }
		T0 <- sum((Y1t - Y20t)^2)
		TA <- sum((Y1t - Y2At)^2)
		agg.test.stats.null.JOFC <- c(agg.test.stats.null.JOFC,T0)
		agg.test.stats.alt.JOFC  <- c(agg.test.stats.alt.JOFC,TA)
		
	}
	}
	
#	
#	if (oos == TRUE) {
#		
#	} else {
#		M0 <- omnibusM(D10A, D20, W=(D10A+D20)/2)
#		MA <- omnibusM(D10A, D2A, W=(D10A+D2A)/2)
#		X0 <- cmdscale(M0, d)
#		XA <- cmdscale(MA, d)
#		T0 <- rowSums((X0[(n+1):(n+m), ] - X0[(n+m+n+1):(n+m+n+m), ])^2)
#		TA <- rowSums((XA[(n+1):(n+m), ] - XA[(n+m+n+1):(n+m+n+m), ])^2)
#	}
	
}

a.obj<-hist2d(discarded.dims.vec.1,discarded.dims.vec.2)
 a.obj$counts[ a.obj$counts>0]

 sum(a.obj$counts[ a.obj$counts>0])



plot(discarded.dims.vec.1, discarded.dims.vec.2)
windows()
hist(discarded.dims.vec.1)
windows()
hist(discarded.dims.vec.2)

x.max<- max(c(agg.test.stats.null,agg.test.stats.alt))
h.breaks<-seq(0,x.max,x.max/40)
windows()
hist(agg.test.stats.null,breaks=h.breaks,col="red")
# Null statistic non-zero if it's hub point and projected to two different faces
# OR If it's a fringe point, and specifically point 4. If it's a fringe point,
# it is not a hub point in both of the projections and as such it has to be point 4,
#since it is not coincident in both projections all other points are coincident
# when they're not the hub point.
# So the probability of choosing number 4 is (exactly, since we're taking all of the points
#as test points, in turn) 1/4. And the probability of choosing two  different projections
#(if they're the same projection the test statistic is 0) (3/4)*(2/4). Expected number of
# counts of null statistic>Large is nmc*3/32




hist(agg.test.stats.alt,breaks=h.breaks,add=TRUE,col="blue")hist(agg.test.stats.null.hub.pt,breaks=h.breaks,add=TRUE,col="purple")
hist(agg.test.stats.null.fringe.pt,breaks=h.breaks,add=TRUE,col="black")
hist(agg.test.stats.null.iden.proj,breaks=h.breaks,add=TRUE,col="darkkhaki")
hist(agg.test.stats.alt.iden.proj,breaks=h.breaks,add=TRUE,col="khaki")
windows()

plot( density(agg.test.stats.null),
		xlim=range( c(agg.test.stats.null,agg.test.stats.alt,
						agg.test.stats.null.true.proj,agg.test.stats.null.iden.proj,
						agg.test.stats.alt.iden.proj,agg.test.stats.alt.hub.pt) ),
		col="red",main="PoM", xlab="" )
lines(density(agg.test.stats.alt), col="blue")
#lines(density(agg.test.stats.null.true.proj), col="green")
#lines(density(agg.test.stats.null.iden.proj), col="pink")
lines(density(agg.test.stats.alt.iden.proj), col="brown")
#lines(density(agg.test.stats.alt.hub.pt), col="orange")
lines(density(agg.test.stats.null.hub.pt), col="purple")
lines(density(agg.test.stats.null.fringe.pt), col="black")

windows()

x.max<- max(c(agg.test.stats.null.JOFC,agg.test.stats.alt.JOFC))
h.breaks<-seq(0,x.max,x.max/40)
hist(agg.test.stats.null.JOFC,breaks=h.breaks,col="red")
hist(agg.test.stats.alt.JOFC,breaks=h.breaks,add=TRUE,col="blue")
windows()

plot( density(agg.test.stats.null.JOFC),
		xlim=range( c(agg.test.stats.null.JOFC,agg.test.stats.alt.JOFC
						) ),
		col="red",main="JOFC stat null.vs.alt", xlab="" )
lines(density(agg.test.stats.alt.JOFC), col="blue")
