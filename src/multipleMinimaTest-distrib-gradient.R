
results.dir<-"./graphs"
meshgrid <- function(a,b) {
	list(
			x=outer(b*0,a,FUN="+"),
			y=outer(b,a*0,FUN="+")
	)
} 
oos.copies<-20
oos.n<-2*oos.copies
n<-5+oos.n
X<-matrix(0,n,2)
X[1:4,]<-cbind(c(0,1,1,0),c(0,0,1,1))

X[5,]<-c(1,0)
X[6,]<-c(0,1)
X[n,]<-c(0.5,0.5)
in.sample.ind<-c(1:4,n)#,8)

indices.test.cluster.5<-seq(5,n-1,2)
indices.test.cluster.6<-seq(6,n-1,2)

X[indices.test.cluster.5,]<-mvrnorm(oos.copies,mu=c(X[5,]),Sigma=0.001*diag(2))

X[indices.test.cluster.6,]<-mvrnorm(oos.copies,mu=c(X[6,]),Sigma=0.001*diag(2))

d.X<- dist(X)


w <- 0.99



i <- 4

for (j in indices.test.cluster.5){
	j <- 5
	
	ind <- (n*(i-1) - i*(i-1)/2 + j-i)
	
	d.X[ind]<-d.X[ind]-1.4
}

for (j in indices.test.cluster.6){
	i<-2
	j<-6
	ind<-(n*(i-1) - i*(i-1)/2 + j-i)
	
	d.X[ind]<- d.X[ind]-1.4
#X[5,]<-c(0,1)
}

i<-5
j<-7
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

#d.X[ind]<- d.X[ind]+0.707

i<-5
j<-8
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

d.X[d.X<0]=0;
#d.X[ind]<- d.X[ind]-1.4

#print(d.X)

init.config<-X
w.vals<-c(0.1,0.45,0.5,0.55,0.99)
w.vals<-c(0.1,0.45,0.5,0.55,0.75,0.9,0.99)


CDF.T0<-array(0,dim=c(length(w.vals),3))
CDF.TA<-array(0,dim=c(length(w.vals),3))
ecdf.T0<-list()
ecdf.TA<-list()


#w.vals<-c(0.45,0.9)
min.config.stress.1.w<- rep(0, length(w.vals))
min.config.stress.2.w<- rep(0, length(w.vals))

  grid.seq.x<-seq(-0.5,1.5,0.1)
	grid.seq.y<-seq(-0.2,1.6,0.1)

T0.1 <-array(0, dim= c(length(grid.seq.x),length(grid.seq.y),oos.copies))
TA.1 <-array(0, dim= c(length(grid.seq.x),length(grid.seq.y),oos.copies))
T0.2 <-array(0, dim= c(length(grid.seq.x),length(grid.seq.y),oos.copies))
TA.2 <-array(0, dim= c(length(grid.seq.x),length(grid.seq.y),oos.copies))
T0.w.1 <-array(0, dim= c(length(w.vals),length(grid.seq.x)*length(grid.seq.y)*oos.copies))
TA.w.1 <-array(0, dim= c(length(w.vals),length(grid.seq.x)*length(grid.seq.y)*oos.copies))
T0.w.2 <-array(0, dim= c(length(w.vals),length(grid.seq.x)*length(grid.seq.y)*oos.copies))
TA.w.2 <-array(0, dim= c(length(w.vals),length(grid.seq.x)*length(grid.seq.y)*oos.copies))




for (w.i in 1:length(w.vals)){
	w <- w.vals[w.i]
	
	
	
	
	W <- matrix(1-w,n,n)
	for (o in indices.test.cluster.5)
		W[o,2] <- W[2,o]<- w
	for (o in indices.test.cluster.6)
		W[o,4] <- W[4,o]<- w
#W[5,7] <- W[7,5]<- 1-w
	diag(W)<-0
	
	
	W.oos <- W
	W.oos[in.sample.ind, in.sample.ind]<-0	
	new.index.order<- c(in.sample.ind,(5:(n-1)))
	W.oos<-W.oos[new.index.order,new.index.order]
	
	
	
	in.sample.Bool<-(1:n %in% in.sample.ind)
	X.embed.1.in <- smacofM(as.matrix(d.X)[in.sample.ind,in.sample.ind],
			ndim    = 2,
			W       = W[in.sample.ind,in.sample.ind],
			init    = init.config[in.sample.Bool,]		,
			verbose = FALSE,
			itmax   = 1000,
			eps     = 1e-6)
	
	X.embed.1.oos <- oosIM(D=as.matrix(d.X),
			X=X.embed.1.in,
			init= init.config[!in.sample.Bool,],
			W=W.oos,
			verbose = FALSE,
			itmax   = 1000,
			eps     = 1e-6,
			bwOos = TRUE,
			isWithin = ifelse(in.sample.Bool,1,0) )
	
	row.names(X.embed.1.oos)<-NULL
	row.names(X.embed.1.in)<-NULL
	X.embed.1.norm <- rbind(X.embed.1.in[1:4,] ,X.embed.1.oos,X.embed.1.in[5,]) 
	print("Init config ")
	print(X.embed.1.norm)
	stress <- sum(as.dist(W)*((dist(X.embed.1.norm)-d.X)^2))
	print("Init config  stress")
	print(stress)
	close.to.init.1<-data.frame(x=numeric(),y=numeric())
	close.to.init.2<-data.frame(x=numeric(),y=numeric())
	far.to.init.1<-data.frame(x=numeric(),y=numeric())
	far.to.init.2<-data.frame(x=numeric(),y=numeric())
	min.stress.1<-100
	min.stress.2<-100

	
	stress.at.loc<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
	
	grid.resp<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
	grid.coords<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
	
	mesh.grid.coords<-meshgrid(grid.seq.x,grid.seq.y)
	
	final.coords.x.5<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
	final.coords.y.5<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
	final.coords.x.6<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
	final.coords.y.6<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
	min.config.1<-matrix(0,n,2)
	min.config.2<-matrix(0,n,2)
	
	for (j in 1:length(grid.seq.y)){
		
		for (i in 1:length(grid.seq.x)){
			init.config <- X
			
			i.x<-grid.seq.x[i]
			j.y<-grid.seq.y[j]
			#i.x<-mesh.grid.coords$x[i,j] 
			#j.y<-mesh.grid.coords$y[i,j] 
			
			for (l in indices.test.cluster.5){
				init.config[l,]<-c(i.x,j.y)
			}
			for (l in indices.test.cluster.6){
				init.config[l,]<-c(1-i.x,1-j.y)
			}
			X.embed.2.in<-smacofM(as.matrix(d.X)[in.sample.ind, in.sample.ind],
					ndim    = 2,
					W       = W[in.sample.ind,in.sample.ind],
					init    = init.config[in.sample.Bool,],
					verbose = FALSE,
					itmax   = 1000,
					eps     = 1e-6)
			
			X.embed.2.oos<-oosIM(D = as.matrix(d.X),
					X = X.embed.2.in,
					init = init.config[!in.sample.Bool,],
					
					verbose =  FALSE,#(i<=(length(grid.seq.x)/2)),
					itmax   = 1000,
					eps     = 1e-6,
					W = W.oos,
					bwOos = FALSE,
					isWithin = ifelse(in.sample.Bool,1,0))
			
			row.names(X.embed.2.oos)<-NULL
			row.names(X.embed.2.in)<-NULL
			
			X.embed.2.norm<- rbind(X.embed.2.in[1:4,] ,X.embed.2.oos,X.embed.2.in[5,])		
			
			final.coords.x.5[i,j]<-X.embed.2.norm[5,1]
			final.coords.y.5[i,j]<-X.embed.2.norm[5,2]
			final.coords.x.6[i,j]<-X.embed.2.norm[6,1]
			final.coords.y.6[i,j]<-X.embed.2.norm[6,2]
			#print("X.embed.2.norm")
			#print(X.embed.2.norm)
			stress <- sum(as.dist(W)*((dist(X.embed.2.norm)-d.X)^2))
			stress.unif.wt <- sum(((dist(X.embed.2.norm)-d.X)^2))
			stress.unif.abs <- sum(abs(dist(X.embed.2.norm)-d.X))
			#print(stress)
			#	print(stress.unif.wt)
			#	print(stress.unif.abs)		
			if (((X.embed.2.norm[6,2]>0.5)) & 
					(((X.embed.2.norm[5,1]>X.embed.2.norm[6,1])))) { 
				#& (sum(X.embed.2.norm[6,]<0.3)==1))
				
				
				grid.resp[i,j]<-1
				#	print("First Test Point")
				#	print(X.embed.2.norm[5,])
				if (stress < min.stress.1) {
					min.stress.1<-stress
					print("Min stress found(real min)")
					print(X.embed.2.norm)
					min.config.1<-X.embed.2.norm
					print(min.stress.1)
					
				}
				close.to.init.1<-rbind(close.to.init.1,X.embed.2.norm[5,])
				close.to.init.2<-rbind(close.to.init.2,X.embed.2.norm[6,])
				
				
			}else{
				grid.resp[i,j]<-0
				#print("First Test Point")
				print(X.embed.2.norm[5,])		
				if (stress < min.stress.2){
					min.stress.2<-stress
					print("Min stress found(second min)")
					print(X.embed.2.norm)
					print(min.stress.2)
					min.config.2<-X.embed.2.norm
					
				}
				far.to.init.1<-rbind(far.to.init.1,X.embed.2.norm[5,])
				far.to.init.2<-rbind(far.to.init.2,X.embed.2.norm[6,])
				
				#} else{
				#	grid.resp[i,j] <- NA
			}
			
			stress.at.loc[i,j]<- stress
			T0.1[i,j,]<- sqrt(rowSums((X.embed.2.norm[indices.test.cluster.5,]-matrix(X.embed.2.norm[2,],oos.copies,2,byrow=TRUE))^2))
			T0.2[i,j,]<- sqrt(rowSums((X.embed.2.norm[indices.test.cluster.6,]-matrix(X.embed.2.norm[4,],oos.copies,2,byrow=TRUE))^2))
			TA.1[i,j,]<- sqrt(rowSums((X.embed.2.norm[indices.test.cluster.5,]-matrix(X.embed.2.norm[4,],oos.copies,2,byrow=TRUE))^2))
			TA.2[i,j,]<- sqrt(rowSums((X.embed.2.norm[indices.test.cluster.6,]-matrix(X.embed.2.norm[2,],oos.copies,2,byrow=TRUE))^2))
			
			
			
		}
		
	}
	
	T0.w.1[w.i,]<-as.vector(T0.1)
	T0.w.2[w.i,]<-as.vector(T0.2)
	TA.w.1[w.i,]<-as.vector(TA.1)
	TA.w.2[w.i,]<-as.vector(TA.2)
}



for (w.i in 1:length(w.vals)){
frac.1 <-sum(T0.w.1[w.i,]<0.45)/(dim(T0.w.1)[2])
frac.2 <-sum(T0.w.1[w.i,]<0.9)/(dim(T0.w.1)[2])-frac.1
frac.3<- 1-frac.1-frac.2
CDF.T0[w.i,]<- c(frac.1,frac.2,frac.3)
ecdf.T0 <- c(ecdf.T0,list(ecdf(T0.w.1[w.i,])))

frac.1 <-sum(TA.w.1[w.i,]<0.45)/(dim(TA.w.1)[2])
frac.2 <-sum(TA.w.1[w.i,]<0.9)/(dim(TA.w.1)[2])-frac.1
frac.3<- 1-frac.1-frac.2
CDF.TA[w.i,]<- c(frac.1,frac.2,frac.3)
ecdf.TA <- c(ecdf.TA,list(ecdf(TA.w.1[w.i,])))
}

print("CDF.T0")
print(CDF.T0)


print("CDF.TA")
print(CDF.TA)

for (w.i in 1:length(w.vals)){
  if (run.in.linux) {X11()} else {windows()}
  par(mfrow=c(2,1))
  plot(ecdf.T0[[w.i]],main="")
  plot(ecdf.TA[[w.i]],main="")
  title(substitute(w == w.i ,list(w.i=w.vals[w.i])))
}


bw.type <- 0.05
if (run.in.linux) {X11()} else {windows()}
#par(mfrow=c(length(w.vals)/2,2))
for (w.i in 1:length(w.vals)){
	h.1<-density(T0.w.1[w.i,],bw=bw.type)
	bw.type <- h.1$bw
	plot(h.1,main="")
 title(substitute(w == w.i ,list(w.i=w.vals[w.i])))
	if (run.in.linux) {X11()} else {windows()}
}

colors.vec<-rainbow(length(w.vals))
if (run.in.linux) {X11()} else {windows()}
for (w.i in 1:length(w.vals)){
	h.1<-density(T0.w.1[w.i,],bw=bw.type)
	bw.type <- h.1$bw
	if (w.i==1) plot(h.1,type="l",col=colors.vec[w.i])
	else {lines(h.1,main="",col=colors.vec[w.i])}
  title(substitute(w == w.i ,list(w.i=w.i)))
}
legend.txt<-w.vals
legend("bottomright", legend=legend.txt, col=colors.vec, lty=rep(1, length(w.vals)))


if (run.in.linux) {X11()} else {windows()}
par(mfrow=c(length(w.vals)/2,2))
for (w.i in 1:length(w.vals)){
	
	h.1<-density(T0.w.2[w.i,],bw=bw.type)
	bw.type <- h.1$bw
	plot(h.1,xlim=c(0,3.5),ylim=c(0,0.8))
}

if (run.in.linux) {X11()} else {windows()}
par(mfrow=c(length(w.vals)/2,2))
for (w.i in 1:length(w.vals)){
	
	h.1<-density(TA.w.1[w.i,],bw=bw.type)
	bw.type <- h.1$bw
  plot(h.1,xlim=c(0,3.5),ylim=c(0,0.8))
}

if (run.in.linux) {X11()} else {windows()}
par(mfrow=c(length(w.vals)/2,2))
for (w.i in 1:length(w.vals)){
	
	h.1<-density(TA.w.2[w.i,],bw=bw.type)
	bw.type <- h.1$bw
plot(h.1,xlim=c(0,3.5),ylim=c(0,0.8))
}
