
meshgrid <- function(a,b) {
  list(
       x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+")
       )
} 

n<-7
X<-matrix(0,n,2)
X[,1]<-c(0,1,1,0,0,0,0.5)#1.01)
X[,2]<-c(0,0,1,1,0,0,0.5)#,0)
X[5,]<-c(1,0)
X[6,]<-c(0,1)
#X[8,]<-c(0,1.0)
in.sample.ind<-c(1:4,7)#,8)
d.X<- dist(X)


w <- 0.99



i <- 4
j <- 5

ind <- (n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<-d.X[ind]-1.4

i<-2
j<-6
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<- d.X[ind]-1.4
#X[5,]<-c(0,1)


i<-5
j<-7
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

#d.X[ind]<- d.X[ind]+0.707

i<-5
j<-8
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

#d.X[ind]<- d.X[ind]-1.4

print(d.X)

init.config<-X
w.vals<-c(0.1,0.45,0.5,0.55,0.99)

min.config.stress.1.w<- rep(0, length(w.vals))
min.config.stress.2.w<- rep(0, length(w.vals))

for (w.i in 1:length(w.vals)){
w <- w.vals[w.i]



W <- matrix(1-w,n,n)
W[5,2] <- W[2,5]<- w
W[6,4] <- W[4,6]<- w
W[5,7] <- W[7,5]<- 1-w
diag(W)<-0


W.oos <- W
W.oos[in.sample.ind, in.sample.ind]<-0	
new.index.order<- c(in.sample.ind,c(5,6))
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
X.embed.1.norm <- rbind(X.embed.1.in[1:4,] ,X.embed.1.oos,X.embed.1.in[5:(n-2),]) 
print("Init config ")
print(X.embed.1.norm)
  stress <- sum(as.dist(W)*((dist(X.embed.1.norm)-d.X)^2))
print("Init config  stress")
print(stress)
close.to.init.1<-c(1.5,1)
close.to.init.2<-c(1.5,1)
far.to.init.1<-c(1.5,1)
far.to.init.2<-c(1.5,1)
min.stress.1<-100
min.stress.2<-100
grid.seq.x<-seq(-0.5,1.5,0.1)
grid.seq.y<-seq(-0.2,1.6,0.1)

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

		init.config[5,]<-c(i.x,j.y)
		init.config[6,]<-c(1-i.x,1-j.y)
		
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
				bwOos = TRUE,
				isWithin = ifelse(in.sample.Bool,1,0))
		
		row.names(X.embed.2.oos)<-NULL
		row.names(X.embed.2.in)<-NULL

		X.embed.2.norm<- rbind(X.embed.2.in[1:4,] ,X.embed.2.oos,X.embed.2.in[5:(n-2),])		
		
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
		if (sum((X.embed.2.norm[6,2]>0.5))==1){ 
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
			print("First Test Point")
			print(X.embed.2.norm[5,])		
			if (stress<min.stress.2){
				min.stress.2<-stress
				print("Min stress found(second min)")
				print(X.embed.2.norm)
				print(min.stress.2)
				min.config.2<-X.embed.2.norm
			
			}
			far.to.init.1<-rbind(far.to.init.1,X.embed.2.norm[5,])
			far.to.init.2<-rbind(far.to.init.2,X.embed.2.norm[6,])

		}
		stress.at.loc[i,j]<- stress

		

	}

}



#plot(x=grid.seq.x,y=grid.seq.y, col=grid.resp)
#x.coords = [x1 x2 x3 x4 ... x1 x2 x3 x4 ...
x.coords <- rep(grid.seq.x,length(grid.seq.y))
#y.coords = [y1 y1 y1 y1 ... y2 y2 y2 y2 ...
y.coords <- rep(grid.seq.y,each=length(grid.seq.x))
#grid.resp<-grid.resp[length(grid.seq.x)A:1,]

plot(x.coords, y.coords,
#plot(unmatrix(mesh.grid.coords$x,byrow=FALSE),
#     unmatrix(mesh.grid.coords$y,byrow=FALSE),
      col=ifelse(unmatrix(grid.resp,byrow=FALSE)==1,"red","black"))
print("For w.value")
print(w)
print("Min stress found(real min)")
print(min.config.1)
print(min.stress.1)
print("Min stress found(second min)")
print(min.config.2)
print(min.stress.2)
min.config.stress.1.w<-min.stress.1
min.config.stress.2.w<-min.stress.2
		
		
if(!is.vector(close.to.init.1)){
windows()

par(pch=1)
plot(x=close.to.init.1[,1],y=close.to.init.1[,2],col="red",
			xlim=c(min(close.to.init.1[,1],close.to.init.2[,1]),max(close.to.init.1[,1],close.to.init.2[,1]))
			,ylim=c(min(close.to.init.1[,2],close.to.init.2[,2]),max(close.to.init.1[,2],close.to.init.2[,2])))

par(pch=3)
points(x=close.to.init.2[,1],y=close.to.init.2[,2],col="blue")
title(paste("Final config Close to true config- w=",w,collapse=TRUE))
legend("bottomright",legend=c(expression("X_5"),expression("X_6")),col=c("red","blue"))

savePlot(paste("true-min-w",w,".pdf",collapse=TRUE))

#select.x<- sort( sample.int(length(grid.seq.x) , 10))
select.x<-1:length(grid.seq.x) 
#select.y <- sort( sample.int(length(grid.seq.y) , 10))
select.y<-1:length(grid.seq.y) 
#The indexing (select.x,select.y) is mixed because mesh.grid function generates a matrix 
#whose columns are for x coordinates, while for final.coords, rows are for x coordinates
#arrows(x0 = unmatrix(mesh.grid.coords$x[select.y,select.x],byrow=FALSE),
#	 y0 = unmatrix(mesh.grid.coords$y[select.y,select.x],byrow=FALSE),
#	 x1 = unmatrix( final.coords.x[select.x,select.y],byrow=FALSE),
#       y1 = unmatrix( final.coords.y[select.x,select.y],byrow=FALSE) ,
#	length=0.1)


#arrows(x0 = mesh.grid.coords$x[t(grid.resp==1)],
#	 y0 = mesh.grid.coords$y[t(grid.resp==1)],
#	 x1 =  final.coords.x[grid.resp==1] ,
#       y1 =  final.coords.y[grid.resp==1] ,
#	length=0.1)

}

if(!is.vector(far.to.init.1)){
windows()

par(pch=1)
plot(x=far.to.init.1[,1],y=far.to.init.1[,2],col="red",
			xlim=c(min(grid.seq.x),max(grid.seq.x))
			,ylim=c(min(grid.seq.y),max(grid.seq.y)))
par(pch=3)
points(x=far.to.init.2[,1],y=far.to.init.2[,2],col="blue")
title(paste("Final config Far to true config- w=",w,collapse=TRUE))
legend("bottomright",legend=c(expression("X_5"),expression("X_6")),col=c("red","blue"))

savePlot(paste("other-min-w",w,".pdf",collapse=TRUE))
#arrows(x0 = mesh.grid.coords$x[t(grid.resp==0)],
#	 y0 = mesh.grid.coords$y[t(grid.resp==0)],
#	 x1 =  final.coords.x.6[grid.resp==0] ,
#       y1 =  final.coords.y.6[grid.resp==0] ,
#	length=0.1)
}

plot3d(x.coords, y.coords,


#plot3d(unmatrix(mesh.grid.coords$x,byrow=FALSE),
#	 unmatrix(mesh.grid.coords$y,byrow=FALSE),
		stress.at.loc)
}





