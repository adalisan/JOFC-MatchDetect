
meshgrid <- function(a,b) {
  list(
       x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+")
       )
} 

n<-7
X<-matrix(0,n,2)
X[,1]<-c(0,1,1,0,0,0,0.5)
X[,2]<-c(0,0,1,1,0,0,0.5)
X[5,]<-c(1.01,1)
X[6,]<-c(0.01,1)

in.sample.ind<-c(1:4,7)
d.X<- dist(X)
W<-matrix(0,n,n)

w<-0.2009
W<-matrix(1-w,7,7)
W[5,3] <- W[3,5]<- w

diag(W)<-0



i<-3
j<-6

ind<-(n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<-d.X[ind]-0.9

i<-4
j<-5
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<- d.X[ind]-0.9
#X[5,]<-c(0,1)


i<-5
j<-7
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<- d.X[ind]+0.9

init.config<-X

X.embed.1.in <- smacofM(as.matrix(d.X)[in.sample.ind,in.sample.ind],
                    ndim    = 2,
                    W       = W[in.sample.ind,in.sample.ind],
                    init    = init.config[in.sample.ind,]		,
                    verbose = TRUE,
                    itmax   = 1000,
                    eps     = 1e-6)
			
X.embed.1.oos <- oosIM(D=as.matrix(d.X),
				X=X.embed.1.in,
				init=init.config,
				W=W,
					verbose = TRUE,
					itmax   = 1000,
					eps     = 1e-6,
					bwOos = FALSE,
					isWithin = (1:7 %in% in.sample.ind) )
			
			
X.embed.1.norm<- rbind(X.embed.1.in[1:4,] ,X.embed.1.oos,X.embed.1.in[5,]) 
print(X.embed.1.norm)


close.to.init.1<-c(0.95,1)
close.to.init.2<-c(0.07,1)
far.to.init.1<-c(0.07,1)
far.to.init.2<-c(0.95,1)
min.stress.1<-100
min.stress.2<-100
grid.seq.x<-seq(-0.2,1.2,0.025)
grid.seq.y<-seq(0.4,1.4,0.05)

stress.at.loc<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))

grid.resp<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
grid.coords<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
mesh.grid.coords<-meshgrid(grid.seq.x,grid.seq.y)
for (j in 1:length(grid.seq.y)){

	for (i in 1:length(grid.seq.x)){

		i.x<-grid.seq.x[i]
		j.y<-grid.seq.y[j]
		#i.x<-mesh.grid.coords$x[i,j] 
		#j.y<-mesh.grid.coords$y[i,j] 

		init.config[5,]<-c(i.x,j.y)
		#init.config[6,]<-c(1-i.x,j.y)
		
		X.embed.2.in<-smacofM(as.matrix(d.X)[in.sample.ind, in.sample.ind],
				ndim    = 2,
				W       = W[in.sample.ind,in.sample.ind],
				init    = init.config[in.sample.ind,]	,
				verbose = FALSE,
				itmax   = 1000,
				eps     = 1e-6)
		X.embed.2.oos<-oosIM(D = as.matrix(d.X),
				X = X.embed.2.in,
				init=init.config,
				W = W,
				verbose = FALSE,
				itmax   = 1000,
				eps     = 1e-6,
				bwOos = FALSE,
				isWithin = (1:7 %in% in.sample.ind))
		
		X.embed.2.norm<- rbind(X.embed.2.in[1:4,] ,X.embed.2.oos,X.embed.2.in[5,])		

		print("X.embed.2.norm")
		#print(X.embed.2.norm)
		stress<- sum(as.dist(W)*((dist(X.embed.2.norm)-d.X)^2))
		print(stress)
		if ((sum(X.embed.2.norm[5,]>0.5)==2)){ 
			#& (sum(X.embed.2.norm[6,]<0.3)==1))
			

			grid.resp[i,j]<-1
			print("First Test Point")
			print(X.embed.2.norm[5,])
			if (stress<min.stress.1) {
				min.stress.1<-stress
				print("Min stress found(real min)")
				print(X.embed.2.norm)
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
			
			}
			far.to.init.1<-rbind(far.to.init.1,X.embed.2.norm[5,])
			far.to.init.2<-rbind(far.to.init.2,X.embed.2.norm[6,])

		}
		stress.at.loc[i,j]<- sum(as.dist(W)*((dist(init.config)-d.X)^2))

		

	}

}


#plot(x=grid.seq.x,y=grid.seq.y, col=grid.resp)

x.coords <- rep(grid.seq.x,length(grid.seq.y))
y.coords <- rep(grid.seq.y,each=length(grid.seq.x))
#grid.resp<-grid.resp[length(grid.seq.x):1,]

plot(x.coords, y.coords,

#plot(unmatrix(mesh.grid.coords$x,byrow=FALSE),
#     unmatrix(mesh.grid.coords$y,byrow=FALSE),
      col=ifelse(unmatrix(grid.resp,byrow=FALSE)==1,"red","black"))
print("Min stress found(real min)")

print(min.stress.1)
print("Min stress found(second min)")

print(min.stress.2)

windows()

par(pch=1)
plot(x=close.to.init.1[,1],y=close.to.init.1[,2],col="red",xlim=c(0,1),ylim=c(0,1.5))
par(pch=3)
points(x=close.to.init.2[,1],y=close.to.init.2[,2],col="blue")
title("Final config Close to true points")

windows()

par(pch=1)
plot(x=far.to.init.1[,1],y=far.to.init.1[,2],col="red",xlim=c(0,1),ylim=c(0,1.5))
par(pch=3)
points(x=far.to.init.2[,1],y=far.to.init.2[,2],col="blue")
title("Final config  Far to true points")


windows()


plot3d(x.coords, y.coords,


#plot3d(unmatrix(mesh.grid.coords$x,byrow=FALSE),
#	 unmatrix(mesh.grid.coords$y,byrow=FALSE),
		stress.at.loc)






