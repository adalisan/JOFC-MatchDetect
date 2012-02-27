


n<-6
X<-matrix(0,n,2)
X[,1]<-c(0,1,1,0,0,0)
X[,2]<-c(0,0,1,1,0,0)
X[5,]<-c(0.95,1)
X[6,]<-c(0.07,1)
d.X<- dist(X)
W<-matrix(0,n,n)

w<-0.2009

W[1,]<-c(0  ,1-w, 1-w, 1-w,1-w,1-w)
W[2,]<-c(1-w,0  ,1-w ,  1-w,1-w,1-w)
W[3,]<-c(1-w,1-w  , 0,  1-w,w,1-w)
W[4,]<-c(1-w  ,1-w  , 1-w,0,1-w,1-w)
W[5,]<-c(1-w  ,1-w  , w,1-w, 0,1-w)
W[6,]<-c(1-w  ,1-w  , 1-w,1-w,1-w,0)

i<-3
j<-6

ind<-(n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<-d.X[ind]-0.92
i<-4
j<-5
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<- d.X[ind]-0.92
init.config<-X
#X[5,]<-c(0,1)
X.embed.1<-smacofM(d.X,
                    ndim    = 2,
                    W       = W,
                    init    = init.config		,
                    verbose = TRUE,
                    itmax   = 1000,
                    eps     = 1e-6)

p.1<-procrustes(X.embed.1[1:4,],X[1:4,],translation=TRUE,dilation=FALSE)

X.embed.1.norm<- X.embed.1 %*%p.1$R + matrix(rep(p.1$tt,n),nrow=n,ncol=2,byrow=TRUE)
print(X.embed.1.norm)

min.stress.1<-100
min.stress.2<-100
grid.seq.x<-seq(-0.2,1.2,0.05)
grid.seq.y<-seq(0.5,1.4,0.05)

grid.resp<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
grid.coords<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
	for (j in 1:length(grid.seq.y)){

for (i in 1:length(grid.seq.x)){

		i.x<-grid.seq.x[i]
		j.y<-grid.seq.y[j] 
		init.config[5,]<-c(i.x,j.y)
		#init.config[6,]<-c(1-i.x,j.y)
		X.embed.2<-smacofM(d.X,
                    ndim    = 2,
                    W       = W,
                    init    = init.config,
                    verbose = FALSE,
                    itmax   = 1000,
                    eps     = 1e-6)

p.2<-procrustes(X.embed.2[1:4,],X[1:4,],translation=TRUE,dilation=FALSE)

X.embed.2.norm<- X.embed.2 %*%p.1$R + matrix(rep(p.2$tt,n),nrow=n,ncol=2,byrow=TRUE)
print("X.embed.2.norm")
#print(X.embed.2.norm)
stress<- sum(as.dist(W)*((dist(X.embed.2.norm)-d.X)^2))
		print(stress)
		if ((sum(X.embed.2.norm[5,]>0.5)==2)) #& (sum(X.embed.2.norm[6,]<0.3)==1))
		{
		grid.resp[i,j]<-1
		print(X.embed.2.norm[5,])
		if (stress<min.stress.1){
		min.stress.1<-stress
		print("Min stress found(real min)")
		print(X.embed.2.norm)
		print(min.stress.1)

		}

		}else{
		grid.resp[i,j]<-0
		print(X.embed.2.norm[5,])		
if (stress<min.stress.2){
		min.stress.2<-stress
		print("Min stress found(second min)")
		print(X.embed.2.norm)
		print(min.stress.2)

		}

		}

	}

}


#plot(x=grid.seq.x,y=grid.seq.y, col=grid.resp)

x.coords <- rep(grid.seq.x,length(grid.seq.y))
y.coords <- rep(grid.seq.y,each=length(grid.seq.x))
#grid.resp<-grid.resp[length(grid.seq.x):1,]
plot(x.coords,y.coords,col=ifelse(unmatrix(grid.resp,byrow=FALSE)==1,"red","black"))
print(min.stress.1)
print(min.stress.2)
