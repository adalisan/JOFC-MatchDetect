
X<-matrix(0,6,2)
X[,1]<-c(0,1,1,0,0,0)
X[,2]<-c(0,0,1,1,0,0)
X[5,]<-c(1,1)
X[6,]<-c(0,1)
d.X<- dist(X)
W<-matrix(0,6,6)

w<-0.2
W[1,]<-c(0  ,1-w, 1-w, 1-w,1-w,1-w)
W[2,]<-c(1-w,0  ,1-w ,  1-w,1-w,1-w)
W[3,]<-c(1-w,1-w  , 0,  1-w,w,1-w)
W[4,]<-c(1-w  ,1-w  , 1-w,0,1-w,w)
W[5,]<-c(1-w  ,1-w  , w,1-w, 0, 1-w)
W[6,]<-c(1-w  ,1-w  , 1-w,w,1-w,0)
init.config<-X
X.embed.1<-smacofM(d.X,
                    ndim    = 2,
                    W       = W,
                    init    = init.config		,
                    verbose = TRUE,
                    itmax   = 1000,
                    eps     = 1e-6)

grid.seq.x<-seq(0,1,0.05)
grid.seq.y<-seq(0.6,1.3,0.05)

grid.resp<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y))
for (i in 1:length(grid.seq.x))

	for (j in length(grid.seq.y)){
i.x<-grid.seq.x[i]
j.y<-grid.seq.y[j] 
init.config[5,]<-c(i.x,j.y)
init.config[6,]<-c(1-i.x,j.y)
X.embed.2<-smacofM(d.X,
                    ndim    = 2,
                    W       = W,
                    init    = init.config,
                    verbose = TRUE,
                    itmax   = 1000,
                    eps     = 1e-6)
X.embed.2.norm<-(X.embed.2+cbind(rep(abs(min(X.embed.2[,1])),6),rep(abs(min(X.embed.2[,2])),6)))
print(X.embed.2.norm)
if ((sum(X.embed.2.norm[5,]>0.92)==2)& (sum(X.embed.2.norm[6,]>0.92)==1))
{
grid.resp[i,j]<-"blue"
}else{
grid.resp[i,j]<-"red
}
}
