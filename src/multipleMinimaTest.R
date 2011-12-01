X<-matrix(0,4,2)
X[,1]<-c(0,1,1,0)
X[,2]<-c(0,0,1,1)
d.X<- dist(X)
W<-matrix(0,4,4)

w<-0.5
W[1,]<-c(0,1-w, w,0)
W[2,]<-c(1-w,0, 0,w)
W[3,]<-c(w,0, 0,1-w)
W[4,]<-c(0,w, 1-w,0)

init.config<-X
X.embed<-smacofM(d.X,
                    ndim    = 2,
                    W       = W,
                    init    = init.config		,
                    verbose = FALSE,
                    itmax   = 1000,
                    eps     = 1e-6)

init.config[3,]<-c(0,1)
init.config[4,]<-c(1,1)
X.embed<-smacofM(d.X,
                    ndim    = 2,
                    W       = W,
                    init    = init.config,
                    verbose = FALSE,
                    itmax   = 1000,
                    eps     = 1e-6)