# TODO: Add comment
# Testing function for Zhiliang's weighted oos
# Author: Sancar
###############################################################################
library(MASS)
library(MCMCpack)
source("./src/smacofM.R")
source("./src/oosIM.R")
source("./src/simulation_math_util_fn.R")
source("./src/oosMDS.R")
X <- mvrnorm(100,mu=rep(0,2), Sigma=diag(2))
dist.X <- as.matrix(dist(X))
X.embed.1<-smacofM(D=as.matrix(dist.X),ndim=2,init=cmdscale(D,k=2))

in.sample <- sample(1:100,4)
in.sample <- sort(in.sample)
D.in <- dist.X[in.sample,in.sample]
X.embed.in.2<-smacofM(D=D.in,ndim=2,init=X.embed.1[in.sample,])
in.sample.ind<-rep(0,100)
in.sample.ind[in.sample]<-1

#X.embed.2.oos<-smacofOos(dist.X,X.embed.in.2,oos.flag=in.sample.ind,init=X.embed.1[-in.sample,])
X.embed.2.oos   <-oosIM(dist.X,X.embed.in.2,isWithin=in.sample.ind,init=X.embed.1[-in.sample,])
X.embed.2<-array(0,dim=c(100,2))
X.embed.2[in.sample,]<-X.embed.in.2
X.embed.2[-in.sample,]<-X.embed.2.oos

proc<-procrustes(X.embed.2,X.embed.1)
config.dist <- norm(X.embed.1-proc$X.new,'F')
config.dist

d2<-dist(X.embed.2)
 norm(dist.X-as.matrix(d2),'F')/10000



X.cmds.1<-cmdscale(dist.X,k=2)
X.cmds.in.2<- cmdscale(D.in,k=2)
X.cmds.oos<- oosMDS(dist.X,X.cmds.in.2,w=in.sample.ind)

X.cmds.embed.2<-array(0,dim=c(100,2))
X.cmds.embed.2[in.sample,]<-X.cmds.in.2
X.cmds.embed.2[-in.sample,]<-X.cmds.oos
proc.2<-procrustes(X.cmds.embed.2,X.cmds.1)
config.dist.2 <- norm(X.cmds.1-proc.2$X.new,'F')
config.dist.2

d2<-dist(X.cmds.embed.2)
 norm(dist.X-as.matrix(d2),'F')/10000







plot(X.embed.1,col="red")
pr.1<-procrustes(X.embed.in.2,X.embed.1[in.sample,])


Weight.Mat<-w.val.to.W.mat(0.2,100,FALSE)
w.embed.1 <- smacofM(dist.X,2,W=Weight.Mat)

w.embed.2.in<- smacofM(D.in,2,W=Weight.Mat)
w.embed.2.oos<-oosIM(dist.X,w.embed.2.in,W= Weight.Mat,isWithin=in.sample.ind)


w.embed.2<-array(0,dim=c(100,2))
w.embed.2[in.sample,]<-w.embed.2.in
w.embed.2[-in.sample,]<-w.embed.2.oos

config.dist.2<-norm(w.embed.1-w.embed.2,'F')
