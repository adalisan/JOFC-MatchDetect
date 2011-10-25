# TODO: Code for investigating limit(n->\infty) of stress per point
# 
# Author: Sancar
###############################################################################
require(MASS)
source("smacofM.R")
n.vals<-c(20,50,100,200,500,1000,2000,5000,10000)
t<-15
avg.stress.avg<-c()
avg.strain.avg<-c()
for (n in n.vals){

	avg.stress.n<-0
d.prime <- 2
d<-1
for (s in 1:t){
X <-mvrnorm(n,mu=rep(0,d.prime),Sigma=diag(d.prime))
D <- dist(X)
init.conf<-NULL
Y<- smacofM(
D,d    ,	W=matrix(1,nrow=n,ncol=n)        ,
init    = init.conf,
verbose = FALSE,
itmax   = 1000,
eps     = 1e-6) 
Y.cmds<-cmdscale(D,d)

stress<- sum((dist(Y)-D)^2)
strain<- sum((dist(Y.cmds)-D)^2)

avg.stress <-stress/(n*(n-1)/2)
avg.stress.n <- (avg.stress.n+avg.stress/t)
avg.strain <-strain/(n*(n-1)/2)
avg.strain.n <- (avg.strain.n+avg.strain/t)



}


avg.stress.avg<-c(avg.stress.avg,avg.stresss.n)
print(avg.stress.avg)

avg.strain.avg<-c(avg.strain.avg,avg.stra.n)

}
