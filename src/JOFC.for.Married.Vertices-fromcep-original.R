require(igraph)
require(optmatch)
source("simulation_math_util_fn.R")
source("smacofM.R")
source("oosIM.R")
source("graph_embedding_fn.R")

 ER = function(n,p)
 {
 A = matrix( rbinom(n^2,1,p) , ncol=n,byrow=T )
 A=A*upper.tri(A)
 A=A+t(A)
 return(A)
 }

 bitflip = function(G,q10,q01,binary=T,symmetric=T,hollow=T)
 # takes graph (adjacency matrix) G and flips 1s to 0s with prob q10 and flips 0s to 1s with prob q01
 # assumes binary=T,symmetric=T,hollow=T
 {
 n=dim(G)[1]
 for(u in 1:(n-1))
 for(v in (u+1):n)
 G[u,v]=G[v,u]=ifelse(G[u,v],1-rbinom(1,1,q10),rbinom(1,1,q01))
 return(G)
 }

n = 20
m = 15 # the first m pairs are known matches ; the last n-m pairs are to-be-matched
Ln=LETTERS[1:n]
colvec=c(rep(1,m),rep(2,(n-m)),rep(3,m),rep(4,(n-m)))
d<-4
npert = 11
nmc = 1000
pert=(0:10)/10
nc = matrix(0,npert,nmc)

 # these three lines are not used, really ...
 # they just set up a dummy D.M which is the correct object for input to fullmatch".
 # i'm just using Sancar's object ... i don't really understand it!
 G<-ER(n,0.5)
 G.comb<-omnibusM(G,G,diag(n))
 Graph.M <- graph.adjacency(G.comb,weighted= NULL ,mode="undirected")
 D.M<-shortest.paths(Graph.M)

seed<-123
for(imc in 1:nmc)
 {
 G<-ER(n,0.5)
 for(ipert in 1:npert)
 {
 Gp<-bitflip(G ,pert[ipert],pert[ipert])
 Graph.1<-graph.adjacency(G)
 Graph.2<-graph.adjacency(Gp)
 D.1<-shortest.paths(Graph.1)
 D.2<-shortest.paths(Graph.2)

 myD.M = D.M # i use just the object D.M ... none of the original entries!
 myD.M[1:n,1:n]=D.1
 myD.M[(n+1):(2*n),(n+1):(2*n)]=D.2
 myD.M[1:n,(n+1):(2*n)]=(D.1+D.2)/2
 myD.M[(n+1):(2*n),1:n]=(D.1+D.2)/2
 for(i in 1:n) for(j in (n+m+1):(2*n)) myD.M[i,j] = myD.M[j,i] = 0

 ccc = cmdscale(myD.M,k=4,eig=T)
 #plot(ccc$eig)
 #pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
 #plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))

 U = as.matrix(dist(ccc$points[ c((m+1):(n),(n+m+1):(n+n)) ,c(2,3)]))
 M = fullmatch(U[1:(n-m),(n-m+1):(2*(n-m))])
 nc[ipert,imc] = present(M)         # returns the number correct in the marriage?
 }
 }

apply(nc,1,mean)
### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot1.pdf")
plot(pert,apply(nc,1,mean)/(n-m),xlab="perturbation parameter",ylab="% correct matches",ylim=c(0,1))
title("marriage o jofc")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()