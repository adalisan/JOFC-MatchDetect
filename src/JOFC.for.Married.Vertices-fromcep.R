

require(igraph)
require(optmatch)
source("./lib/simulation_math_util_fn.R")
source("./lib/smacofM.R")
source("./lib/oosIM.R")
source("./lib/diffusion_distance.R")
source("./lib/graph_embedding_fn.R")

cep=TRUE
verbose= FALSE

n = 20
m = 15 # the first m pairs are known matches ; the last n-m pairs are to-be-matched

npert = 11
nmc = 100
pert=(0:10)/10
nc.jofc = matrix(0,npert,nmc)
nc.cmds = matrix(0,npert,nmc)
matched.cost<-0.01

w.vals<-0.8

matched.cost<-0.01 #If matched.cost is equal to 1, consider an unweighted graph, with edges between matched vertices
#If matched.cost is  between 0 and 1, the graph is weighted with edges between matched vertices with weights equal to matched.cost. Edges between 
# vertices of the same condition have either weight 1 or 2 according to whether they're connected according to given adjacency matrix or not.
debug.mode<-FALSE
if (debug.mode){
	n.v<-100
	n.np<- 5
	m<- 5
}
d.dim<-4
T.diff<-100


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
		oos.sampling<-sample.int(n,size=n-m,replace=FALSE)
		in.sample.ind<-rep(TRUE,2*n)
		in.sample.ind[oos.sampling]<-FALSE
		in.sample.ind[n+oos.sampling]<-FALSE
    
    
    J =jofc.diffusion.dist(G,Gp,
  	in.sample.ind,
		d.dim=d.dim,
		graph.is.directed=FALSE,	
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE)
		
		M = solveMarriage(J)
		nc.jofc[ipert,imc] = present(M)         # returns the number correct in the marriage?
		if (cep){
			myD.M = D.M # i use just the object D.M ... none of the original entries!
			myD.M[1:n,1:n]=D.1
			myD.M[(n+1):(2*n),(n+1):(2*n)]=D.2
			myD.M[1:n,(n+1):(2*n)]=(D.1+D.2)/2
			myD.M[(n+1):(2*n),1:n]=(D.1+D.2)/2
			for(i in 1:n) for(j in (n+m+1):(2*n)) myD.M[i,j] = myD.M[j,i] = 0
			for(i in (m+1):n) for(j in (n+1):(2*n)) myD.M[i,j] = myD.M[j,i] = 0
			
			ccc = cmdscale(myD.M,k=4,eig=T)
			#plot(ccc$eig)
			#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
			#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
			
			U = as.matrix(dist(ccc$points[ c((m+1):(n),(n+m+1):(n+n)) ,c(2,3)]))
			M = fullmatch(U[1:(n-m),(n-m+1):(2*(n-m))])
			nc.cmds[ipert,imc] = present(M)         # returns the number correct in the marriage?
		}
		
		
		
	}
}

apply(nc.jofc,1,mean)
### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot1.pdf")
plot(pert,apply(nc.jofc,1,mean)/(n-m),xlab="perturbation parameter",ylab="% correct matches",ylim=c(0,1),col="red")
points(pert,apply(nc.cmds,1,mean)/(n-m),xlab="perturbation parameter",ylab="% correct matches",ylim=c(0,1),pch=2,col="blue")
#legend("bottomright")
title("marriage o jofc")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()
