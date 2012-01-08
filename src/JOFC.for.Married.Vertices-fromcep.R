

require(igraph)
require(optmatch)
source("./lib/simulation_math_util_fn.R")
source("./lib/smacofM.R")
source("./lib/oosIM.R")
source("./lib/oosMDS.R")
source("./lib/diffusion_distance.R")
source("./lib/graph_embedding_fn.R")

cep=TRUE
verbose= FALSE
oos=TRUE



n = 20
m = 15 # the first m pairs are known matches ; the last n-m pairs are to-be-matched

npert = 11
nmc = 100
pert=(0:10)/10
nc.jofc.diff = matrix(0,npert,nmc)
nc.jofc.weighted = matrix(0,npert,nmc)
nc.jofc.unweighted = matrix(0,npert,nmc)



nc.cmds = matrix(0,npert,nmc)

matched.cost<-0.01

w.vals<-0.8

matched.cost<-0.01 #If matched.cost is equal to 1, consider an unweighted graph, with edges between matched vertices
#If matched.cost is  between 0 and 1, the graph is weighted with edges between matched vertices with weights equal to matched.cost. Edges between 
# vertices of the same condition have either weight 1 or 2 according to whether they're connected according to given adjacency matrix or not.

d.dim<-4
T.diff<-2


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
    Y.emb<-NULL
		Gp<-bitflip(G ,pert[ipert],pert[ipert])
		Graph.1<-graph.adjacency(G)
		Graph.2<-graph.adjacency(Gp)
		D.1<-shortest.paths(Graph.1)
		D.2<-shortest.paths(Graph.2)
		oos.sampling<-sample(1:n, size=(n-m), replace=FALSE)
		in.sample.ind<-rep(TRUE,2*n)
		in.sample.ind[oos.sampling]<-FALSE
		in.sample.ind[n+oos.sampling]<-FALSE
    #if (imc==1) print(in.sample.ind)
    
    J.1 =jofc.diffusion.dist(G,Gp,
  	in.sample.ind,
		d.dim=d.dim,
		graph.is.directed=FALSE,
    oos=oos,
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE)
		
		M = solveMarriage(J.1)
		nc.jofc.diff[ipert,imc] = present(M)         # returns the number correct in the marriage?
    
     
    J.2 =jofc(G,Gp,
    in.sample.ind,
		d.dim=d.dim,
		graph.is.directed=FALSE,	
     oos=oos,          
              
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE)
		
		M = solveMarriage(J.2)
		nc.jofc.weighted[ipert,imc] = present(M)         # returns the number correct in the marriage?
    
     J.3 =jofc(G,Gp,
    in.sample.ind,
  	d.dim=d.dim,
		graph.is.directed=FALSE,
    oos=oos,
    use.weighted.graph=FALSE,
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE)
		
		M = solveMarriage(J.3)
		nc.jofc.unweighted[ipert,imc] = present(M)         # returns the number correct in the marriage?
    
    
    
    
    

			myD.M = D.M # i use just the object D.M ... none of the original entries!
			myD.M[1:n,1:n]=D.1
			myD.M[(n+1):(2*n),(n+1):(2*n)]=D.2
			myD.M[1:n,(n+1):(2*n)]=(D.1+D.2)/2
			myD.M[(n+1):(2*n),1:n]=(D.1+D.2)/2
			for(i in 1:n) for(j in (n+m+1):(2*n)) myD.M[i,j] = myD.M[j,i] = 0
			for(i in (m+1):n) for(j in (n+1):(2*n)) myD.M[i,j] = myD.M[j,i] = 0
			
      
       if (oos){
         if (imc==1) print(in.sample.ind)
       myD.M.in<- myD.M[in.sample.ind,in.sample.ind]
  			ccc = cmdscale(myD.M.in,k=d.dim,eig=T)
 			#plot(ccc$eig)
 			#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
 			#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
       Y.emb<-oosMDS(myD.M,X=ccc$points, w=ifelse(in.sample.ind,1,0),init="gower")
       
# 			
 			 U = as.matrix(dist(Y.emb[,c(2,3)]))
       } else {
      
      
  		ccc = cmdscale(myD.M,k=d.dim,eig=T)
			#plot(ccc$eig)
			#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
			#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
			
			U = as.matrix(dist(ccc$points[ c((m+1):(n),(n+m+1):(n+n)) ,c(2,3)]))
 		  }

      
      
      
      
      
			M = fullmatch(U[1:(n-m),(n-m+1):(2*(n-m))])
			nc.cmds[ipert,imc] = present(M)         # returns the number correct in the marriage?
		
		
		
		
	}
}


### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot1.pdf")
colors.vec<-c( "red","blue","orange","green")
plot(pert,apply(nc.jofc.diff,1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[1])
points(pert,apply(nc.cmds,1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=2,col=colors.vec[2])

points(pert,apply(nc.jofc.weighted,1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=3,col=colors.vec[3])

points(pert,apply(nc.jofc.unweighted,1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=4,col=colors.vec[4])

legend.txt<- c("diff dist","CMDS","weighted.graph","unweighted.graph")
legend(x="topright",legend=legend.txt, col =colors.vec,pch=1:4)
title("marriage o jofc")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()




d.dim.vec<-c(4,6,10)

nc.jofc.diff.d    = array(0,dim=c(npert,nmc,length(d.dim.vec)))
nc.jofc.weighted.d = array(0,dim=c(npert,nmc,length(d.dim.vec)))
nc.jofc.unweighted.d = array(0,dim=c(npert,nmc,length(d.dim.vec)))



nc.cmds.d = array(0,dim=c(npert,nmc,length(d.dim.vec)))

# these three lines are not used, really ...
# they just set up a dummy D.M which is the correct object for input to fullmatch".
# i'm just using Sancar's object ... i don't really understand it!
G<-ER(n,0.5)
G.comb<-omnibusM(G,G,diag(n))
Graph.M <- graph.adjacency(G.comb,weighted= NULL ,mode="undirected")
D.M<-shortest.paths(Graph.M)

seed<-123

for (d.dim in d.dim.vec){
  d.i<-which(d.dim.vec==d.dim)
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
    
    
    J.1 =jofc.diffusion.dist(G,Gp,
  	in.sample.ind,
		d.dim=d.dim,
		graph.is.directed=FALSE,
    oos=oos,
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE)
		
		M = solveMarriage(J.1)
		nc.jofc.diff.d[ipert,imc,d.i] = present(M)         # returns the number correct in the marriage?
    
     
    J.2 =jofc(G,Gp,
    in.sample.ind,
		d.dim=d.dim,
		graph.is.directed=FALSE,
     oos=oos,          
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE)
		
		M = solveMarriage(J.2)
		nc.jofc.weighted.d[ipert,imc,d.i] = present(M)         # returns the number correct in the marriage?
    
     J.3 =jofc(G,Gp,
    in.sample.ind,
  	d.dim=d.dim,
		graph.is.directed=FALSE,
     oos=oos,
               use.weighted.graph=FALSE,
		wt.matrix.1=NULL,
		wt.matrix.2=NULL,
		sep.graphs=TRUE)
		
		M = solveMarriage(J.3)
		nc.jofc.unweighted.d[ipert,imc,d.i] = present(M)         # returns the number correct in the marriage?
    
    
    
    
    
		if (cep){
      U<-matrix()
      
			myD.M = D.M # i use just the object D.M ... none of the original entries!
			myD.M[1:n,1:n]=D.1
			myD.M[(n+1):(2*n),(n+1):(2*n)]=D.2
			myD.M[1:n,(n+1):(2*n)]=(D.1+D.2)/2
			myD.M[(n+1):(2*n),1:n]=(D.1+D.2)/2
			for(i in 1:n) for(j in (n+m+1):(2*n)) myD.M[i,j] = myD.M[j,i] = 0
			for(i in (m+1):n) for(j in (n+1):(2*n)) myD.M[i,j] = myD.M[j,i] = 0
			
      
      if (oos){
         if (imc==1) print(in.sample.ind)
       myD.M.in<- myD.M[in.sample.ind,in.sample.ind]
    		ccc = cmdscale(myD.M.in,k=d.dim,eig=T)
 			#plot(ccc$eig)
 			#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
 			#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
       Y.emb<-oosMDS(myD.M,X=ccc$points, w=ifelse(in.sample.ind,1,0),init="gower")
       
# 			
 			 U = as.matrix(dist(Y.emb[,c(2,3)]))
      } else {
      
      
			ccc = cmdscale(myD.M,k=d.dim,eig=T)
			#plot(ccc$eig)
			#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
			#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
			
			U = as.matrix(dist(ccc$points[ c((m+1):(n),(n+m+1):(n+n)) ,c(2,3)]))
		  }
			M = fullmatch(U[1:(n-m),(n-m+1):(2*(n-m))])
			nc.cmds.d[ipert,imc,d.i] = present(M)         # returns the number correct in the marriage?
		}
		
		
		
	}
}
}

### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot2.pdf")
colors.vec<-c( "red","blue","orange","green")
plot(pert,apply(nc.jofc.diff.d[,,1],1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[1])
for (d.i in 2:length(d.dim.vec)){
points(pert,apply(nc.jofc.diff.d[,,d.i],1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=1,col=colors.vec[d.i])
}
for (d.i in 1:length(d.dim.vec)){
points(pert,apply(nc.cmds.d[,,d.i],1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=2,col=colors.vec[d.i])
}
for (d.i in 1:length(d.dim.vec)){
points(pert,apply(nc.jofc.weighted.d[,,d.i],1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=3,col=colors.vec[d.i])
}
for (d.i in 1:length(d.dim.vec)){
points(pert,apply(nc.jofc.unweighted.d[,,d.i],1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=4,col=colors.vec[d.i])
}
legend.txt<- c("diff dist","CMDS","weighted.graph","unweighted.graph")
legend.txt <- c(legend.txt,d.dim.vec)
legend(x="topright",legend=legend.txt, col =c(rep(1,4),colors.vec),pch=c(1:4,rep(1,4)))
title("marriage o jofc")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()







