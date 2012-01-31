

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
oos.cep = TRUE

a.cep <-20


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


w.vals.vec <- c(0.5,0.7,0.9,0.95)

w.max.index<-length(w.vals.vec)

matched.cost<-0.01 #If matched.cost is equal to 1, consider an unweighted graph, with edges between matched vertices
#If matched.cost is  between 0 and 1, the graph is weighted with edges between matched vertices with weights equal to matched.cost. Edges between 
# vertices of the same condition have either weight 1 or 2 according to whether they're connected according to given adjacency matrix or not.

d.dim<-4
T.diff<-1

dims.for.dist<-2:d.dim
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
				w.vals.vec,
				graph.is.directed=FALSE,
				oos=oos,
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE
		)
		
		
		M = solveMarriage(J.1[[1]])
		nc.jofc.diff[ipert,imc] = present(M)         # returns the number correct in the marriage?
		
		
		J.2 =jofc(G,Gp,
				in.sample.ind,
				d.dim=d.dim,
				w.vals.vec,
				graph.is.directed=FALSE,	
				oos=oos,          
				
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE)
		
		M = solveMarriage(J.2[[1]])
		nc.jofc.weighted[ipert,imc] = present(M)         # returns the number correct in the marriage?
		
		J.3 =jofc(G,Gp,
				in.sample.ind,
				d.dim=d.dim,
				w.vals.vec,	
				graph.is.directed=FALSE,
				oos=oos,
				use.weighted.graph=FALSE,
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE)
		
		M = solveMarriage(J.3[[1]])
		nc.jofc.unweighted[ipert,imc] = present(M)         # returns the number correct in the marriage?
		
		
		
		
		
		
		myD.M = D.M
		myD.M[1:n,1:n]=D.1
		myD.M[(n+1):(2*n),(n+1):(2*n)]=D.2
		myD.M[1:n,(n+1):(2*n)]=(D.1+D.2)/2
		myD.M[(n+1):(2*n),1:n]=(D.1+D.2)/2
		for(i in 1:n) for(j in (n+m+1):(2*n)) myD.M[i,j] = myD.M[j,i] = a.cep
		for(i in (m+1):n) for(j in (n+1):(2*n)) myD.M[i,j] = myD.M[j,i] = a.cep
		
		
		if (oos.cep){
			oos.in.sample.vec<- c(1:m,(n+1):(n+m))
			oos.in.sample.ind<- rep(FALSE,2*n)
			oos.in.sample.ind[oos.in.sample.vec]<-TRUE
			myD.M.in<- myD.M[oos.in.sample.ind, oos.in.sample.ind]
			ccc = cmdscale(myD.M.in,k=d.dim,eig=T)
			#plot(ccc$eig)
			#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
			#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
			Y.emb<-oosMDS(myD.M,X=ccc$points, w=ifelse(oos.in.sample.ind,1,0),init="gower")			 			
			U = as.matrix(dist(Y.emb[,dims.for.dist]))
		} else {

			ccc = cmdscale(myD.M,k=d.dim,eig=T)
			#plot(ccc$eig)
			#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
			#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
			
			U = as.matrix(dist(ccc$points[ c((m+1):(n),(n+m+1):(n+n)) ,dims.for.dist]))
		}
		
		
		
		
		
		
		M = fullmatch(U[1:(n-m),(n-m+1):(2*(n-m))])
		nc.cmds[ipert,imc] = present(M)         # returns the number correct in the marriage?
		
		
		
		
	}
}


### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot1.pdf")
colors.vec<-c( "red","blue","orange","green")
colors.vec<-colors.vec[1:4]
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

nc.jofc.diff.d       = array(0,dim=c(npert,nmc,length(d.dim.vec)))
nc.jofc.weighted.d   = array(0,dim=c(npert,nmc,length(d.dim.vec)))
nc.jofc.unweighted.d = array(0,dim=c(npert,nmc,length(d.dim.vec)))



nc.cmds.d           = array(0,dim=c(npert,nmc,length(d.dim.vec)))

# these three lines are not used, really ...
# they just set up a dummy D.M which is the correct object for input to fullmatch".
# i'm just using Sancar's object ... i don't really understand it!
G<-ER(n,0.5)
G.comb<-omnibusM(G,G,diag(n))
Graph.M <- graph.adjacency(G.comb,weighted= NULL ,mode="undirected")
D.M<-shortest.paths(Graph.M)

seed<-123

for (d.dim in d.dim.vec){
	dims.for.dist<-2:d.dim
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
			oos.sampling<-sample(1:n,size=n-m,replace=FALSE)
			in.sample.ind<-rep(TRUE,2*n)
			in.sample.ind[oos.sampling]<-FALSE
			in.sample.ind[n+oos.sampling]<-FALSE
			
			
			J.1 =jofc.diffusion.dist(G,Gp,
					in.sample.ind,
					d.dim=d.dim,
					w.vals.vec,
					graph.is.directed=FALSE,
					oos=oos,
					wt.matrix.1=NULL,
					wt.matrix.2=NULL,
					sep.graphs=TRUE)
			
			M = solveMarriage(J.1[[1]])
			nc.jofc.diff.d[ipert,imc,d.i] = present(M)         # returns the number correct in the marriage?
			
			
			J.2 =jofc(G,Gp,
					in.sample.ind,
					
					d.dim=d.dim,
					w.vals.vec,
					graph.is.directed=FALSE,
					oos=oos,          
					wt.matrix.1=NULL,
					wt.matrix.2=NULL,
					sep.graphs=TRUE)
			
			M = solveMarriage(J.2[[1]])
			nc.jofc.weighted.d[ipert,imc,d.i] = present(M)         # returns the number correct in the marriage?
			
			J.3 =jofc(G,Gp,
					in.sample.ind,
					d.dim=d.dim,
					w.vals.vec,
					graph.is.directed=FALSE,
					oos=oos,
					use.weighted.graph=FALSE,
					wt.matrix.1=NULL,
					wt.matrix.2=NULL,
					sep.graphs=TRUE)
			
			M = solveMarriage(J.3[[1]])
			nc.jofc.unweighted.d[ipert,imc,d.i] = present(M)         # returns the number correct in the marriage?
			
			
			
			
			
			if (cep){
				U<-matrix()
				
				myD.M = D.M # i use just the object D.M ... none of the original entries!
				myD.M[1:n,1:n]=D.1
				myD.M[(n+1):(2*n),(n+1):(2*n)]=D.2
				myD.M[1:n,(n+1):(2*n)]=(D.1+D.2)/2
				myD.M[(n+1):(2*n),1:n]=(D.1+D.2)/2
				for(i in 1:n) for(j in (n+m+1):(2*n)) myD.M[i,j] = myD.M[j,i] = a.cep
				for(i in (m+1):n) for(j in (n+1):(2*n)) myD.M[i,j] = myD.M[j,i] = a.cep
				
				
				if (oos.cep){
					
					oos.in.sample.vec<- c(1:m,(n+1):(n+m))
					oos.in.sample.ind<- rep(FALSE,2*n)
					oos.in.sample.ind[oos.in.sample.vec]<-TRUE
					myD.M.in<- myD.M[oos.in.sample.ind,oos.in.sample.ind]
					
					ccc = cmdscale(myD.M.in,k=d.dim,eig=T)
					#plot(ccc$eig)
					#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
					#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
					Y.emb<-oosMDS(myD.M,X=ccc$points, w=ifelse(oos.in.sample.ind,1,0),init="gower")					
# 			
					U = as.matrix(dist(Y.emb[,dims.for.dist]))
				} else {
					
					
					ccc = cmdscale(myD.M,k=d.dim,eig=T)
					#plot(ccc$eig)
					#pairs(ccc$points , col=colvec,pch=c(Ln,Ln))
					#plot(ccc$points[,c(2,3)] , col=colvec,pch=c(Ln,Ln))
					
					U = as.matrix(dist(ccc$points[ c((m+1):(n),(n+m+1):(n+n)) ,dims.for.dist]))
				}
				M = fullmatch(U[1:(n-m),(n-m+1):(2*(n-m))])
				nc.cmds.d[ipert,imc,d.i] = present(M)         # returns the number correct in the marriage?
			}
			
						
		}
	}
}

### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot2.pdf")
colors.vec<-c( "red","blue","orange","green","purple")

colors.vec<-colors.vec[1:length(d.dim.vec)]
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
legend(x="topright",legend=legend.txt, 
		col =c(rep(1,4),colors.vec),pch=c(1:4,rep(1,length(d.dim.vec))))
title("marriage o jofc")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()


pdf("plot3.pdf")
ccolors.vec<-c( "red","blue","orange","green","purple")

colors.vec<-colors.vec[1:length(d.dim.vec)]
plot(pert,apply(nc.jofc.diff.d[,,1],1,mean)/(n-m),
		xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[1],type="l")
for (d.i in 2:length(d.dim.vec)){
	lines(pert,apply(nc.jofc.diff.d[,,d.i],1,mean)/(n-m),xlab="perturbation parameter",
			ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[d.i])
}
for (d.i in 1:length(d.dim.vec)){
	lines(pert,apply(nc.cmds.d[,,d.i],1,mean)/(n-m),xlab="perturbation parameter",
			ylab="Fraction of  correct matches",ylim=c(0,1),lty=2,col=colors.vec[d.i])
}
for (d.i in 1:length(d.dim.vec)){
	lines(pert,apply(nc.jofc.weighted.d[,,d.i],1,mean)/(n-m), lty=4,lwd=2,
			xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[d.i])
}

legend.txt<- c("diff dist","CMDS","weighted.graph")
legend.txt <- c(legend.txt,d.dim.vec)
legend(x="topright",legend=legend.txt, col =c(rep(1,3),colors.vec),
		lty=c(c(1,2,4),rep(1,length(d.dim.vec))))
title("marriage o jofc")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()







nc.jofc.diff.w       = array(0,dim=c(npert,nmc,length(w.vals.vec)))
nc.jofc.weighted.w   = array(0,dim=c(npert,nmc,length(w.vals.vec)))
nc.jofc.unweighted.w = array(0,dim=c(npert,nmc,length(w.vals.vec)))



# these three lines are not used, really ...
# they just set up a dummy D.M which is the correct object for input to fullmatch".
# i'm just using Sancar's object ... i don't really understand it!
G<-matrix(0,n,n)
G.comb<-matrix(0,2*n,2*n)
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
		oos.sampling<-sample(1:n,size=n-m,replace=FALSE)
		in.sample.ind<-rep(TRUE,2*n)
		in.sample.ind[oos.sampling]<-FALSE
		in.sample.ind[n+oos.sampling]<-FALSE
		
		
		J.1 =jofc.diffusion.dist(G,Gp,
				in.sample.ind,
				d.dim=d.dim,
				w.vals.vec,
				graph.is.directed=FALSE,
				oos=oos,
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE)
		for (w.i in 1:w.max.index){
			M = solveMarriage(J.1[[w.i]])
			nc.jofc.diff.w[ipert,imc,w.i] = present(M)         # returns the number correct in the marriage?
		}
		
		J.2 =jofc(G,Gp,
				in.sample.ind,
				d.dim=d.dim,
				w.vals.vec,
				graph.is.directed=FALSE,
				oos=oos,          
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE)
		for (w.i in 1:w.max.index){
			M = solveMarriage(J.2[[w.i]])
			nc.jofc.weighted.w[ipert,imc,w.i] = present(M)         # returns the number correct in the marriage?
		}
		J.3 =jofc(G,Gp,
				in.sample.ind,
				d.dim=d.dim,
				w.vals.vec,
				graph.is.directed=FALSE,
				oos=oos,
				use.weighted.graph=FALSE,
				wt.matrix.1=NULL,
				wt.matrix.2=NULL,
				sep.graphs=TRUE)
		for (w.i in 1:w.max.index){
			M = solveMarriage(J.3[[w.i]])
			nc.jofc.unweighted.w[ipert,imc,w.i] = present(M)         # returns the number correct in the marriage?
		}
		
		
		
	}
}


### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot4.pdf")
colors.vec<-c( "red","blue","orange","black","purple")
colors.vec<- colors.vec[1:w.max.index]
plot(pert,apply(nc.jofc.diff.w[,,1],1,mean)/(n-m),xlab="perturbation parameter",
		ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[1],type="l")
for (w.i in 2:w.max.index){
	lines(pert,apply(nc.jofc.diff.w[,,w.i],1,mean)/(n-m),
			xlab="perturbation parameter",ylab="Fraction of  correct matches",
			ylim=c(0,1),lty=1,col=colors.vec[w.i])
}

for (w.i in 1:w.max.index){
	lines(pert,apply(nc.jofc.weighted.w[,,w.i],1,mean)/(n-m),
			xlab="perturbation parameter",ylab="Fraction of  correct matches",
			ylim=c(0,1),lty=3,col=colors.vec[w.i])
}
for (w.i in 1:w.max.index){
	lines(pert,apply(nc.jofc.unweighted.w[,,w.i],1,mean)/(n-m),
			xlab="perturbation parameter",ylab="Fraction of  correct matches",
			ylim=c(0,1),lty=4,col=colors.vec[w.i])
}
legend.txt<- c("diff dist","weighted.graph","unweighted.graph")
legend.txt <- c(legend.txt,w.vals.vec)
legend(x="topright",legend=legend.txt, col =c(rep(1,3),colors.vec),
		lty=c(c(1,3,4),rep(1,w.max.index)))
title("marriage o jofc")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()



T.diff.vec<- c(0.5,1,2,5,10)
nc.jofc.diff.T       = array(0,dim=c(npert,nmc,length(d.dim.vec)))




# these three lines are not used, really ...
# they just set up a dummy D.M which is the correct object for input to fullmatch".
# i'm just using Sancar's object ... i don't really understand it!
G<-ER(n,0.5)
G.comb<-omnibusM(G,G,diag(n))
Graph.M <- graph.adjacency(G.comb,weighted= NULL ,mode="undirected")
D.M<-shortest.paths(Graph.M)

seed<-123

for (T.i in 1:length(T.diff.vec)){
	T.diff<-T.diff.vec[T.i]
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
			oos.sampling<-sample(1:n,size=n-m,replace=FALSE)
			in.sample.ind<-rep(TRUE,2*n)
			in.sample.ind[oos.sampling]<-FALSE
			in.sample.ind[n+oos.sampling]<-FALSE
			
			
			J.1 =jofc.diffusion.dist(G,Gp,
					in.sample.ind,
					d.dim=d.dim,
					w.vals.vec,
					graph.is.directed=FALSE,
					oos=oos,
					wt.matrix.1=NULL,
					wt.matrix.2=NULL,
					sep.graphs=TRUE)
			
			M = solveMarriage(J.1[[1]])
			nc.jofc.diff.T[ipert,imc,T.i] = present(M)         # returns the number correct in the marriage?							
		}
	}
}

### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot5.pdf")
colors.vec<-c( "red","blue","orange","green","purple")

colors.vec<-colors.vec[1:length(T.diff.vec)]
plot(pert,apply(nc.jofc.diff.T[,,1],1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[1])
for (d.i in 2:length(d.dim.vec)){
	points(pert,apply(nc.jofc.diff.T[,,T.i],1,mean)/(n-m),xlab="perturbation parameter",ylab="Fraction of  correct matches",ylim=c(0,1),pch=1,col=colors.vec[T.i])
}
legend.txt <- c(T.diff.vec)
legend(x="topright",legend=legend.txt, 
		col =c(colors.vec),pch=rep(1,length(T.diff.vec)))
title("marriage o jofc (different diffusion param T)")
abline(h=1/(n-m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()




















