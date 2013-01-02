
debug.mode<-TRUE
debug.mode.fine <-FALSE
require(igraph)
require(optmatch)
require(MASS)
source("./lib/simulation_math_util_fn.R")
source("./lib/smacofM.R")
source("./lib/oosIM.R")
source("./lib/oosMDS.R")
source("./lib/diffusion_distance.R")
source("./lib/graph_embedding_fn_many.R")

cep=TRUE
verbose= FALSE
oos=TRUE
oos.cep = TRUE

a.cep <-20


n<- 100
m = 20 # number of test nodes in the second graph that  are to-be-matched
m_vals <- c(10,20,50,80,90)
m_len <- length(m_vals)

pert=(0:5)/10
pert= c(0,0.1,0.3,0.5)

npert <-  length(pert)

nc.jofc.dice.wt.p <- nc.jofc.dice.wt.r <- nc.jofc.dice.wt.f <- array(0,dim=c(npert,nmc,m_len))


nc.cmds = matrix(0,npert,nmc)

matched.cost<-0.01


#w.vals.vec <- c(0.5,0.7,0.9,0.95)
w.vals.vec <- c(0.7)

w.max.index<-length(w.vals.vec)

matched.cost<-0.01 #If matched.cost is equal to 1, consider an unweighted graph, with edges between matched vertices
#If matched.cost is  between 0 and 1, the graph is weighted with edges between matched vertices with weights equal to matched.cost. Edges between 
# vertices of the same condition have either weight 1 or 2 according to whether they're connected according to given adjacency matrix or not.

d.start <- 12
T.diff<-2

dims.for.dist <- 1:d.dim

seed<-123


gen.1.to.k.matched.graphs <- function(n,pert,repeat.counts) {
  npert <-  length(pert)
	G.orig<-ER(n,0.5)
	diag(G.orig)<-1
	k<-1
  
	int.end.indices<-cumsum(repeat.counts)
	int.start.indices<-c(1,int.end.indices+1)
	corr.list<-vector("list",n)
  G<-matrix(0,new.n <- sum(repeat.counts),new.n)
	for (i in 1:n){
		for (j in 1:repeat.counts[i]){
			G[k,]<-rep(G.orig[i,],times=repeat.counts)
			G[,k]<-rep(G.orig[i,],times=repeat.counts)
			#	G<-perturbG(G,0.1)
			k <- k+1
		}
		corr.list[[i]] <- list(a=(int.start.indices[i]:int.end.indices[i]),b=i)
	}
	
	diag(G.orig)<-0
	diag(G)<-0
  G.list <- list()
  Gp.list <- list()
	for(ipert in 1:npert)
	{
    
    #Gp<-bitflip(G.orig ,pert[ipert],pert[ipert])
    Gp<-G.orig
    Gp.list<-c(Gp.list,list(Gp))
    G.t<-bitflip(G ,pert[ipert],pert[ipert])
    G.list<-c(G.list,list(G.t))
  }
  return(list(G=G.list,Gp.list=Gp.list,   corr.list=corr.list,
              int.start.indices = int.start.indices,
              int.end.indices = int.end.indices))
}



nmc <- 25



for(imc in 1:nmc)
{
  repeat.counts <-1+rgeom(n,0.2)
  repeat.counts[repeat.counts>10]=10;
  #  repeat.counts <-rep(1,n)
  new.n <- sum(repeat.counts)
  
  gen.graph.pair <- gen.1.to.k.matched.graphs(n,pert,repeat.counts)
  G.list <- gen.graph.pair$G
  Gp.list <- gen.graph.pair$Gp.list
  corr.list <- gen.graph.pair$corr.list
  int.start.indices <- gen.graph.pair$int.start.indices
  int.end.indices   <- gen.graph.pair$int.end.indices
  
  for (m_it in 1:m_len) {
    m<- m_i <- m_vals[m_it]
    
    oos.sampling<-sample(1:n, size=m_i, replace=FALSE)
		in.sample.ind.1<-rep(TRUE,new.n)
		for ( s in 1:m){
			a<-int.start.indices[oos.sampling[s]]
			b<-int.end.indices[oos.sampling[s]]
			in.sample.ind.1[a:b]<-FALSE
		}
		
		in.sample.ind.2<-rep(TRUE,n)
		in.sample.ind.2[oos.sampling]<-FALSE
		
		#if (imc==1) print(in.sample.ind)
		
    for (pert_i in  1:npert) {
      Gp <- Gp.list[[pert_i]]
      G.rep <- G.list[[pert_i]]
      J.1 =JOFC.graph.custom.dist.many (G.rep, Gp, corr.list,
				in.sample.ind.1,in.sample.ind.2,
                                        d.dim=d.start,
                                        w.vals.vec=w.vals.vec,
				graph.is.directed=FALSE,
                                        vert_diss_measure  =  'C_dice_weighted',
                                        T.param  =  NULL,
		
                                        graph.is.weighted=TRUE)
		
      #print(head(J.1[[1]]))
      #print(diag(J.1[[1]]))
      M = solveMarriage.many(J.1[[1]],10)
      match.perf.eval <- present.many(M,corr.list)
      nc.jofc.dice.wt.p[pert_i,imc,m_it] = mean(match.perf.eval$P)
      nc.jofc.dice.wt.r[pert_i,imc,m_it] = mean(match.perf.eval$R)
      nc.jofc.dice.wt.f[pert_i,imc,m_it] <- mean(match.perf.eval$F)
      print(dim(nc.jofc.dice.wt.p))
      if (pert_i>1 && imc>1){
        print("Precision")
        print(apply(nc.jofc.dice.wt.p[,1:imc,m_it],1,mean))
        print("Recall")
        print(apply(nc.jofc.dice.wt.r[,1:imc,m_it],1,mean))
        print("F-measure")
        print(apply(nc.jofc.dice.wt.f[,1:imc,m_it],1,mean))
      }
      else{
        print("Precision")
        print(nc.jofc.dice.wt.p[pert_i,imc,m_it])
        print("Recall")
        print(nc.jofc.dice.wt.r[pert_i,imc,m_it])
        print("F-measure")
        print(nc.jofc.dice.wt.f[pert_i,imc,m_it])
      } 
		
    }
		
				
  }
}
		
		
		
pdf("./graphs/plot-PRF2.pdf")
colors.vec<-c( "red","blue","orange","green")
par(lty=1)
for (m_v_i in 1:m_len){
par(lty=m_v_i)
plot.graph.with.CI.error(t(nc.jofc.dice.wt.f[,,m_v_i]),plot.title="",plot.col=colors.vec[1],
                         conf.int=TRUE,add=m_v_i>1,fp.points=pert,
                         customx.labels=NULL,customy.labels=NULL,
                         ispowercurve=FALSE,ylim=c(0,1),xlab="perturbation parameter",
                         ylab="Average F-measure/Precision/Recall")
par(lty=m_v_i)
plot.graph.with.CI.error(t(nc.jofc.dice.wt.p[,,m_v_i]),plot.title="",plot.col=colors.vec[3],
                         conf.int=TRUE,add=TRUE,fp.points=pert,
                         customx.labels=NULL,customy.labels=NULL,
                         ispowercurve=FALSE,ylim=c(0,1),xlab="perturbation parameter",
                         ylab="Average F-measure/Precision/Recall")
par(lty=m_v_i)
plot.graph.with.CI.error(t(nc.jofc.dice.wt.r[,,m_v_i]),plot.title="",plot.col=colors.vec[4],
                         conf.int=TRUE,add=TRUE,fp.points=pert,
                         customx.labels=NULL,customy.labels=NULL,
                         ispowercurve=FALSE,ylim=c(0,1),xlab="perturbation parameter",
                         ylab="Average F-measure/Precision/Recall")
}
legend.txt<- c("F-measure","Precision","Recall")
title("1-to-k matching o jofc")
abline(h=1/(m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
legend.txt<-c(legend.txt,m_vals)
ltypes<- c(rep(1,3),1:m_len)
legend(x="topright",legend=legend.txt, col =colors.vec[c(c(1,3,4),rep(1,m_len))],lty=ltypes)
		
		
				
dev.off()
		
pdf("./graphs/F-measure_1_to_k_match.pdf")
colors.vec<-c( "red","blue","orange","green","magenta")
par(lty=1)
for (m_v_i in 1:m_len){
#  par(lty=m_v_i)
  plot.graph.with.CI.error(t(nc.jofc.dice.wt.f[,,m_v_i]),plot.title="",plot.col=colors.vec[m_v_i],
                           conf.int=TRUE,add=m_v_i>1,fp.points=pert,
                           customx.labels=NULL,customy.labels=NULL,
                           ispowercurve=FALSE,ylim=c(0,1),xlab="perturbation parameter",
                           ylab="Average F-measure")
	}
  
legend.txt<- c()
title("1-to-k matching o jofc")
abline(h=1/(m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
legend.txt<-c(legend.txt,n-m_vals)
ltypes<- c(1:m_len)
legend(x="topright",legend=legend.txt, col =colors.vec[1:m_len],lty=1,title="Num. of Hard Seeds")

dev.off()


### notice that the %correctmatches decreases as perturbation parameter increases! :-)

pdf("plot1.pdf")
colors.vec<-c( "red","blue","orange","green")
colors.vec<-colors.vec[1:4]
plot(pert,apply(nc.jofc.dice.wt.f,1,mean),xlab="perturbation parameter",ylab="Average F-measure",ylim=c(0,1),col=colors.vec[1])
#points(pert,apply(nc.cmds,1,mean)/(n-m),xlab="perturbation parameter",ylab="Average F-measure",ylim=c(0,1),pch=2,col=colors.vec[2])

points(pert,apply(nc.jofc.weighted.f,1,mean),xlab="perturbation parameter",ylab="Average F-measure",ylim=c(0,1),pch=3,col=colors.vec[3])

points(pert,apply(nc.jofc.unweighted.f,1,mean),xlab="perturbation parameter",ylab="Average F-measure",ylim=c(0,1),pch=4,col=colors.vec[4])

legend.txt<- c("diff dist","weighted.graph","unweighted.graph")
legend(x="topright",legend=legend.txt, col =colors.vec,pch=1:4)
title("1-to-k matching o jofc")
abline(h=1/(m),lty=2) ### chance?  apparently not!?
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()



pdf("plot2.pdf")
colors.vec<-c( "red","blue","orange","green")
colors.vec<-colors.vec[1:4]
plot(pert,apply(nc.jofc.dice.wt.p,1,mean),xlab="perturbation parameter",
		ylab="Average precision/recall",ylim=c(0,1),col=colors.vec[1],type="l")
#points(pert,apply(nc.cmds,1,mean)/(n-m),xlab="perturbation parameter",
#     ylab="Average precision/recall",ylim=c(0,1),pch=2,col=colors.vec[2])


lines(pert,apply(nc.jofc.weighted.p,1,mean),xlab="perturbation parameter",
		ylab="Average precision/recall",ylim=c(0,1),col=colors.vec[3])

lines(pert,apply(nc.jofc.unweighted.p,1,mean),xlab="perturbation parameter",
		ylab="Average precision/recall",ylim=c(0,1),col=colors.vec[4])


lines(pert,apply(nc.jofc.dice.wt.r,1,mean),xlab="perturbation parameter",
		ylab="Average precision/recall",ylim=c(0,1),lty=2,col=colors.vec[1])


lines(pert,apply(nc.jofc.weighted.r,1,mean),xlab="perturbation parameter",
		ylab="Average precision/recall",ylim=c(0,1),lty=2,col=colors.vec[3])

lines(pert,apply(nc.jofc.unweighted.r,1,mean),xlab="perturbation parameter",
		ylab="Average precision/recall",ylim=c(0,1),lty=2,col=colors.vec[4])


legend.txt<- c("diff dist","weighted.graph","unweighted.graph")
legend(x="topright",legend=legend.txt, col =colors.vec,pch=1:4)
title("marriage o jofc")
abline(v=1/2,lty=2) ### chance?  gotta be!?
dev.off()





