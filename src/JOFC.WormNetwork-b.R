#GRAPH TYPES
# 1. a)Undirected b)Directed
# 2. a)JOFC  embed   b) CMDS  embed
# 3. a)From Adjacency Matrix to Unweighted Minimum Distance 
#    b)From Adjacency Matrix to Weighted Minimum Distance
#    c) From Weight Matrix
# 4.  a) sep.graphs= TRUE, CEP's way oftreat two graphs separately to compute dissimilarities
# and impute W (off-diagonalblock matrix)
# b) sep.graphs= FALSE  join the graphs and compute dissimilarities from joint graph


require(igraph)
require(optmatch)
source("./lib/simulation_math_util_fn.R")
source("./lib/smacofM.R")
source("./lib/oosIM.R")
source("./lib/graph_embedding_fn.R")
load("./data/celegansGraph.Rd")
cep=FALSE
verbose= FALSE

Ac<-Ac*upper.tri(Ac)
Ag<-Ag*upper.tri(Ag)
Ac<-(Ac+t(Ac))
Ag<-(Ag+t(Ag))
diag(Ac)<-0
diag(Ag)<-0
Ac.w<- Ac
Ag.w<- Ag


Ac<- ifelse(Ac>0,1,0)
Ag<- ifelse(Ag>0,1,0)

n = nrow(Ag)
m = n-10 # the first m pairs are known matches ; the last n-m pairs are to-be-matched

T.diff<-2
d<-4
npert = 11
nmc = 100
pert=(0:10)/10

nc.worms.jofc<-rep(0,nmc)
w.vals<-0.8

matched.cost<-0.01 #If matched.cost is equal to 1, consider an unweighted graph, with edges between matched vertices
#If matched.cost is  between 0 and 1, the graph is weighted with edges between matched vertices with weights equal to matched.cost. Edges between 
# vertices of the same condition have either weight 1 or 2 according to whether they're connected according to given adjacency matrix or not.
debug.mode<-FALSE

d.dim<-10




run.jofc.replicate.batch <- function(i, mc.rep.per.batch,method="jofc"){

      num.matches<-rep(0,mc.rep.per.batch)

      for ( j in 1:mc.rep.per.batch) {
	oos.sampling<-sample.int(n,size=n-m,replace=FALSE)
	in.sample.ind<-rep(TRUE,2*n)
	in.sample.ind[oos.sampling]<-FALSE
	in.sample.ind[n+oos.sampling]<-FALSE
	
	   if (method=="jofc"){
	
		J = jofc(Ac,Ag,in.sample.ind,d.dim=d.dim,w.vals.vec=w.vals,
				wt.matrix.1=NULL,wt.matrix.2=NULL,
				use.weighted.graph=TRUE,
				sep.graphs=TRUE) 
            } else if (method=="jofc.wt"){
	
		J = jofc(Ac.w,Ag,in.sample.ind,d.dim=d.dim,w.vals.vec=w.vals,
				wt.matrix.1=Ac.w,wt.matrix.2=Ag.w,
				use.weighted.graph=TRUE,
				sep.graphs=TRUE) 
            } else if (method=="jofc.diff.dist"){
		J = jofc.diffusion.dist(Ac.w,Ag.w,in.sample.ind,d.dim=d.dim,,w.vals.vec=w.vals,
				wt.matrix.1=Ac.w,wt.matrix.2=Ag.w,
				sep.graphs=TRUE) 

		}

		M = solveMarriage(J[[1]])
		num.matches[j] = present(M)  
       }
      return (num.matches)
}

#require(doSMP)
#workers <- startWorkers(num.cpus) 
#registerDoSMP(workers)
 times <- 4	# times to run the loop
run.per.batch <-12 

#run.results<-foreach(run.mc=1:(nmc/run.per.batch),.combine=c) %dopar% run.jofc.replicate.batch(run.mc,run.per.batch,method="jofc")


sfInit( parallel=TRUE, cpus=num.cpus )




sfClusterSetupRNG()   
sfLibrary (igraph)
sfLibrary(optmatch)
sfSource("./lib/simulation_math_util_fn.R")
sfSource("./lib/smacofM.R")
sfSource("./lib/oosIM.R")
sfSource("./lib/graph_embedding_fn.R")
sfSource("./lib/diffusion_distance.R")

#p.prime.cond = p+q


print("Starting parallelization in gaussian_simulation_jofc_tradeoff_sf") 
run.results <- sfLapply( 1:nmc, run.jofc.replicate.batch,mc.rep.per.batch=run.per.batch,method="jofc")

sfStop()


# stop workers
#stopWorkers(workers)


