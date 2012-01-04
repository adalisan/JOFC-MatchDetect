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

Ac<-(Ac+t(Ac))/2
Ag<-(Ag+t(Ag))/2
diag(Ac)<-0
diag(Ag)<-0
Ac.w<- Ac
Ag.w<- Ag


Ac<- ifelse(Ac>0,1,0)
Ag<- ifelse(Ag>0,1,0)

n = nrow(Ag)
m = n-10 # the first m pairs are known matches ; the last n-m pairs are to-be-matched

T.diff<-100
d<-4
npert = 11
nmc = 100
pert=(0:10)/10
nc.jofc = matrix(0,npert,nmc)
nc.cmds = matrix(0,npert,nmc)
matched.cost<-0.01
nc.worms.jofc<-rep(0,nmc)
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
d<-4



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
	oos.sampling<-sample.int(n,size=n-m,replace=FALSE)
	in.sample.ind<-rep(TRUE,2*n)
	in.sample.ind[oos.sampling]<-FALSE
	in.sample.ind[n+oos.sampling]<-FALSE
	
	
	
		J = jofc(Ac.w,Ag.w,in.sample.ind,d.dim=8,
				wt.matrix.1=Ac.w,wt.matrix.2=Ag.w,
				use.weighted.graph=TRUE,
				sep.graphs=TRUE,
				use.diff.distance=TRUE) 
		M = solveMarriage(J)
		nc.worms.jofc[imc] = present(M)         # returns the number correct in the marriage?
		
}

nc.worms.jofc
### notice that the %correctmatches decreases as perturbation parameter increases! :-)

