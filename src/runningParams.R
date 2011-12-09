
results.dir <- "./graphs/"



run.for.Sweave <-TRUE
par.compute <- FALSE
run.in.linux<- (.Platform$OS.type=="unix")
compare.pom.cca<-TRUE
#c.vals<-c(0.01)
c.vals<-c(0.01,0)

verbose<-TRUE
oos <-TRUE


gauss.sim <-T
dirichlet.sim <- T
run.mcnemars.test <- F
cca.reg <- F
power.comparison.test<- F
add.plot.title<- F



n.vals <-c(50,100,150,200,300,500)
#
#n.vals<-c(150)
nmc <-  150
s<- 150

profile.mode <- FALSE
debug.mode<- FALSE

if(debug.mode){
  n.vals<-c(70)
	nmc<-3
	s<-10
}


params<-list(
  run.for.Sweave = run.for.Sweave,
  	nmc=nmc,
		coincid.vec.dotpr.thres =0.9,
		eigen.spectrum =F,
		grassmannian.dist = F,
		use.Euc.points.x = F,
		use.Euc.points.y = F,
		p = 3,
		r = 10,
		
		q = 15,
		
		d=2,
		n=n.vals[1],
		s=s,
		
		
		Wchoice="NA+diag(0)",
		pre.scaling=TRUE,
		oos =oos,
#
# assume oos observations in different conditions are matched for oos-embedding
#
		assume.matched.for.oos = TRUE,
		fp.points = seq(0,1,0.01),
		
		
#Ignore separability error, Set weights to 0 for dissimilarities
#between different objects under different conditions
# Weights corresponding to off-diagonal entries of L
		separability.entries.w = FALSE,
#
#Use imputed dissimilarities between in-sample and out-of sample objects of different conditions
# weights corresponding to V_{12}, V_{21} in the omnibus matrix
#
		oos.use.imputed = FALSE,
		proc.dilation=FALSE,
		
		plot.title = "",
		old.gauss.model = F,
		verbose=verbose,
		c.val=0,
		#w.vals = c(0.001,0.1,0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999),
		w.vals = c(0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999),
    #w.vals= c(0.5,0.999),
		wt.equalize=FALSE,
		rival.w = 0.999,
		power.comparison.test = F
)







