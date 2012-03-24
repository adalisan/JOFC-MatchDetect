# TODO: Add comment
# TODO: Classical MDS version of wMDS
# TODO: Read parameters from txt file
# TODO: Create Rnw file
# Author: Sancar
###############################################################################
#setwd(paste(Sys.getenv("R_HOME") ,"./../projects/",collapse=""))



source("./src/runningParams.R")



print(debug.mode)
if (par.compute){
	if (par.compute.sf){
		
	} else{
		if( run.in.linux) {
			require(doMC)
			require(foreach)
			registerDoMC(2)
			
		}
		else {require(doSMP)	
			require(foreach)
			#	setMKLthreads(1)
			workers<-startWorkers(4,FORCE=TRUE)
			registerDoSMP(workers)
		}
	}
}







#Rprof(filename = "Rprof.out", append = FALSE, interval = 0.02,
#       memory.profiling=FALSE)


Fid1.List<-list()
Fid2.List<-list()
Comm.List<-list()

avg.FCratios.List<-c()
avg.wtFCratios.List<-c()
FCratios.List<-list()
wtFCratios.List<-list()
avg.Fid1<-c()
avg.Fid2<-c()
avg.Comm<-c()

avg.Fid1.alt<-c()
avg.Fid2.alt<-c()
avg.Comm.alt<-c()

estim.wstar<-c()



if (gauss.sim)
	sim.res.g<-list()
if (dirichlet.sim)
	sim.res.d<-list()

if (profile.mode) Rprof(filename = "Rprof.out", append = FALSE, interval = 0.02,
       memory.profiling=FALSE)

for (n.v in n.vals){

for (c.val in c.vals) {
	
	#coincid.vec.dotpr.thres <-0.9
	#eigen.spectrum <- F
	#grassmannian.dist <- F
	#use.Euc.points.x <-F
	#use.Euc.points.y <-F
	#p<-6
	#r<-3
	
	#q<-14
	
	#d<-2
	
	#oos <- TRUE
	#
	# assume oos observations in different conditions are matched for oos-embedding
	#
	#assume.matched.for.oos <-TRUE
	#fp.points<- seq(0,1,0.01)
	#w.vals<-c(0.001,seq(.1,.9,0.1),0.999)
	#rival.w<- 0.999
	
	#Ignore separability error, Set weights to 0 for dissimilarities
	#between different objects under different conditions
	# Weights corresponding to off-diagonal entries of L
	#separability.entries.w<- FALSE
	#
	#Use imputed dissimilarities between in-sample and out-of sample objects of different conditions
	# weights corresponding to V_{12}, V_{21} in the omnibus matrix
	#
	#oos.use.imputed <- FALSE
	
	
	#plot.title <- ""
	#old.gauss.model <- F
	
	
	methods.vec<-c("jofc","pom","cca")
	

	run.time.g<-0
	
	params$c.val <- c.val
	
	params$n <- n.v
	if (gauss.sim){
		
		params$p<-8
		params$r<-10
		
		params$q<-12
		params$oos <- TRUE
		
		params$d<-2
		
		params.text.1 <- bquote(p==.(params$p) | r==.(params$r) |q==.(params$q) |c==.(params$c.val)|d==.(params$d))
		params.text.3 <- bquote(n==.(params$n) |assume.matched.oos==.(params$assume.matched.for.oos) | nmc==.(params$nmc)| oos==.(params$oos) )
		model.letter<-"G"
		plot.title.G <-paste(model.letter,deparse(params.text.1),"\n,Wchoice=",params$Wchoice,
				" pre.scaling = ",params$pre.scaling,#"\n",		
				deparse(params.text.3) ,"\n")
		params$plot.title <- plot.title.G
		#parameter logging start
		#sink(file=file.path('reports',c(model.letter,"c",params$c.val,"params.txt")))
		print(plot.title.G)
		print(params)
		#sink()
		#parameter logging end
		info(logger,"Gaussian Setting Simulation Starting")
		sim.res.g<-simulate.generate.test.model.plot.Kcond("MVN",params,par.compute,K)
			print("Gaussian Setting Simulation Ended")
		par(lty=1)
		
	}
} #end c.vals


} #end n.vals

if (profile.mode) Rprof(NULL)


#sink(file.path("logs","traceback.txt"))
traceback()
#sink()
warnings()

#Rprof(NULL)
if ((!run.in.linux) & par.compute &(!par.compute.sf)) 
	stopWorkers(workers)
