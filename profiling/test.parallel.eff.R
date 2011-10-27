c.val <-0.01
p<-10
r<-30
q<-10
nmc<-2
run.in.linux <-FALSE
verbose=TRUE
nonpar.time <- system.time(
gaussian_simulation_jofc_tradeoff(
    p, r, q, c.val,
  	d           = p-1,
		pprime1     = p+q,   # cca arguments
		pprime2     = p+q,   # cca arguments
		Wchoice     = "avg", 
		pre.scaling = TRUE,
		oos         = TRUE,
		alpha       = NULL,
		n = 100, m = 100, nmc =nmc,
		sim.grass=FALSE,
		eigen.spectrum=FALSE,
		old.gauss.model.param=FALSE,
		separability.entries.w=FALSE,
		compare.pom.cca=TRUE,
		oos.use.imputed=FALSE,
		level.mcnemar=0.01,
		def.w=0.5,
		rival.w=NULL,
		proc.dilation=FALSE,
		assume.matched.for.oos=TRUE,
		w.vals=c(0.5,0.9),
		wt.equalize=FALSE,
		verbose) 
)


if (par.compute){
  if( run.in.linux) {
		require(doMC)
		require(foreach)
		registerDoMC(2)
		
	}
	else {require(doSMP)	
		require(foreach)
		setMKLthreads(1)
		workers<-startWorkers(4,FORCE=TRUE)
		registerDoSMP(workers)
	}
}






par.time <-system.time(
gaussian_simulation_jofc_tradeoff_par (
    p, r, q, c.val,
    d           = p-1,
		pprime1     = p+q,   # cca arguments
		pprime2     = p+q,   # cca arguments
		Wchoice     = "avg", 
		pre.scaling = TRUE,
		oos         = TRUE,
		alpha       = NULL,
		n = 100, m = 100, nmc =nmc,
		
		old.gauss.model.param=FALSE,
		separability.entries.w=FALSE,
		compare.pom.cca=TRUE,
		oos.use.imputed=FALSE,
		level.mcnemar=0.01,
		def.w=0.5,
		rival.w=NULL,
		proc.dilation=FALSE,
		assume.matched.for.oos=TRUE,
		w.vals=c(0.5,0.9),
		wt.equalize=FALSE,
		verbose)
)

#Rprof(NULL)
if ((!run.in.linux) & par.compute) 
  stopWorkers(workers)



print("non-parallel version")
print(nonpar.time)

print("parallel version")
print(par.time)

