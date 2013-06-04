# TODO: Add comment
# 
# Author: Sancar
###############################################################################

debug.mode<-F
run.for.Sweave<-F


if (!run.for.Sweave){
  source("./src/runningParams.R")
} else {
  #	sink(paste("run at ",Sys.Date()))
}
print(debug.mode)
if (par.compute){
  
  
  
  
  num.cores<-parallel::detectCores()
  iter_per_core <- ceiling(num_iter/num.cores)
  require(foreach)
  
  
  if (!par.compute.sf && .Platform$OS.type != "windows" && require("multicore")) {
    require(doMC)
    registerDoMC(parallel::detectCores())
  } else if (!par.compute.sf && FALSE &&                     # doSMP is buggy
               require("doSMP")) {
    workers <- startWorkers(num.cores,FORCE=TRUE) # My computer has 4 cores
    on.exit(stopWorkers(workers), add = TRUE)
    registerDoSMP(workers)
  } else if (require("doSNOW")) {
    cl <- snow::makeCluster(num.cores, type = "SOCK")
    on.exit(snow::stopCluster(cl), add = TRUE)
    registerDoSNOW(cl)
  } else {
    registerDoSEQ()
  }
  
  
  
} else {
  registerDoSEQ()
}

model <- NULL   
if ( gauss.sim) {
  model="MVN"
} else if (dirichlet.sim) {
  model="Dirichlet"
} else{
  stop("invalid model name")}

num.cpus<-num.cores

gaussian_simulation_jofc_tradeoff
gaussian_simulation_jofc_tradeoff_par
gaussian_simulation_jofc_tradeoff_sf
dirichlet_simulation_jofc_tradeoff
dirichlet_simulation_jofc_tradeoff_par
dirichlet_simulation_jofc_tradeoff_sf


if ((model=="MVN") && par.compute)              call.func<-gaussian_simulation_jofc_tradeoff_par
if ((model=="MVN") && !par.compute)          call.func<-gaussian_simulation_jofc_tradeoff
if ((model=="Dirichlet") && par.compute)        call.func<-dirichlet_simulation_jofc_tradeoff_par
if ((model=="Dirichlet") && !par.compute)    call.func<-dirichlet_simulation_jofc_tradeoff
if ((model=="MVN")       && par.compute.sf)           call.func<-gaussian_simulation_jofc_tradeoff_sf
if ((model=="Dirichlet") && par.compute.sf)     call.func<-dirichlet_simulation_jofc_tradeoff_sf
sim.res<-list()

if (model=="MVN") real.dim<- params$p.g+params$q.g
if (model=="Dirichlet") real.dim <- params$p.d+params$q.d+2

params$compare.pom.cca<-TRUE


if (model=="MVN"){
  params$p <-params$p.g
  params$q <- params$q.g
  params$r <- params$r.g
} else{
  params$p <- params$p.d
  params$q <- params$q.d
  params$r <- params$r.d
}








if (run.for.Sweave) print("c")

if (run.for.Sweave) print(c.val)
if (run.for.Sweave) print(p)
if (run.for.Sweave) print("w values")
if (run.for.Sweave)  print(w.vals)














begin.time.g <-Sys.time()
args.for.func.call<-with(params,list(p=p, r=r, q=q, 
                                     c.val=c.val,d=d,
                                     
                                     Wchoice     = "avg", 
                                     pre.scaling = TRUE,
                                     oos         = oos,
                                     alpha       = NULL,
                                     n = n, m = s, nmc = nmc,
                                     
                                     old.gauss.model.param=old.gauss.model,
                                     separability.entries.w=separability.entries.w,
                                     compare.pom.cca=compare.pom.cca,
                                     oos.use.imputed=oos.use.imputed,
                                     
                                     
                                     rival.w=rival.w,
                                     proc.dilation=FALSE,
                                     assume.matched.for.oos =assume.matched.for.oos,
                                     w.vals=w.vals,
                                     wt.equalize=wt.equalize,
                                     verbose=verbose, 
                                     power.comparison.test=power.comparison.test,
                                     cca.reg=TRUE))


if (!par.compute.sf)
  args.for.func.call<- c(args.for.func.call,list(pprime1= real.dim, pprime2= real.dim))





bootstrap.res <- with(params, do.call(test.bootstrapped.JOFC,args=c(list(model="gaussian"),args.for.func.call)))


############################################
#  Running simulation function 
############################################
sim.res <- do.call(call.func,args=args.for.func.call)

end.time.g<-Sys.time()
run.time.g <- end.time.g-begin.time.g

print(run.time.g)

save.image(file= paste("JOFC_MVN_Dir_Sim_",format(Sys.time(), "%b %d %H:%M:%S"),".RData"))
write.csv(sim.res$auc.meas,"auc_meas.txt")




source("./lib/JOFC-MCSim-plot_fn.R")
#draw.plots(sim.res,"MVN",params,plot.w.vals=1:length(params$w.vals),TRUE,TRUE,FALSE)

#source("./src/color-setting.R")


w.vals <- params$w.vals

#+
#guides(colour=guide_legend( title =expression(rho),title.hjust=1,title.vjust=-1,label.hjust=1))



draw.plots(sim.res,"MVN",params,plot.w.vals=c(1:6,9,12,14,15),FALSE,FALSE,FALSE)
draw.auc.error.bar.argmax.bar(w.vals)

auc.avg <-apply(auc.meas,2,mean)
auc.sd <-apply(auc.meas,2,sd)
auc.mat<- matrix(auc.avg,1,length(w.vals))
auc.mat<- rbind(auc.mat,auc.sd)
rownames(auc.mat)<-c("mean","SE")
colnames(auc.mat) <- w.vals
xtable(auc.mat,digits=4,caption= "mean and standard error of  AUC($w$) for $nmc=400$ MC replicates" ,label= "tab:AUCW")


