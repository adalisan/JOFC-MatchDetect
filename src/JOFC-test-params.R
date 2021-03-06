# TODO: Add comment
# 
# Author: Sancar
###############################################################################



print(debug.mode)
if (par.compute){
	if (par.compute.sf){
		
	} else{
		if( run.in.linux) {
			require(doMC)
			require(foreach)
			registerDoMC(2)
			
		}
		else {
			require(doSMP)	
			require(foreach)
			#	setMKLthreads(1)
			workers<-startWorkers(4,FORCE=TRUE)
			registerDoSMP(workers)
		}
	}
}

model="MVN"


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


call.func <- test_params
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


params$compare.pom.cca<-TRUE


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

############################################
#  Running simulation function 
############################################
sim.res <- do.call(call.func,args=args.for.func.call)




draw.plots<-function(sim.res,model,params){

#
# Plotting of ROC curves
#
sim.res.pow<-sim.res$power
w.val.len <-length(params$w.vals)

lty.i.vec<-c()
for (i in 1:w.val.len){
	lty.i <- 1+((i-1)%%10)
	
	lty.i.vec <- c(lty.i.vec,lty.i)
	par(lty=lty.i)
	plot.ROC.with.CI(sim.res.pow[i,,],plot.title="",plot.col = colors.vec[i],
			conf.int=FALSE,add=(i>1),ylim=1)
	
}

if (compare.pom.cca) {
	i<- 1+((w.val.len)%%10)
	par(lty=i)
	plot.ROC.with.CI(sim.res$power.cmp$pom,plot.title="",plot.col = colors.vec[w.val.len+1],
			conf.int=FALSE,add=TRUE,ylim=1)
	lty.i.vec <- c(lty.i.vec,i)			
	i<- 1+((w.val.len+1)%%10)
	par(lty=i)
	plot.ROC.with.CI(sim.res$power.cmp$cca,plot.title="",plot.col = colors.vec[w.val.len+2],
			conf.int=FALSE,add=TRUE,ylim=1)
	lty.i.vec <- c(lty.i.vec,i)	
	if (cca.reg){
		i<- 1+((w.val.len+2)%%10)
		par(lty=i) 
		plot.ROC.with.CI(sim.res$power.cmp$reg.cca,plot.title="",plot.col = colors.vec[w.val.len+3],
				conf.int=FALSE,add=TRUE,ylim=1)
		lty.i.vec <- c(lty.i.vec,i)	
		#lwd.i.vec <- c(lwd.i.vec,line.width)
	}
	
	
}
if (compute.bound) {
	i<- 1+((w.val.len+3)%%10)
	par(lty=i)
	lines(sim.res$optim.power,col=colors.vec[w.val.len+3])
	lty.i.vec <- c(lty.i.vec,i)
}

legend.txt <- params$w.vals
if (compare.pom.cca)
	legend.txt <-c(legend.txt ,"pom","cca")
if (compute.bound) legend.txt <-c(legend.txt ,"bound")
if (cca.reg)
	legend.txt <-c(legend.txt ,"cca.reg")
if (compute.bound) {legend.txt <-c(legend.txt ,"bound")}

legend("bottomright",legend=legend.txt,
		col=colors.vec,lty=lty.i.vec)
title(plot.title)
par(lty=1)
fname<- file.path('graphs',paste(c(model,"-FC-Tradeoff-",ifelse(oos,"OOS","noOOS"),"c",c.val),collapse=""))


if(!run.in.linux&(!run.for.Sweave))  savePlot(paste(fname,".pdf",sep="",collapse=""),type="pdf")
if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste(fname,".ps",sep="",collapse=""),type="ps")
if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste(fname,".png",collapse="",sep=""),type="png")

dev.off()






#
#  Power vs w plotting

#
if( run.in.linux) X11() else {windows()}


#temp
par(lty=1)
lwd.old <- par()$lwd
par(lwd=3)
#
# Power vs w plots for fixed alpha
avg.power.w <-rep(0,length(w.val.len))
for (i in 1:w.val.len){
	beta.w<-sim.res$power[i,,2]
	avg.power.w[i]<-mean(beta.w)		
}


#print(which.max(avg.power.w))
#print(params$w.vals)	

#w
max.w <- which.max(avg.power.w)
sim.res$wstar.estim <- params$w.vals [max.w]
if (verbose) print("Estimate of wstar for average power curve")
if (verbose) print(sim.res$wstar.estim)
#sink("debug-wstar.txt")
#print(sim.res$power[,,6])

#Which w value was the best w for the highest number of mc replicates
#Note that multiple w.values might have the best power
sim.res$wstar.idx.estim.mc<- apply(sim.res$power[,,6],2,function(x) which(x==max(x)))

print(sim.res$wstar.idx.estim.mc)
#sink()
avg.power.w.2 <-rep(0,length(w.val.len))
for (i in 1:w.val.len){
	beta.w<-sim.res$power[i,,6]
	avg.power.w.2[i]<-mean(beta.w)
}

avg.power.w.3 <-rep(0,length(w.val.len))
for (i in 1:w.val.len){
	beta.w<-sim.res$power[i,,11]
	avg.power.w.3[i]<-mean(beta.w)			
}

x.vals<- 1:w.val.len

plot(x=x.vals,y=avg.power.w,col="blue",type="l",ylim=c(0,1),log='x',xlab=expression(w),ylab=expression(beta),xaxt="n",ps=10)
lines(x=x.vals,y=avg.power.w.2,col="red",xlog=TRUE)
lines(x=x.vals,y=avg.power.w.3,col="green",xlog=TRUE)
par(cex.axis=0.9)
axis(side=1,at=1:length(params$w.vals), labels = params$w.vals)
legend("bottomright",expression(alpha==0.01,alpha==0.05,
				alpha==0.1),col=c("blue","red","green"),lty=rep(1,3))


#title(paste(model," model: power vs w plot"))
fname <- file.path(paste(results.dir,ifelse(oos,"OOS","noOOS"),model,"-power-w-c",c.val,sep="", collapse=""))
if(!run.in.linux&(!run.for.Sweave))	savePlot(paste(fname,".pdf",sep="",collapse=""),"pdf")
if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(fname,".png",sep="",collapse=""),"png")

if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(fname,".ps",sep="",collapse=""),"ps")

par(lwd=lwd.old)

#dev.print(png,file=fname,height=600,width=600)
dev.off()



}











