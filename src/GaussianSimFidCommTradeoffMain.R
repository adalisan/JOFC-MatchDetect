# TODO: Add comment
# TODO: Classical MDS version of wMDS
# TODO: Read parameters from txt file
# TODO: Create Rnw file
# Author: Sancar
###############################################################################
#setwd(paste(Sys.getenv("R_HOME") ,"./../projects/",collapse=""))

#setwd("./DataFusion_Priebe/")





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
		color.file<- read.csv("./data/Cat_12.csv",header=FALSE,skip=2,as.is=TRUE,sep=";")
		tmp.col<-dim(color.file)[2]
		rgb.vals<-as.matrix(color.file[,(tmp.col-2):tmp.col])
		colors.vec.alt<- apply(rgb.vals,1,function(x) (rgb(x[1],x[2],x[3],255,maxColorValue=255)))
		colors.vec <- c("red","green","aquamarine","purple",colors.vec.alt[1],
				"darkblue",colors.vec.alt[7],"salmon","rosybrown","magenta","orange",
				"darkorange4")
		colors.vec<-colors.vec.alt
		colors.vec[3]<-"gold4"
		colors.vec[2]<-"darkblue"
		colors.vec[4]<-"darkorange4"
		colors.vec[9]<-"red"
		colors.vec.len<-length(colors.vec)
		colors.vec[colors.vec.len+1]<-"cornflowerblue"
		colors.vec[colors.vec.len+2]<-"azure3"
		colors.vec.len<-length(colors.vec)
		run.time.g<-0
		
		params$c.val <- c.val
		
		params$n <- n.v
		if (gauss.sim){
			
			params$p <- 3
			params$r <- 10
			
			params$q <- 15
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
			sink(paste(results.dir,model.letter,"c",params$c.val,"params",Sys.Date(),".txt"))
			print(plot.title.G)
			print(params)
			sink()
			#parameter logging end
			print("Gaussian Setting Simulation Starting")
			sim.res.g<-simulate.generate.test.model.plot("MVN",params,par.compute)
			print("Gaussian Setting Simulation Ended")
			
			Fid1<-sim.res.g$FidComm.Terms$F1
			Fid2<-sim.res.g$FidComm.Terms$F2
			Comm<-sim.res.g$FidComm.Terms$C
			
			Fid1.List <- c(Fid1.List,list(Fid1))
			Fid2.List <-c(Fid2.List,list(Fid2))
			Comm.List <-c(Comm.List,list(Comm))
			
			avgFid1.n<-colMeans(Fid1)
			avgFid2.n<-colMeans(Fid2)
			avgComm.n<-colMeans(Comm)
			
			avg.Fid1.alt<-rbind(avg.Fid1.alt,avgFid1.n)
			avg.Fid2.alt<-rbind(avg.Fid2.alt,avgFid2.n)
			avg.Comm.alt<-rbind(avg.Comm.alt,avgComm.n)
			
			FC.ratio.n<-colMeans(sim.res.g$FC.ratios$f.c)
			wtFC.ratio.n<-colMeans(sim.res.g$FC.ratios$wtf.c)
			
			
			#avg.FCRatios.List<-rbind(FCRatios.List,FC.ratio.n)
			#avg.wtFCRatios.List<-rbind(wtFCRatios.List,wtFC.ratio.n)
			estim.wstar<-c(estim.wstar,sim.res.g$wstar)	
			
			
			
			
			
			
			save.image(paste(results.dir,"MVN","-c",params$c.val,"-n",params$n,Sys.Date(),".RData",collapse=""))
			
			avg.cont.table<- (sim.res.g$conting.table+0.001)/nmc
			print("Aggregate Cont Table: ")
			print(avg.cont.table)
			p.val <-mcnemar.test(avg.cont.table)$p.value
			print("Aggregate Cont Table: McNemar's p-val")
			print(p.val)
			
			sink(paste(results.dir,model.letter,"-n",params$n,"c",params$c.val,"sign-tests.txt",collapse=""))
			sign.test <- try(sign.test.cont.table(sim.res.g$conting.table.list))
			
			if (inherits(sign.test,"try-error")) {
				print(paste("error in ",model,collapse=""))
				print(sim.res.g$conting.table.list)
			}		else{
				print("Cont Table List: sign test p-val")
				print(sign.test$p.value)
			}
			
			
			
			sign.rank.sum.test<-sign.rank.sum.test.cont.table(sim.res.g$conting.table.list)
			print("Cont Table List: signed rank sum test p-val")
			print(sign.rank.sum.test$p.value)
			sink()
			
			
			
			
			
			
			
			conting.table.list.g <- lapply(sim.res.g$conting.table.list,function(x) x+0.001)
			
			p.vals.list <- lapply(conting.table.list.g,mcnemar.test)
			p.vals<-c()
			for (t in p.vals.list)
				p.vals<- c(p.vals,t$p.value)
			p.vals <- sort(p.vals)
			print("MC Replicate Cont Tables: McNemar's p-val")
			
			
			save.image(paste(results.dir,"MVN","c",params$c.val,"-n",params$n,Sys.Date(),".RData",collapse=""))
			
			
			
			windows()
			#colors.vec[3]<-"blue"
			lty.i.vec<-c()
			lwd.i.vec <- c()
			line.width<-1.5
			w.val.len<- length(params$w.vals)
			plot.w.vals<-params$w.vals[1:w.val.len]
			
			# w.val.indices<-c(1:4,7,9,10)
			
			w.val.indices<- 1:w.val.len
			#w.val.indices<-c(1,5,8)
			
			
			plot.w.vals <- plot.w.vals[w.val.indices]
			num.w.vals.plotted <- length(plot.w.vals)
			
			for (i in 1:length(plot.w.vals)){
				print(i)
				w.val.i <- plot.w.vals[i] 
				lty.i <- 1+((i-1)%%10)
				line.width<-1.5
				lty.i.vec <- c(lty.i.vec,lty.i)
				lwd.i.vec <- c(lwd.i.vec,line.width)
				par(lty=lty.i)
				plot.ROC.with.CI(sim.res.g$power[w.val.indices[i],,],plot.title="",plot.col = colors.vec[i],
						conf.int=FALSE,add=(i>1),linewd=line.width,ylim=1)
				
			}
			
			if (compare.pom.cca) {
				line.width<-3
				i<- 1+((num.w.vals.plotted)%%10)
				par(lty=i)
				plot.ROC.with.CI(sim.res.g$power.cmp$pom,plot.title="",plot.col = colors.vec[num.w.vals.plotted+1],
						conf.int=FALSE,add=TRUE,linewd=line.width,ylim=1)
				lty.i.vec <- c(lty.i.vec,i)	
				lwd.i.vec <- c(lwd.i.vec,line.width)		
				i<- 1+((num.w.vals.plotted+1)%%10)
				par(lty=i) 
				plot.ROC.with.CI(sim.res.g$power.cmp$cca,plot.title="",plot.col = colors.vec[num.w.vals.plotted+2],
						conf.int=FALSE,add=TRUE,linewd=line.width,ylim=1)
				lty.i.vec <- c(lty.i.vec,i)	
				lwd.i.vec <- c(lwd.i.vec,line.width)
				
				if (cca.reg){
					i<- 1+((num.w.vals.plotted+2)%%10)
					par(lty=i) 
					plot.ROC.with.CI(sim.res.g$power.cmp$cca.reg,plot.title="",plot.col = colors.vec[num.w.vals.plotted+3],
							conf.int=FALSE,add=TRUE,linewd=line.width,ylim=1)
					lty.i.vec <- c(lty.i.vec,i)	
					lwd.i.vec <- c(lwd.i.vec,line.width)
				}
				
			}
			
			
			i<- 1+((num.w.vals.plotted+3)%%10)
			par(lty=i)
			lines(sim.res.g$optim.power,col=colors.vec[num.w.vals.plotted+4])
			lty.i.vec <- c(lty.i.vec,i)	
			
			
			
			
			legend.txt <- plot.w.vals
			if (compare.pom.cca)
				legend.txt <-c(legend.txt ,"pom","cca")
			
			if (cca.reg)
				legend.txt <-c(legend.txt ,"cca.reg")
			legend.txt <-c(legend.txt ,"bound")
			legend("bottomright",legend=legend.txt,
					col=colors.vec,lty=lty.i.vec,lwd=lwd.i.vec)
			
			fname<- paste(c(results.dir,"MVN","-FC-Tradeoff-",ifelse(oos,"OOS","noOOS"),"c",params$c.val,"-n",params$n),sep="",collapse="")
			
			if(!run.in.linux&(!run.for.Sweave))  savePlot(paste(fname,".pdf",sep="",collapse=""),type="pdf")
			
			if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste(fname,".ps",sep="",collapse=""),type="ps")
						
			if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste(fname,".png",sep="",collapse=""),type="png")
			
			# Zoomed in plot of ROC curves
			windows()
			lty.i.vec<-c()
			lwd.i.vec <- c()
			line.width<-1.5
			w.val.len<- length(params$w.vals)
			plot.w.vals<-params$w.vals[1:w.val.len]
			
			#w.val.indices<-c(1:4,7,9,10)
			w.val.indices<- 1:w.val.len
			#w.val.indices<-c(1,5,8)
			plot.w.vals <- plot.w.vals[w.val.indices]
			
			num.w.vals.plotted <- length(plot.w.vals)
			
			for (i in 1:length(plot.w.vals)){
				print(i)
				w.val.i <- plot.w.vals[i] 
				lty.i <- 1+((i-1)%%10)
				line.width<-1.5
				lty.i.vec <- c(lty.i.vec,lty.i)
				lwd.i.vec <- c(lwd.i.vec,line.width)
				par(lty=lty.i)
				plot.ROC.with.CI(sim.res.g$power[w.val.indices[i],,],plot.title="",plot.col = colors.vec[i],
						conf.int=FALSE,add=(i>1),linewd=line.width,xlim=0.2,ylim=1)
				
			}
			
			if (compare.pom.cca) {
				line.width<-3
				i<- 1+((num.w.vals.plotted)%%10)
				par(lty=i)
				plot.ROC.with.CI(sim.res.g$power.cmp$pom,plot.title="",plot.col = colors.vec[num.w.vals.plotted+1],
						conf.int=FALSE,add=TRUE,linewd=line.width,xlim=0.2,ylim=1)
				lty.i.vec <- c(lty.i.vec,i)	
				lwd.i.vec <- c(lwd.i.vec,line.width)		
				
				
				i<- 1+((num.w.vals.plotted+1)%%10)
				
				par(lty=i) 
				plot.ROC.with.CI(sim.res.g$power.cmp$cca,plot.title="",plot.col = colors.vec[num.w.vals.plotted+2],
						conf.int=FALSE,add=TRUE,linewd=line.width,xlim=0.2,ylim=1)
				lty.i.vec <- c(lty.i.vec,i)	
				lwd.i.vec <- c(lwd.i.vec,line.width)
			}
			
			
			i<- 1+((num.w.vals.plotted+2)%%10)
			par(lty=i)
			lines(sim.res.g$optim.power,col=colors.vec[num.w.vals.plotted+3])
			lty.i.vec <- c(lty.i.vec,i)	
			legend.txt <- plot.w.vals
			if (compare.pom.cca)
				legend.txt <-c(legend.txt ,"pom","cca","bound")
			
			legend("bottomright",legend=legend.txt,
					col=colors.vec[1:(num.w.vals.plotted+2)],lty=lty.i.vec,lwd=lwd.i.vec)
			
			
			fname<- paste(c(results.dir,"MVN","-FC-Tradeoff-",ifelse(oos,"OOS","noOOS"),"c",params$c.val,"-n",params$n,"-zoomed"),sep="",collapse="")
					
			if(!run.in.linux&(!run.for.Sweave)) savePlot(paste(fname,".pdf",sep="",collapse=""),type="pdf")
			
			if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste(fname,".ps",sep="",collapse=""),type="ps")
			
			if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste(fname,".png",sep="",collapse=""),type="png")
			
		}
		
		
		begin.time <-Sys.time()
		
		run.time <-0
		
		if (dirichlet.sim){
		params$p<-5
		params$r<-20
			
			params$q<-22
			params$oos <- TRUE
			
			params$d<-2
			
			
			params.text.1 <- bquote(p==.(params$p) | r==.(params$r) |q==.(params$q) |c==.(params$c.val)|d==.(params$d))
			params.text.3 <- bquote(n==.(params$n) |assume.matched.oos==.(params$assume.matched.for.oos) | nmc==.(params$nmc)| oos==.(params$oos) )
			model.letter <-"D"
			plot.title.D <-paste(model.letter,deparse(params.text.1),"\n,",
					#"\n",		
					deparse(params.text.3) ,"\n")
			
			params$plot.title <-plot.title.D
			#parameter logging start
			sink(paste(results.dir,model.letter,"c",params$c.val,"params",Sys.Date(),".txt"))
			
			print(plot.title.D)
			print(params)
			sink()
			#parameter end
			
			print("Dirichlet Setting Simulation Starting")
			
			sim.res.d<-simulate.generate.test.model.plot("Dirichlet",params,par.compute)
			run.time<-Sys.time()-begin.time
			print("Dirichlet Setting Simulation Ended")
			
			
			Fid1<-sim.res.d$FidComm.Terms$F1
			Fid2<-sim.res.d$FidComm.Terms$F2
			Comm<-sim.res.d$FidComm.Terms$C
			
			Fid1.List <- c(Fid1.List,list(Fid1))
			Fid2.List <-c(Fid2.List,list(Fid2))
			Comm.List <-c(Comm.List,list(Comm))
			
			avgFid1.n<-colMeans(Fid1)
			avgFid2.n<-colMeans(Fid2)
			avgComm.n<-colMeans(Comm)
			
			avg.Fid1.alt<-rbind(avg.Fid1.alt,avgFid1.n)
			avg.Fid2.alt<-rbind(avg.Fid2.alt,avgFid2.n)
			avg.Comm.alt<-rbind(avg.Comm.alt,avgComm.n)
			
			FC.ratio.n<-colMeans(sim.res.d$FC.ratios$f.c)
			wtFC.ratio.n<-colMeans(sim.res.d$FC.ratios$wtf.c)
			
			
			#avg.FCRatios.List<-rbind(FCRatios.List,FC.ratio.n)
			#avg.wtFCRatios.List<-rbind(wtFCRatios.List,wtFC.ratio.n)
			estim.wstar<-c(estim.wstar,sim.res.d$wstar)	
			
			
			
			
			
			
			save.image(paste(results.dir,"Dirichlet","c",params$c.val,"-n",params$n,Sys.Date(),".RData",collapse=""))
			
			
			avg.cont.table<- (sim.res.d$conting.table+0.001)/nmc
			print("Aggregate Cont Table: ")
			print(avg.cont.table)
			p.val <-mcnemar.test(avg.cont.table)$p.value
			print("Aggregate Cont Table: McNemar's p-val")
			print(p.val)
			
			sink(paste(results.dir,model.letter,"c",params$c.val,"-n",params$n,"sign-tests.txt",collapse=""))
			sign.test <- try(sign.test.cont.table(sim.res.d$conting.table.list))
			
			if (inherits(sign.test,"try-error")) {
				print(paste("error in ",model,collapse=""))
				print(sim.res.d$conting.table.list)
			}
			else{
				print("Cont Table List: sign test p-val")
				print(sign.test$p.value)
			}
			
			
			sign.rank.sum.test<-sign.rank.sum.test.cont.table(sim.res.d$conting.table.list)
			print("Cont Table List: signed rank sum test p-val")
			print(sign.rank.sum.test$p.value)
			sink()
			
			sim.res.d$conting.table.list <- lapply(sim.res.d$conting.table.list,function(x) x+0.001)
			
			p.vals.list <- lapply(sim.res.d$conting.table.list,mcnemar.test)
			p.vals<-c()
			for (t in p.vals.list)
				p.vals<- c(p.vals,t$p.value)
			p.vals <- sort(p.vals)
			print("MC Replicate Cont Tables: McNemar's p-val")
			
			
			save.image(paste(results.dir,"Dirichlet","c",params$c.val,"-n",params$n,Sys.Date(),".RData",collapse=""))
			
			lty.i.vec<-c()
			lwd.i.vec <- c()
			line.width<-1.5
			w.val.len<- length(params$w.vals)
			plot.w.vals<-params$w.vals[1:w.val.len]
			
			#w.val.indices<-c(1:4,7,9,10)
			w.val.indices<- 1:w.val.len
			plot.w.vals <- plot.w.vals[w.val.indices]
			num.w.vals.plotted <- length(plot.w.vals)
			windows()
			
			for (i in 1:length(plot.w.vals)){
				print(i)
				w.val.i <- plot.w.vals[i] 
				lty.i <- 1+((i-1)%%10)
				
				line.width<-1.5
				lty.i.vec <- c(lty.i.vec,lty.i)
				lwd.i.vec <- c(lwd.i.vec,line.width)
				par(lty=lty.i)
				plot.ROC.with.CI(sim.res.d$power[w.val.indices[i],,],plot.title="",plot.col = colors.vec[i],
						conf.int=FALSE,add=(i>1),linewd=line.width,ylim=1)
				
			}
			
			if (compare.pom.cca) {
				line.width<-3
				i<- 1+((num.w.vals.plotted)%%10)
				par(lty=i)
				plot.ROC.with.CI(sim.res.d$power.cmp$pom,plot.title="",plot.col = colors.vec[num.w.vals.plotted+1],
						conf.int=FALSE,add=TRUE,linewd=line.width,ylim=1)
				lty.i.vec <- c(lty.i.vec,i)	
				lwd.i.vec <- c(lwd.i.vec,line.width)		
				i<- 1+((num.w.vals.plotted+1)%%10)
				par(lty=i) 
				plot.ROC.with.CI(sim.res.d$power.cmp$cca,plot.title="",plot.col = colors.vec[num.w.vals.plotted+2],
						conf.int=FALSE,add=TRUE,linewd=line.width,ylim=1)
				lty.i.vec <- c(lty.i.vec,i)	
				lwd.i.vec <- c(lwd.i.vec,line.width)
				if (cca.reg) {
					i<- 1+((num.w.vals.plotted+2)%%10)
					par(lty=i) 
					plot.ROC.with.CI(sim.res.d$power.cmp$cca.reg,plot.title="",plot.col = colors.vec[num.w.vals.plotted+3],
							conf.int=FALSE,add=TRUE,linewd=line.width,ylim=1)
					lty.i.vec <- c(lty.i.vec,i)	
					lwd.i.vec <- c(lwd.i.vec,line.width)
				}
				
				
			}
			
			
			
			i<- 1+((num.w.vals.plotted+3)%%10)
			par(lty=i)
			lines(sim.res.d$optim.power,col=colors.vec[num.w.vals.plotted+4])
			lty.i.vec <- c(lty.i.vec,i)	
			legend.txt <- plot.w.vals
			if (compare.pom.cca)
				legend.txt <-c(legend.txt ,"pom","cca")
			if (cca.reg)
				legend.txt <-c(legend.txt ,"cca.reg")
			
			legend.txt <-c(legend.txt ,"bound")
			legend("bottomright",legend=legend.txt,
					col=colors.vec,lty=lty.i.vec,lwd=lwd.i.vec)
			
			#Write plot to file
			
			fname<- paste(c(results.dir,"Dirichlet","-FC-Tradeoff-",ifelse(oos,"OOS","noOOS"),"c",params$c.val,"-n",params$n),sep="",collapse="")
			
			
			if(!run.in.linux&(!run.for.Sweave))  savePlot(paste(fname,".pdf",sep="",collapse=""),type="pdf")
			
			if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste(fname,".ps",sep="",collapse=""),type="ps")
			
			if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste(fname,".png",sep="",collapse=""),type="png")
			
			
			
			
			# Zoomed in plot of ROC curves
			windows()
			lty.i.vec<-c()
			lwd.i.vec <- c()
			line.width<-1.5
			w.val.len<- length(params$w.vals)
			plot.w.vals<-params$w.vals[1:w.val.len]
			
			w.val.indices<- 1:w.val.len
			plot.w.vals <- plot.w.vals[w.val.indices]
			num.w.vals.plotted <- length(plot.w.vals)
			
			for (i in 1:length(plot.w.vals)){
				print(i)
				w.val.i <- plot.w.vals[i] 
				lty.i <- 1+((i-1)%%10)
				
				line.width<-1.5
				lty.i.vec <- c(lty.i.vec,lty.i)
				lwd.i.vec <- c(lwd.i.vec,line.width)
				par(lty=lty.i)
				plot.ROC.with.CI(sim.res.d$power[w.val.indices[i],,],plot.title="",plot.col = colors.vec[i],
						conf.int=FALSE,add=(i>1),linewd=line.width,xlim=0.2,ylim=1)
				
			}
			
			if (compare.pom.cca) {
				line.width<-3
				i<- 1+((num.w.vals.plotted)%%10)
				par(lty=i)
				plot.ROC.with.CI(sim.res.d$power.cmp$pom,plot.title="",plot.col = colors.vec[num.w.vals.plotted+1],
						conf.int=FALSE,add=TRUE,linewd=line.width,xlim=0.2,ylim=1)
				lty.i.vec <- c(lty.i.vec,i)	
				lwd.i.vec <- c(lwd.i.vec,line.width)		
				i<- 1+((num.w.vals.plotted+1)%%10)
				par(lty=i) 
				plot.ROC.with.CI(sim.res.d$power.cmp$cca,plot.title="",plot.col = colors.vec[num.w.vals.plotted+2],
						conf.int=FALSE,add=TRUE,linewd=line.width,xlim=0.2,ylim=1)
				lty.i.vec <- c(lty.i.vec,i)	
				lwd.i.vec <- c(lwd.i.vec,line.width)
				
				
			}
			
			i<- 1+((num.w.vals.plotted+2)%%10)
			par(lty=i)
			lines(sim.res.d$optim.power,col=colors.vec[num.w.vals.plotted+3])
			lty.i.vec <- c(lty.i.vec,i)	
			
			
			
			legend.txt <- plot.w.vals
			if (compare.pom.cca)
				legend.txt <-c(legend.txt ,"pom","cca")
			
			
			
			legend.txt <-c(legend.txt ,"bound")
			legend("bottomright",legend=legend.txt,
					col=colors.vec,lty=lty.i.vec,lwd=lwd.i.vec)
			
			fname<- paste(c(results.dir,"Dirichlet","-FC-Tradeoff-",ifelse(oos,"OOS","noOOS"),"c",params$c.val,"-n",params$n,"-zoomed"),sep="",collapse="")
			
			if(!run.in.linux&(!run.for.Sweave))  savePlot(paste(fname,".pdf",sep="",collapse=""),type="pdf")
			
			if(!run.in.linux&(!run.for.Sweave))  savePlot(filename=paste(fname,".ps",sep="",collapse=""),type="ps")
			
			if(!run.in.linux&(!run.for.Sweave)) savePlot(filename=paste(fname,".png",sep="",collapse=""),type="png")
			
			
			par(lty=1)
			
		}
	} #end c.vals
	
	
	
	avg.Fid1.alt<-rbind(avg.Fid1.alt,avgFid1.n)
	avg.Fid2.alt<-rbind(avg.Fid2.alt,avgFid2.n)
	avg.Comm.alt<-rbind(avg.Comm.alt,avgComm.n)
	
} #end n.vals

Rprof(NULL)



FC.results.mat <-list(n.vals,F1=avg.Fid1,F2=avg.Fid2,C=avg.Comm,
		avg.Fid1.alt,avg.Fid2.alt,avg.Comm.alt,FC.rat=avg.FCratios.List,wtFC.rat=avg.wtFCratios.List)
sink("FC.results.txt")
print(FC.results.mat)
sink()


sink("FC.results-debug.txt")
Fid1.List
Fid2.List
Comm.List
sink()



sink("wstar-vals")

print(estim.wstar)
sink()


sink("traceback.txt")
traceback()
sink()
warnings()

#Rprof(NULL)
if ((!run.in.linux) & par.compute & (!par.compute.sf)) 
	stopWorkers(workers)
