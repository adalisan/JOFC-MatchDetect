## Time-stamp: <wiki_classification_SA&DM.R zma 2010-09-23 00:06>








#
#
#run.jofc.graph <- function(D1, D2, D10A,D20,D2A,
#					D.oos.1,
#					D.oos.2.null ,
#					D.oos.2.alt ,
#					
#					ideal.omnibus.0  ,
#					ideal.omnibus.A ,
#	
#				n,m,
#				d,c.val,
#				model,oos,Wchoice,separability.entries.w,wt.equalize,assume.matched.for.oos,oos.use.imputed,
#				w.vals,
#				verbose=FALSE)   {
#	
#
#	
#	w.max.index <- length(w.vals)
#	T0 <- matrix (0,w.max.index, m)
#	TA <- matrix (0,w.max.index, m)
#	## ==== jofc ====
#	
#	# Impute "between-condition" dissimilarities from different objects  
#	if (Wchoice == "avg") {
#		L <- (D1 + D2)/2
#	} else if (Wchoice == "sqrt") {
#		L <- sqrt((D1^2 + D2^2)/2)
#	} else if (Wchoice == "NA+diag(0)") {
#		L <- matrix(NA,n,n)
#		diag(L)<- 0
#	}
#	
#	
#	if (oos == TRUE) {
#		
#		#In sample embedding
#		# Form omnibus dissimilarity matrix
#		M <- omnibusM(D1, D2, L)
#		init.conf<-NULL
#		
#		
#		# Embed in-sample using different weight matrices (differentw values)
#		X.embeds<-JOFC.Insample.Embed(M,d,w.vals,separability.entries.w,init.conf=init.conf,wt.equalize=wt.equalize)
#		if (0) {
#		Fid.Err.Term.1 <- X.embeds[[w.max.index+2]]
#		Fid.Err.Term.2 <- X.embeds[[w.max.index+3]]
#		Comm.Err.Term  <- X.embeds[[w.max.index+4]]
#		
#		Fid.Err.Sum.Term.1 <- X.embeds[[w.max.index+5]]
#		Fid.Err.Sum.Term.2 <- X.embeds[[w.max.index+6]]
#		Comm.Err.Sum.Term  <- X.embeds[[w.max.index+7]]
#		FC.ratio  <- X.embeds[[w.max.index+8]]
#		FC.ratio.2  <- X.embeds[[w.max.index+9]]
#		FC.ratio.3  <- X.embeds[[w.max.index+10]]
#		
#		print("Fid.Err.Term.1" )
#		print(Fid.Err.Term.1 )
#		print("Comm.Err.Term ")
#		print(Comm.Err.Term )
#		min.stress.for.w.val <- X.embeds[[w.max.index+1]]
#		}
#		if (verbose) print("JOFC embeddings complete\n")
#		
#		
#		#
#		# OOS Dissimilarity matrices
#		#
#		
#		
#		
#	
#			
#			
#			
#		#Imputing dissimilarity  entries for OOS
#		if (Wchoice == "avg") {
#			L.tilde.null <- (D.oos.1 + D.oos.2.null)/2
#			L.tilde.alt  <- (D.oos.1 + D.oos.2.alt)/2
#		} else if (Wchoice == "sqrt") {
#			L.tilde.null <- sqrt((D.oos.1^2 + D.oos.2.null^2)/2)
#			L.tilde.alt  <- sqrt((D.oos.1^2 + D.oos.2.alt^2)/2)
#			
#		} else if (Wchoice == "NA+diag(0)") {
#			L.tilde.null <- matrix(NA,m,m)
#			L.tilde.alt <- matrix(NA,m,m)
#			diag(L.tilde.null)<- 0
#			diag(L.tilde.alt)<- 0
#		}
#		#Form OOS omnibus matrices
#		M.oos.0 <- omnibusM(D.oos.1,D.oos.2.null, L.tilde.null)
#		M.oos.A <- omnibusM(D.oos.1,D.oos.2.alt,  L.tilde.alt)
#		
#		
#		
#		for (l in 1:w.max.index){
#			if (verbose) print("OOS embedding for JOFC for w= \n")
#			if (verbose) print(w.vals[l])
#			
#			w.val.l <- w.vals[l]
#			X <- X.embeds[[l]]
#			
#			
#			oos.obs.flag<- c(rep(1,2*n),rep(0,2*m))
#			
#			#Compute Weight matrix corresponding in-sample  entries
#			oos.Weight.mat.1<-w.val.to.W.mat(w.val.l,(2*n),separability.entries.w,wt.equalize)
#			
#			#Compute Weight matrix corresponding OOS  entries
#			oos.Weight.mat.2<-w.val.to.W.mat(w.val.l,(2*m),separability.entries.w,wt.equalize)
#			
#			# If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
#			# they are matched for the matched pairs, but unmatched for the unmatched pairs)
#			# If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched 
#			# pairs
#			if (!assume.matched.for.oos){
#				oos.Weight.mat.2[1:m,m+(1:m)]<-0
#				oos.Weight.mat.2[m+(1:m),(1:m)]<-0
#			}
#			# if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
#			# from different conditions like fidelity terms
#			# otherwise they are ignored
#			if (oos.use.imputed){
#				oos.Weight.mat.w <- matrix(1-w.val.l,2*n,2*m)
#			} else{
#				oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n,m), matrix(0,n,m) ),
#						cbind(matrix(0,n,m),matrix(1-w.val.l,n,m))
#				)
#			}
#			oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
#			# Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
#			# We are using previous in-sample embeddings, anyway
#			oos.Weight.mat[1:(2*n),1:(2*n)]<-0
#			if (verbose) print("dim(M.oos.0)")
#			if (verbose) print(dim(M.oos.0))
#			if (verbose) print("dim(M.oos.A)")
#			if (verbose) print(dim(M.oos.A))
#			if (verbose) print("dim(oos.Weight.mat)")
#			if (verbose) print(dim(oos.Weight.mat))
#			if (verbose) print("dim(X)")
#			if (verbose) print(dim(X))
#			#if (verbose) {print("oos.obs.flag")
#			
#			
#			
#			omnibus.oos.D.0 <- omnibusM(M,M.oos.0, ideal.omnibus.0[1:(2*n),(2*n)+(1:(2*m))])
#			omnibus.oos.D.A <- omnibusM(M,M.oos.A, ideal.omnibus.A[1:(2*n),(2*n)+(1:(2*m))])
#			oos.Weight.mat[is.na(omnibus.oos.D.0)]<-0
#			omnibus.oos.D.0[is.na(omnibus.oos.D.0)]<-1
#			omnibus.oos.D.A[is.na(omnibus.oos.D.A)]<-1
#			if (verbose) print("JOFC null omnibus OOS embedding \n")
##if (profile.mode)			Rprof("profile-oosIM.out",append=TRUE)
#			Y.0t<-oosIM(D=omnibus.oos.D.0,
#					X=X,
#					init     = "random",
#					verbose  = FALSE,
#					itmax    = 1000,
#					eps      = 1e-8,
#					W        = oos.Weight.mat,
#					isWithin = oos.obs.flag,
#					bwOos    = TRUE)
#			
#			if (verbose) print("JOFC alternative omnibus OOS embedding \n")
#			Y.At<-oosIM(D=omnibus.oos.D.A,
#					X=X,
#					init     = "random",
#					verbose  = FALSE,
#					itmax    = 1000,
#					eps      = 1e-8,
#					W        = oos.Weight.mat,
#					isWithin = oos.obs.flag,
#					bwOos    = TRUE)
##if (profile.mode)				Rprof(NULL)
#			Y1t<-Y.0t[1:m,]
#			Y2t<-Y.0t[m+(1:m),]
#			Y1t.A<-Y.At[1:m,]
#			Y2At<-Y.At[m+(1:m),]
#		
#			
#			T0[l,] <- rowSums((Y1t - Y2t)^2)
#			TA[l,] <- rowSums((Y1t.A - Y2At)^2)
#			}
#			
#	
#	
#}
#
#
#return(list(T0=T0,TA=TA))
#
#}

source("./src/wiki_hypTest_2cond_Params.R")

run.wiki.JOFC.sim.mc.replicate <- function(m.i,N, test.samp.size, w.val.len, m, Diss.E, Diss.F
                                           , n, d, model, oos, Wchoice, separability.entries.w
                                           , wt.equalize, assume.matched.for.oos
                                           , oos.use.imputed, w.vals, size, verbose
                                           , level.mcnemar) {
  source("./lib/simulation_math_util_fn.R")
  source("./lib/smacofM.R")
  source("./lib/oosIM.R")
  
  
  power.mc <-array(0,dim=c(w.val.len,length(size)))
  left.out.samp<- sample(1:N,2*test.samp.size)
  test.matched<- left.out.samp[1:test.samp.size] 
  test.unmatched<- left.out.samp[test.samp.size+(1:test.samp.size)] 
  test.matched   <- sort(test.matched)
  test.unmatched <- sort(test.unmatched)
  
  #sample.for.unmatched<- sample(1:n,test.samp.size)
  orig.indices<-1:N
  #sample.for.unmatched.orig.index<-orig.indices[-left.out.samp][sample.for.unmatched]
  
  T0 <- matrix(0,w.val.len,m)   #Test statistics for JOFC under null
  TA <- matrix(0,w.val.len,m)    #Test statistics for JOFC under alternative
  
  D1<-Diss.E[-left.out.samp,-left.out.samp]
  D2<-Diss.F[-left.out.samp,-left.out.samp]
  train.test.0<-orig.indices[-left.out.samp]
  train.test.0<-c(train.test.0,test.matched)
  train.test.A<-orig.indices[-left.out.samp]
  train.test.A<-c(train.test.A,test.unmatched)
  
  D10A<- Diss.E[train.test.0]
  D20<-  Diss.F[train.test.0]
  
  D2A<-  Diss.F[train.test.A]
  
  D.oos.1      <- Diss.E[test.matched,test.matched]
  D.oos.2.null <- Diss.F[test.matched,test.matched]
  D.oos.2.alt  <- Diss.F[test.unmatched,test.unmatched]
  
  
  L.in.oos.0 <- omnibusM.inoos(Diss.E[-left.out.samp,][,test.matched],Diss.F[-left.out.samp,][,test.matched],matrix(0,n,m))
  L.in.oos.A <- omnibusM.inoos(Diss.E[-left.out.samp,][,test.matched],Diss.F[-left.out.samp,][,test.unmatched],matrix(0,n,m))
  
  ideal.omnibus.0 <- omnibusM(omnibusM (D1,D2,matrix(0,n,n)),omnibusM(D.oos.1,D.oos.2.null,matrix(0,m,m)),L.in.oos.0)
  ideal.omnibus.A <- omnibusM(omnibusM (D1,D2,matrix(0,n,n)),omnibusM(D.oos.1,D.oos.2.alt,matrix(0,m,m)),L.in.oos.A)
  
  
  print(str(T0))
  print(str(TA))
  
  print(str(D1))
  print(str(D2))
  print(str(D10A))
  print(str(D20))
  
  
  power.w.star<- 0
  print("starting JOFC embedding ")
  
  
  
  JOFC.results <- run.jofc(	D1, D2, D10A,D20,D2A,
                            D.oos.1, D.oos.2.null ,		D.oos.2.alt ,
                            
                            L.in.oos.0 ,	L.in.oos.A, 	n,m,	d,
                            model,oos,Wchoice,separability.entries.w,wt.equalize,
                            assume.matched.for.oos,oos.use.imputed,
                            pom.config=NULL,
                            w.vals, 	size,	verbose=verbose) 
  
  
  if (verbose) print("JOFC test statistic complete \n")
  
  T0 <- JOFC.results$T0 
  TA <- JOFC.results$TA 
  for (l in 1:w.val.len){
    
    w.val.l <- w.vals[l]					
    
    power.l <- get_power(T0[l,],TA[l,],size)
    power.mc[l,]<-power.l
    power.mcnemar.l <- get_power(T0[l,],TA[l,],level.mcnemar)
    if (power.mcnemar.l>power.w.star){
      rival.w <- w.vals[l]
      power.w.star <- power.mcnemar.l
      w.val.rival.idx <- l
    }
  }
  return(list(T0=T0,TA=TA,power.mc=power.mc))
}



require(parallel)
require(foreach)

num.cores<-parallel::detectCores()
iter_per_core <- ceiling(nmc/num.cores)

cl <- NULL
use.snow<-TRUE




if(!par.compute){
  registerDoSEQ()
} else if ( !use.snow && .Platform$OS.type != "windows" && require("multicore") ) {
  require(doMC)
  registerDoMC()
  
} else if (use.snow && require("doSNOW")) {
  cl <- snow::makeCluster(num.cores,type= "SOCK")
  registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl), add = TRUE)

} else if (                     # doSMP is buggy
  require("doSMP")&& !use.snow) {
  workers <- startWorkers(num.cores,FORCE=TRUE) # My computer has 4 cores

  registerDoSMP(workers)
  on.exit(stopWorkers(workers), add = TRUE)

} else {
  registerDoSEQ()
}  



if (par.compute){
  
  JOFC.wiki.res<-parLapply(cl=cl,1:nmc, run.wiki.JOFC.sim.mc.replicate, N = N, test.samp.size = test.samp.size,
                           w.val.len = w.val.len, m = m, Diss.E = GE, Diss.F = GF, n = n, d=d,
                           model = "gaussian", oos = oos, Wchoice = Wchoice,
                           separability.entries.w = separability.entries.w, wt.equalize = wt.equalize,
                           assume.matched.for.oos = assume.matched.for.oos, oos.use.imputed = oos.use.imputed,
                           w.vals = w.vals, size = size, verbose = verbose,  level.mcnemar = level.mcnemar			
  )
  
  
}  else {	
  JOFC.wiki.res<-lapply(1:nmc, run.wiki.JOFC.sim.mc.replicate, N = N, test.samp.size = test.samp.size,
                        w.val.len = w.val.len, m = m, Diss.E = Diss.E, Diss.F = Diss.F, n = n, d=d,
                        model = "gaussian", oos = oos, Wchoice = Wchoice,
                        separability.entries.w = separability.entries.w, wt.equalize = wt.equalize,
                        assume.matched.for.oos = assume.matched.for.oos, oos.use.imputed = oos.use.imputed,
                        w.vals = w.vals, size = size, verbose = verbose,  level.mcnemar = level.mcnemar			
  )
}




for (mc.i in 1:nmc){
  power.nmc[,mc.i,]<- JOFC.wiki.res[[mc.i]]$power.mc
  
}





#   
#   
# 	JOFC.wiki.res <- foreach(i=1:nmc, .combine="c",.export=c("bitflip_MC_rep","run.experiment.JOFC")) %dopar% {
# 	  
# 	  
# 	  
# 	  
# 	  
# 	  
# 	  mc.rep.result<-run.wiki.JOFC.sim.mc.replicate(m.i=i, N = N, test.samp.size = test.samp.size,
# 	                                                w.val.len = w.val.len, m = m, TE = TE, TF = TF, n = n,
# 	                                                model = "gaussian", oos = oos, Wchoice = Wchoice,
# 	                                                separability.entries.w = separability.entries.w, wt.equalize = wt.equalize,
# 	                                                assume.matched.for.oos = assume.matched.for.oos, oos.use.imputed = oos.use.imputed,
# 	                                                w.vals = w.vals, size = size, verbose = verbose,  level.mcnemar = level.mcnemar			
# 	  )
# 	  list(mc.rep.result)
# 	}	
# 	
# 	for (mc.i in 1:nmc){
# 	  power.nmc[,mc.i,]<- JOFC.wiki.res[[mc.i]]$power.mc
# 	}   

save.image(file= paste("JOFC_Wiki_Exp_HypTest",format(Sys.time(), "%b %d %H:%M:%S"),".RData"))



colors.vec <- c("red","green","gold4","purple","aquamarine",
                "darkblue","azure3","salmon","rosybrown","magenta","orange",
                "darkorange4")


colors.vec.len<-length(colors.vec)
#	colors.vec[colors.vec.len+1]<-"cornflowerblue"
#	colors.vec[colors.vec.len+2]<-"azure3"
colors.vec.len<-length(colors.vec)
par(lty=1)

lty.i.vec<-c()
for (i in 1:w.val.len){
  lty.i <- 1+((i-1)%%10)
  
  lty.i.vec <- c(lty.i.vec,lty.i)
  #par(lty=lty.i)
  plot.ROC.with.CI(power.nmc[i,,],plot.title="",plot.col = colors.vec[i],
                   conf.int=FALSE,add=(i>1),ylim=1)
  
}

dev.copy2pdf(file= paste("JOFC_Wiki_Exp_HypTest",format(Sys.time(), "%b %d %H:%M:%S"),".pdf"))
dev.copy(device=png,file= paste("JOFC_Wiki_Exp_HypTest",format(Sys.time(), "%b %d %H:%M:%S"),".png"))

legend.txt <- w.vals

legend("bottomright",legend=legend.txt,
       col=colors.vec,lty=1)
#							plot.title<-paste("TA",vary.param,unlist(params.list[[param.index]][vary.param]))
#							title(plot.title)




