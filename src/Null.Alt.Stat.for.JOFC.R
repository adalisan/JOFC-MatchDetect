# TODO: Add comment
# 
# Author: Sancar
###############################################################################

results.dir <- "graphs"

par.compute <- TRUE
run.in.linux<- FALSE
compare.pom.cca<-FALSE


verbose<-FALSE
oos <-TRUE

gauss.sim <-T
dirichlet.sim <- F
run.mcnemars.test <- F
cca.reg <- F
power.comparison.test<- F








#n.vals <-c(50,100,150,200,300,500)
#
n.vals<-c(150)
nmc <-  200
s<- 1

debug.mode<-FALSE
if(debug.mode){
	n.vals<-c(50)
	nmc<-3
	s<-1
}





params<-list(
		nmc.param=nmc,
		p.param = 5,
		r.param = 10,
		
		q.param =5,
		
		d.param=3,
		n.param=n.vals[1],
		s.param=s,
		
		c.val.param=0.1,
		#w.vals.param = c(0.001,0.1,0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999),
		#w.vals.param = c(0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999),
		#w.vals.param = c(0.5,0.75,0.95,0.999),
		w.vals.param = 0.75,
		
		
		
		Wchoice="NA+diag(0)",
		pre.scaling=TRUE,
		oos =oos,
#
# assume oos observations in different conditions are matched for oos-embedding
#
		assume.matched.for.oos = FALSE,
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
		wt.equalize=FALSE,
		rival.w = 0.999,
		power.comparison.test = F
)


w.max.index <- length(params$w.vals.param)	
size <- seq(0, 1, 0.01)
len <- length(size)



model="gaussian"

pre.scaling = TRUE  #Make the measurement spaces have the same scale

alpha       = NULL

sim.grass = FALSE
eigen.spectrum = FALSE



level.mcnemar =0.01  #At what alpha, should unweighted(w=0.5) and optimal w^* be compared
def.w = 0.5          #The null hypothesis is that power(def.w) >= power(rival.w) (by default ,def.w is the w for the unweighted case which equals 0.5)

proc.dilation=FALSE #when investigating convergence of JOFC to PoM, should Procrustes analysis of configurations include the dilation component?

old.gauss.model.param <- FALSE

param.index.count<-20
params.list <- rep(list(params),param.index.count)
p.vals<-rep(c(5,12,19,25,35),each=4)
r.vals<-rep(c(5,20,50,150),5)

for (i in (1:length(p.vals))){
	params.list[[i]]$p.param <- p.vals[i]
#params.list[[6]]$p <- 25
	
}

for (i in (1:length(r.vals))){
	params.list[[i]]$r.param <- r.vals[i]
	params.list[[i]]$p.param <- p.vals[i]
}

#params.list[[5]]$r <- 100

#q.vals<-rep(c(5,10,30,40),5)
q.vals<-rep(c(5,10,30,40))
vary.param<-"p_r"
#params.list[[1]]$q.param <- 22


#Function runs JOFC embedding with the given parameters to get test statistics 
run.one.repl.with.params<-function(mc.counter,params){
	require(MASS)
	results<-with(params, {
				require(MASS)
				w.max.index <-length(w.vals.param)
				oos.diss.mc.0<- array(0,dim=c(w.max.index ))
				oos.diss.mc.A<- array(0,dim=c(w.max.index ))
				power.w.star <- 0
				
				m<-s.param
				sink(paste(mc.counter,"c-check.txt"))
				print(c.val.param)
				sink()
				T0 <- matrix(0,w.max.index,m)   #Test statistics for JOFC under null
				TA <- matrix(0,w.max.index,m)    #Test statistics for JOFC under alternative
				Y.oos.matched <- c()
				Y.oos.unmatched<-c()
				
				
				cont.table  <- matrix(0,2,2)
				
				Fid.Err.Term.1 <- array(0,dim=c(w.max.index))
				Fid.Err.Term.2 <- array(0,dim=c(w.max.index))
				Comm.Err.Term <- array(0,dim=c(w.max.index))
				
				sigma <- matrix(0,p.param,p.param)
				means <- array(0 , dim=c(w.max.index,2*d.param))
				
				
				
				
					if (model=="gaussian"){
						sigma<- diag(p.param)
						if (old.gauss.model.param) sigma <-Posdef(p.param,r.param)
						alpha.mc <- mvrnorm(n.param+(2*m), rep(0,p.param),sigma)
					} else if (model=="dirichlet"){
						alpha.mc <- rdirichlet(n.param+(2*m), rep(1,p.param+1))
					} else stop("unknown model")
					
					
				
				
				## n.param pairs of matched points
				if (model=="gaussian"){
					xlist <- matched_rnorm(n.param, p.param, q.param, c.val.param, r.param, alpha=alpha.mc[1:n.param, ],sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)
				} else{
					xlist <- matched_rdirichlet(n.param, p.param, r.param, q.param, c.val.param, alpha.mc[1:n.param, ])
				}
				
				## m pairs of matched points
				if (model=="gaussian"){
					ylist <- matched_rnorm(m, p.param, q.param, c.val.param, r.param, alpha=as.array(alpha.mc[n.param+(1:m), ,drop=FALSE]),sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)
				} else{
					ylist <- matched_rdirichlet(m, p.param, r.param, q.param, c.val.param, alpha=as.array(alpha.mc[n.param+(1:m), ,drop=FALSE]))
				}
				if (model=="gaussian"){
					y.alt <- matched_rnorm(m, p.param, q.param, c.val.param, r.param, alpha=as.array((alpha.mc[n.param+m+(1:m), ,drop=FALSE])),sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)
				} else{
					y.alt <- matched_rdirichlet(m, p.param, r.param, q.param, c.val.param,alpha= as.array(alpha.mc[n.param+m+(1:m), ,drop=FALSE]))
				}
				
				
				X1 <- xlist$X1
				X2 <- xlist$X2
				Y1<- ylist$X1
				Y20<- ylist$X2
				Y2A<- y.alt$X2
				
				if (model=="gaussian")
					sigma.mc<-xlist$sigma.beta
				
				D1 <- dist(X1)
				D2 <- dist(X2)
				
				D.whole <- dist(rbind(X1,X2))
				
				
				if (verbose) print("random matched pairs generated\n")
				
#prescaling
				if (pre.scaling) {
					s <- lm(as.vector(D1) ~ as.vector(D2) + 0)$coefficients
				} else {
					s <- 1
				}
				
				
				D1<-as.matrix(D1)
				D2<-as.matrix(D2)
				pom.config<-c()
				cca.config<-c()
				D2<-D2*s
				
				## ==== jofc ====
				
# Impute "between-condition" dissimilarities from different objects  
				if (Wchoice == "avg") {
					L <- (D1 + D2)/2
				} else if (Wchoice == "sqrt") {
					L <- sqrt((D1^2 + D2^2)/2)
				} else if (Wchoice == "NA+diag(0)") {
					L <- matrix(NA,n.param,n.param)
					diag(L)<- 0
				}
				
				
				#In sample embedding
				# Form omnibus dissimilarity matrix
				M <- omnibusM(D1, D2, L)
				init.conf<-NULL
				
				if (compare.pom.cca) init.conf<- pom.config
				
				# Embed in-sample using different weight matrices (differentw values)
				X.embeds<-JOFC.Insample.Embed(M, d.param, w.vals.param, separability.entries.w, init.conf=init.conf,
						wt.equalize=wt.equalize)
				
				Pert.mat.mc<-(as.matrix(dist(X.embeds[[1]])^2))-as.matrix(D.whole^2)
				
				#Form OOS omnibus matrices
				M.oos.0<- matrix(0,2,2)#omnibusM(D.oos.1,D.oos.2.null, L.tilde.null)
				M.oos.A<- matrix(0,2,2)#omnibusM(D.oos.1,D.oos.2.alt,  L.tilde.alt)
				
				
				
				for (l in 1:w.max.index){
					if (verbose) print("OOS embedding for JOFC for w= \n")
					if (verbose) print(w.vals.param[l])
					
					w.val.l <- w.vals.param[l]
					X <- X.embeds[[l]]
					
					
					oos.obs.flag<- c(rep(1,2*n.param),rep(0,2*m))
					
					#Compute Weight matrix corresponding in-sample  entries
					oos.Weight.mat.1<-w.val.to.W.mat(w.val.l,(2*n.param),separability.entries.w,wt.equalize)
					
					#Compute Weight matrix corresponding OOS  entries
					#oos.Weight.mat.2<-w.val.to.W.mat(w.val.l,(2*m),separability.entries.w,wt.equalize)
					oos.Weight.mat.2<-matrix(c(0,w.val.l,w.val.l,0),2,2)
					
					# If assume.matched.for.oos is true, we assume OOS dissimilarities are matched(in reality,
					# they are matched for the matched pairs, but unmatched for the unmatched pairs)
					# If assume.matched.for.oos is true, we ignore the dissimilarities between matched/unmatched 
					# pairs
					if (!assume.matched.for.oos){
						oos.Weight.mat.2[1:m,m+(1:m)]<-0
						oos.Weight.mat.2[m+(1:m),(1:m)]<-0
					}
					# if (oos.use.imputed is true) we treat the dissimiilarities between  in-sample and out-of-sample measurements
					# from different conditions like fidelity terms
					# otherwise they are ignored
					if (oos.use.imputed){
						oos.Weight.mat.w <- matrix(1-w.val.l,2*n.param,2*m)
					} else{
						oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n.param,m), matrix(0,n.param,m) ),
								cbind(matrix(0,n.param,m),matrix(1-w.val.l,n.param,m))
						)
					}
					print(dim(oos.Weight.mat.1))
					print(dim(oos.Weight.mat.2))
					print(dim(oos.Weight.mat.w))
					oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
					# Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
					# We are using previous in-sample embeddings, anyway
					oos.Weight.mat[1:(2*n.param),1:(2*n.param)]<-0
					if (verbose) print("dim(M.oos.0)")
					if (verbose) print(dim(M.oos.0))
					if (verbose) print("dim(M.oos.A)")
					if (verbose) print(dim(M.oos.A))
					if (verbose) print("dim(oos.Weight.mat)")
					if (verbose) print(dim(oos.Weight.mat))
					if (verbose) print("dim(X)")
					if (verbose) print(dim(X))
					#if (verbose) {print("oos.obs.flag")}
					
					
					ideal.omnibus.0  <- as.matrix(dist(rbind(X1,X2,Y1,Y20)))
					ideal.omnibus.A  <- as.matrix(dist(rbind(X1,X2,Y1,Y2A)))
					omnibus.oos.D.0 <- omnibusM(M,M.oos.0,ideal.omnibus.0[1:(2*n.param),(2*n.param)+(1:(2*m))])
					omnibus.oos.D.A <- omnibusM(M,M.oos.A, ideal.omnibus.A[1:(2*n.param),(2*n.param)+(1:(2*m))])
					oos.Weight.mat [is.na(omnibus.oos.D.0)]<-0
					omnibus.oos.D.0[is.na(omnibus.oos.D.0)]<-1
					omnibus.oos.D.A[is.na(omnibus.oos.D.A)]<-1
					if (verbose) print("JOFC null omnibus OOS embedding \n")
#if (profile.mode)			Rprof("profile-oosIM.out",append=TRUE)
					Y.0t<-oosIM(D=omnibus.oos.D.0,
							X=X,
							init     = "random",
							verbose  = FALSE,
							itmax    = 1000,
							eps      = 1e-8,
							W        = oos.Weight.mat,
							isWithin = oos.obs.flag,
							bwOos    = FALSE)
					
					
					embedded.dist.0<- sum((Y.0t[1,]-Y.0t[2,])^2)
					if (verbose) print("JOFC alternative omnibus OOS embedding \n")
					Y.At<-oosIM(D=omnibus.oos.D.A,
							X=X,
							init     = "random",
							verbose  = FALSE,
							itmax    = 1000,
							eps      = 1e-8,
							W        = oos.Weight.mat,
							isWithin = oos.obs.flag,
							bwOos    = FALSE)
#if (profile.mode)				Rprof(NULL)
					embedded.dist.A<- sum((Y.At[1,]-Y.At[2,])^2)	
					
					oos.diss.mc.0[l]<- embedded.dist.0
					oos.diss.mc.A[l]<- embedded.dist.A
				}
				list(oos.0 = oos.diss.mc.0,oos.A = oos.diss.mc.A
				#p1=Pert.mat.mc,p2=as.matrix(D.whole)
				)
			}	
	)
	
}

compute.stats.for.params<-function(params.list,vary.param){
	
	require(MASS)
	T0.agg<-array(0,dim=c(nmc,w.max.index,params.list[[1]]$s))
	TA.agg<-array(0,dim=c(nmc,w.max.index,params.list[[1]]$s))
	
	
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
	
	means.0<-c()
	vars.0<-c()
	skews.0<-c()
	kurts.0 <- c()
	
	means.A<-c()
	vars.A<-c()
	skews.A<-c()
	kurts.A <- c()
	
	means.0.agg<-c()
	means.0.agg<-c()
	means.A.agg<-c()
	vars.A.agg<-c()

	
	Y.comm.dist.par.agg<- c()
	
	Y.comm.dist.par.g.agg<-c()
	Y.comm.dist.par.log.n.agg<-c()
	Y.sep.dist.par.agg<- c()
	
#laplace(llocation = "identity", lscale = "loge", elocation = list(),
#		es<cale = list(), ilocation = NULL, iscale = NULL,
#		imethod = 1, zero = 2)
#fit = vglm(y ~ 1, laplace, lddat, trace = TRUE, crit = "l")
	
	
	
	r.comm.dist.par.agg<-c()
	
	r.comm.dist.par.log.n.agg<-c()
	r.sep.dist.par.agg<- c()
	
	
	
	oos.diss.0.exp.dist.par.agg<-c()
	
	
	oos.diss.A.exp.dist.par.agg<- c()
	oos.diss.0.par.log.n.agg<- c()
	oos.diss.A.par.log.n.agg<-c()
	
	
	
	
	oos.diss.all.params.0 <-array(0,dim=c(param.index.count,w.max.index,nmc))
	oos.diss.all.params.A <-array(0,dim=c(param.index.count,w.max.index,nmc))
#param.index<-1
    

    
	for (param.index in 1:param.index.count){
		
		
		params<- params.list[[param.index]]
		
		
		
		#sink(paste("log-params-",vary.param,param.index,".txt",sep="",collapse=""))

		
		sfStop()
		sfInit(parallel=TRUE,cpus=num.cpus)
		sfExport("w.max.index",
				"size",
				"len",
				"model",
				
				
				"pre.scaling", 
	
#		alpha       = NULL
				"sim.grass",
				"eigen.spectrum",
#		
#		
#		
				"level.mcnemar",
				
				"def.w",          #The null hypothesis is that power(def.w) >= power(rival.w) (by default ,def.w is the w for the unweighted case which equals 0.5)
#		
				"proc.dilation", #when investigating convergence of JOFC to PoM, should Procrustes analysis of configurations include the dilation component?
#		
				"old.gauss.model.param",
				"compare.pom.cca"
						
		)
		sfLibrary(MASS)
		sfSource("./lib/oosIM.R")
		sfSource("./lib/simulation_math_util_fn.R")
		sfSource("./lib/smacofM.R")
#		attach(params)		
#		pprime1     = ifelse(model=="gaussian",p+q,p+q+2)   # cca arguments , signal+noise dimension
#		pprime2     = ifelse(model=="gaussian",p+q,p+q+2)   # cca arguments, signal+noise dimension
#		
#		pre.scaling = TRUE  #Make the measurement spaces have the same scale
#		
#		alpha       = NULL
#		
#		sim.grass=FALSE
#		eigen.spectrum=FALSE
#		
#		
#		
#		level.mcnemar=0.01  #At what alpha, should unweighted(w=0.5) and optimal w^* be compared
#		def.w=0.5          #The null hypothesis is that power(def.w) >= power(rival.w) (by default ,def.w is the w for the unweighted case which equals 0.5)
#		
#		proc.dilation=FALSE #when investigating convergence of JOFC to PoM, should Procrustes analysis of configurations include the dilation component?
#		
#		old.gauss.model.param <- FALSE
#		detach(params)
#		
		
		
		
		results.agg<-sfLapply(1:nmc,run.one.repl.with.params,params)   
		
		for (mc in 1:nmc){
			if (length(params$w.vals.param)==1) {
				oos.diss.all.params.0[param.index,1,mc]<-results.agg[[mc]]$oos.0
				oos.diss.all.params.A[param.index,1,mc]<-results.agg[[mc]]$oos.A
				
			} else{			
			oos.diss.all.params.0[param.index,,mc]<-results.agg[[mc]]$oos.0
			oos.diss.all.params.A[param.index,,mc]<-results.agg[[mc]]$oos.A
			}
		} 
		sfStop()
		
		
	}
	
	return(list(stats.0=oos.diss.all.params.0,stats.A=oos.diss.all.params.A))
}