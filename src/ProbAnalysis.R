# TODO: Add comment
# 
# Author: Sancar
###############################################################################





results.dir <- "graphs"





par.compute <- FALSE
run.in.linux<- FALSE
compare.pom.cca<-FALSE
#c.vals<-c(0.01)
c.vals<-c(0.1)

verbose<-FALSE
oos <-TRUE

gauss.sim <-T
dirichlet.sim <- F
run.mcnemars.test <- F
cca.reg <- F
power.comparison.test<- F
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








#n.vals <-c(50,100,150,200,300,500)
#
n.vals<-c(150)
nmc <-  200
s<- 1

debug.mode<-FALSE
if(debug.mode){
	n.vals<-c(50)
	nmc<-3
	s<-30
}





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


if (gauss.sim)
	sim.res.g<-list()
if (dirichlet.sim)
	sim.res.d<-list()

params<-list(
		nmc=nmc,
		coincid.vec.dotpr.thres =0.9,
		eigen.spectrum =F,
		grassmannian.dist = F,
		use.Euc.points.x = F,
		use.Euc.points.y = F,
		p = 5,
		r = 10,
		
		q =5,
		
		d=3,
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
		c.val=0.1,
		#w.vals = c(0.001,0.1,0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999),
		#w.vals = c(0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999),
		w.vals = c(0.5),
		wt.equalize=FALSE,
		rival.w = 0.999,
		power.comparison.test = F
)


w.max.index <- length(params$w.vals)	
size <- seq(0, 1, 0.01)
len <- length(size)


T0.agg<-array(0,dim=c(nmc,w.max.index,params$s))
TA.agg<-array(0,dim=c(nmc,w.max.index,params$s))


model="gaussian"

attach(params)		
pprime1     = ifelse(model=="gaussian",p+q,p+q+2)   # cca arguments , signal+noise dimension
pprime2     = ifelse(model=="gaussian",p+q,p+q+2)   # cca arguments, signal+noise dimension

pre.scaling = TRUE  #Make the measurement spaces have the same scale

alpha       = NULL

sim.grass=FALSE
eigen.spectrum=FALSE



level.mcnemar=0.01  #At what alpha, should unweighted(w=0.5) and optimal w^* be compared
def.w=0.5          #The null hypothesis is that power(def.w) >= power(rival.w) (by default ,def.w is the w for the unweighted case which equals 0.5)

proc.dilation=FALSE #when investigating convergence of JOFC to PoM, should Procrustes analysis of configurations include the dilation component?

old.gauss.model.param <- FALSE
detach(params)

param.index.count<-8
params.list <- rep(list(params),param.index.count)
p.vals<-rep(c(5,19),each=4)
for (i in (1:length(p.vals))){
params.list[[i]]$p <- p.vals[i]
#params.list[[6]]$p <- 25

}
r.vals<-rep(c(1.5,3,10,30),2)
for (i in (1:length(r.vals))){
params.list[[i]]$r <- r.vals[i]
#params.list[[6]]$p <- 25

}

#params.list[[5]]$r <- 100
vary.param<-"r"
#params.list[[1]]$q <- 22







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





#param.index<-1
for (param.index in 1:param.index.count){

oos.dist.agg.0 <- c()
oos.dist.agg.A <- c()

	n<-n.vals
	params<- params.list[[param.index]]
	
	attach(params)
	D.nmc<-array(0,dim=c(2*n,2*n,nmc))
	
	Pert.mat<-array(0,dim=c(2*n,2*n,nmc))
	
	
	
	Ideal.Gamma.matched<- rgamma  (150,shape=p/2,scale=4*(1-c.val)^2/r)+ 
			rgamma(150,shape=q/2,scale=4*c.val^2*(1+1/r))
	Ideal.Gamma.unmatched<- rgamma(150,shape=p/2,scale=4*(1-c.val)^2*(1+1/r))+ 
			rgamma(150,shape=q/2,scale=4*(c.val)^2*(1+1/r))
	
	windows()
	hist(Ideal.Gamma.matched)
	savePlot(file.path(results.dir,"Ideal.Gamma.matched.hist.pdf"),"pdf")
	dev.off()
	windows()
	plot(density(Ideal.Gamma.matched))
	savePlot(file.path(results.dir,"Ideal.Gamma.matched.dens.pdf"),"pdf")
	dev.off()
	windows()
	hist(Ideal.Gamma.unmatched)
	savePlot(file.path(results.dir,"Ideal.Gamma.unmatched.hist.pdf"),"pdf")
	dev.off()
	windows()
	plot(density(Ideal.Gamma.unmatched))
	savePlot(file.path(results.dir,"Ideal.Gamma.unmatched.dens.pdf"),"pdf")
	dev.off()
	windows()
	
	Pert.log.normal<-rgamma(150,shape=p/2, scale=0.006)
	
	Pert.Gamma.matched   <- Pert.log.normal+Ideal.Gamma.matched
	Pert.Gamma.unmatched <-   Pert.log.normal+Ideal.Gamma.matched
	oos.diss.0 <-c()
	oos.diss.A <-c()
	detach(params)
	
	
	#sink(paste("log-params-",vary.param,param.index,".txt",sep="",collapse=""))
	
	for (mc in 1:nmc){
		
		results<-with(params, {	
					power.w.star <- 0
					
					m<-s
					
					T0 <- matrix(0,w.max.index,m)   #Test statistics for JOFC under null
					TA <- matrix(0,w.max.index,m)    #Test statistics for JOFC under alternative
					Y.oos.matched <- c()
					Y.oos.unmatched<-c()
					
					
					cont.table  <- matrix(0,2,2)
					
					Fid.Err.Term.1 <- array(0,dim=c(w.max.index))
					Fid.Err.Term.2 <- array(0,dim=c(w.max.index))
					Comm.Err.Term <- array(0,dim=c(w.max.index))
					
					sigma <- matrix(0,p,p)
					means <- array(0 , dim=c(w.max.index,2*d))
					
					
					
					if (is.null(alpha)) {
						if (model=="gaussian"){
							sigma<- diag(p)
							if (old.gauss.model.param) sigma <-Posdef(p,r)
							alpha.mc <- mvrnorm(n+(2*m), rep(0,p),sigma)
						} else if (model=="dirichlet"){
							alpha.mc <- rdirichlet(n+2*m, rep(1,p+1))
						} else stop("unknown model")
						
						
					} else {
						alpha.mc <- alpha[[mc]]
					}
					
					## n pairs of matched points
					if (model=="gaussian"){
						xlist <- matched_rnorm(n, p, q, c.val, r, alpha=alpha.mc[1:n, ],sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)
					} else{
						xlist <- matched_rdirichlet(n, p, r, q, c.val, alpha.mc[1:n, ])
					}
					
					## n pairs of matched points
					if (model=="gaussian"){
						ylist <- matched_rnorm(m, p, q, c.val, r, alpha=as.array(alpha.mc[n+(1:m), ,drop=FALSE]),sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)
					} else{
						ylist <- matched_rdirichlet(m, p, r, q, c.val, alpha=as.array(alpha.mc[n+(1:m), ,drop=FALSE]))
					}
					if (model=="gaussian"){
						y.alt <- matched_rnorm(m, p, q, c.val, r, alpha=as.array((alpha.mc[n+m+(1:m), ,drop=FALSE])),sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param)
					} else{
						y.alt <- matched_rdirichlet(m, p, r, q, c.val,alpha= as.array(alpha.mc[n+m+(1:m), ,drop=FALSE]))
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
						L <- matrix(NA,n,n)
						diag(L)<- 0
					}
					
					
					#In sample embedding
					# Form omnibus dissimilarity matrix
					M <- omnibusM(D1, D2, L)
					init.conf<-NULL
					
					if (compare.pom.cca) init.conf<- pom.config
					
					# Embed in-sample using different weight matrices (differentw values)
					X.embeds<-JOFC.Fid.Commens.Tradeoff(M,d,w.vals,separability.entries.w,init.conf=init.conf,wt.equalize=wt.equalize)
					
					Pert.mat.mc<-(as.matrix(dist(X.embeds[[1]])^2))-as.matrix(D.whole^2)
					
					#Form OOS omnibus matrices
					M.oos.0<- matrix(0,2,2)#omnibusM(D.oos.1,D.oos.2.null, L.tilde.null)
					M.oos.A<- matrix(0,2,2)#omnibusM(D.oos.1,D.oos.2.alt,  L.tilde.alt)
					
					
					
					for (l in 1:w.max.index){
						if (verbose) print("OOS embedding for JOFC for w= \n")
						if (verbose) print(w.vals[l])
						
						w.val.l <- w.vals[l]
						X <- X.embeds[[l]]
						
						
						oos.obs.flag<- c(rep(1,2*n),rep(0,2*m))
						
						#Compute Weight matrix corresponding in-sample  entries
						oos.Weight.mat.1<-w.val.to.W.mat(w.val.l,(2*n),separability.entries.w,wt.equalize)
						
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
							oos.Weight.mat.w <- matrix(1-w.val.l,2*n,2*m)
						} else{
							oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n,m), matrix(0,n,m) ),
									cbind(matrix(0,n,m),matrix(1-w.val.l,n,m))
							)
						}
						oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
						# Since we are going to oos-embedding, set the weights  of in-sample embedding of stress
						# We are using previous in-sample embeddings, anyway
						oos.Weight.mat[1:(2*n),1:(2*n)]<-0
						if (verbose) print("dim(M.oos.0)")
						if (verbose) print(dim(M.oos.0))
						if (verbose) print("dim(M.oos.A)")
						if (verbose) print(dim(M.oos.A))
						if (verbose) print("dim(oos.Weight.mat)")
						if (verbose) print(dim(oos.Weight.mat))
						if (verbose) print("dim(X)")
						if (verbose) print(dim(X))
						#if (verbose) {print("oos.obs.flag")
						
						
						ideal.omnibus.0  <- as.matrix(dist(rbind(X1,X2,Y1,Y20)))
						ideal.omnibus.A  <- as.matrix(dist(rbind(X1,X2,Y1,Y2A)))
						omnibus.oos.D.0 <- omnibusM(M,M.oos.0,ideal.omnibus.0[1:(2*n),(2*n)+(1:(2*m))])
						omnibus.oos.D.A <- omnibusM(M,M.oos.A, ideal.omnibus.A[1:(2*n),(2*n)+(1:(2*m))])
						oos.Weight.mat[is.na(omnibus.oos.D.0)]<-0
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
								bwOos    = TRUE)
						
						
						embedded.dist.0<- (Y.0t[1]-Y.0t[2])^2
						if (verbose) print("JOFC alternative omnibus OOS embedding \n")
						Y.At<-oosIM(D=omnibus.oos.D.A,
								X=X,
								init     = "random",
								verbose  = FALSE,
								itmax    = 1000,
								eps      = 1e-8,
								W        = oos.Weight.mat,
								isWithin = oos.obs.flag,
								bwOos    = TRUE)
#if (profile.mode)				Rprof(NULL)
						embedded.dist.A<- (Y.At[1]-Y.At[2])^2	
						
						oos.dist.agg.0 <-c(oos.dist.agg.0, embedded.dist.0)
						oos.dist.agg.A <-c(oos.dist.agg.A, embedded.dist.A)
					}
					list(p1=Pert.mat.mc,p2=as.matrix(D.whole),p3=oos.dist.agg.0,p4=oos.dist.agg.A)
				}	
		)
		Pert.mat[,,mc]<-results[[1]]
		D.nmc   [,,mc]<-results[[2]]
		oos.diss.0<-c(oos.diss.0,results[[3]])
		oos.diss.A<-c(oos.diss.A,results[[4]])
		
	} 
	
	comm.terms <- cbind ((1:n)+n,1:n)
	comm.terms <- rbind(comm.terms,cbind ((1:n),(1:n)+n))
	sep.terms <- cbind ((1:n)+n,1:n)
	sep.terms <- rbind(comm.terms,cbind ((1:n),(1:n)+n))
	diag.terms<- cbind ((1:n)+n,1:n)
	diag.terms<- rbind(diag.terms,cbind ((1:n),1:n))
	
	
	Y.comm.agg<-c()
	r.comm.agg<-c()
	Y.sep.agg<-c()
	r.sep.agg<-c()
	
	
	#sink()
	for (mc in 1:nmc){
		Pert.mat.mc<-Pert.mat[,,mc]
		D.mc<-D.nmc[,,mc]
		Y.comm<- Pert.mat.mc[comm.terms]
		r.comm<-Y.comm/D.mc[comm.terms]
		Y.sep<- rbind(Pert.mat.mc[1:n,1:n],Pert.mat.mc[n+(1:n),n+(1:n)])
		
		r.sep<-Y.sep/rbind(D.mc[1:n,1:n],D.mc[n+(1:n),n+(1:n)])
		
		Y.sep[diag.terms]<-NA
		r.sep[diag.terms]<-NA
		
		Y.sep <- as.vector(Y.sep)
		r.sep <- as.vector(r.sep)
		Y.sep<-Y.sep[!is.na(Y.sep)]
		r.sep<-r.sep[!is.na(r.sep)]
		
		
		
		Y.comm.agg <- c(Y.comm.agg, as.vector(Y.comm))
		r.comm.agg <- c(r.comm.agg, as.vector(r.comm))
		Y.sep.agg  <- c(Y.sep.agg,Y.sep)
		r.sep.agg  <- c(r.sep.agg,r.sep)
	}
	
	
	
	
	windows()
	
	plot(density(as.vector(Y.comm.agg)))
	
	Y.comm.dist.par<- fitdistr((Y.comm.agg+abs(min(Y.comm.agg))+0.01),"exponential")
	print(Y.comm.dist.par)
	Y.comm.dist.par.g<- fitdistr((Y.comm.agg+abs(min(Y.comm.agg))+0.01),"gamma")
	Y.comm.dist.par.log.n<-fitdistr((Y.comm.agg+abs(min(Y.comm.agg))+0.01),"lognormal")
	
	x<-seq(min(Y.comm.agg)-0.1,max(Y.comm.agg)+abs(min(Y.comm.agg))+0.1,0.1)
	#lines(x,dexp(x,rate=Y.comm.dist.par$estimate["rate"]), col="blue")
	lines(x,dgamma(x,shape=Y.comm.dist.par.g$estimate["shape"],rate=Y.comm.dist.par.g$estimate["rate"]),col="red")
	lines(x,dlnorm(x,meanlog=Y.comm.dist.par.log.n$estimate["meanlog"],sdlog=Y.comm.dist.par.log.n$estimate["sdlog"]),col="blue")
	
  savePlot(file.path(results.dir,"Y.comm.perturb.dist"),"pdf")
  
	
	
	
	
	windows()
	plot(density(as.vector(Y.sep.agg)))
	
	

	Y.sep.dist.par<- fitdistr(Y.sep.agg,"normal")
	print(Y.sep.dist.par)
	x<-seq(min(Y.sep.agg)-0.1,max(Y.sep.agg)+0.1,0.1)
	lines(x,dnorm(x,mean=Y.sep.dist.par$estimate["mean"],sd=Y.sep.dist.par$estimate["sd"]),col="blue")
	
	
	#laplace(llocation = "identity", lscale = "loge", elocation = list(),
	#		escale = list(), ilocation = NULL, iscale = NULL,
	#		imethod = 1, zero = 2)
	#fit = vglm(y ~ 1, laplace, lddat, trace = TRUE, crit = "l")
	savePlot(file.path(results.dir,"Y.sep.perturb.dist"),"pdf")
  


     windows()
	plot(density(as.vector(r.comm.agg)))
	
	x<-seq(min(r.comm.agg)-0.1,max(r.comm.agg)+abs(min(r.comm.agg))+0.1,0.1)
	r.comm.dist.par<- fitdistr(r.comm.agg+abs(min(r.comm.agg))+0.01,"exponential")
	lines(x,dexp(x,rate=r.comm.dist.par$estimate["rate"]),col="red")
	r.comm.dist.par.log.n<-fitdistr((r.comm.agg+abs(min(r.comm.agg))+0.01),"lognormal")
	
	lines(x,dlnorm(x,meanlog=r.comm.dist.par.log.n$estimate["meanlog"],
					sdlog=r.comm.dist.par.log.n$estimate["sdlog"]),col="blue")

  savePlot(file.path(results.dir,"r.ratio.comm.perturb.dist"),"pdf")
 
	
	print(r.comm.dist.par)
	print(r.comm.dist.par.log.n)
	
	
	
	
	
	windows()
	plot(density(as.vector(r.sep.agg)))
	x<-seq(min(r.sep.agg)-0.1,max(r.sep.agg)+0.1,0.1)
	r.sep.dist.par<- fitdistr(r.sep.agg,"normal")
	lines(x,dnorm(x,mean=r.sep.dist.par$estimate["mean"],sd=r.sep.dist.par$estimate["sd"]),col="blue")
	print(r.sep.dist.par)
	  savePlot(file.path(results.dir,"r.ratio.sep.perturb.dist"),"pdf")
 
	
	
	plot(density(as.vector(oos.diss.0)))
	oos.diss.0.exp.dist.par<- fitdistr(as.vector(oos.diss.0+abs(min(oos.diss.0))+0.01),"exponential")
	
	
	x<-seq(min(oos.diss.0)-0.1,max(oos.diss.0)+abs(min(oos.diss.0))+0.1,0.1)
	lines(x,dexp(x,rate=oos.diss.0.exp.dist.par$estimate["rate"]),col="red")
	oos.diss.0.par.log.n<-fitdistr((oos.diss.0+abs(min(oos.diss.0))+0.01),"lognormal")
	
	lines(x,dlnorm(x,meanlog=oos.diss.0.par.log.n$estimate["meanlog"],
					sdlog=oos.diss.0.par.log.n$estimate["sdlog"]),col="blue")
	  savePlot(file.path(results.dir,"oos.perturb.for.nulldist"),"pdf")
 
	
	
	
	windows()
	em.result<-em("E",oos.diss.0,parameters=list(pro=c(0.9,0.1),mean=c(0,0.2),variance=list(modelName="V",d=1,G=2,sigmasq=c(.04))))
	oos.diss.0.dist.par<- em.result$parameters
	print(oos.diss.0.dist.par)
	
	
	
	plot(density(as.vector(oos.diss.A)))
	oos.diss.A.exp.dist.par<- fitdistr(as.vector(oos.diss.A+abs(min(oos.diss.A))+0.01),"exponential")
	x<-seq(min(oos.diss.A)-0.1,max(oos.diss.A)+0.1,0.1)
	lines(x,dexp(x,rate=oos.diss.A.exp.dist.par$estimate["rate"]),col="red")
	oos.diss.A.par.log.n<-fitdistr((oos.diss.A+abs(min(oos.diss.A))+0.01),"lognormal")
	
	lines(x,dlnorm(x,meanlog=oos.diss.A.par.log.n$estimate["meanlog"],
					sdlog=oos.diss.A.par.log.n$estimate["sdlog"]),col="blue")

  savePlot(file.path(results.dir,"oos.perturb.for.altdist"),"pdf")
 
	windows()
	em.result<-em("V",oos.diss.A,parameters=list(pro=c(0.9,0.1),mean=c(0,8),variance=list(modelName="V",d=1,G=2,sigmasq=c(3,3))))
	oos.diss.A.dist.par<- em.result$parameters
	print(oos.diss.A.dist.par)

	
	
	
	Y.comm.dist.par.agg<- c(Y.comm.dist.par.agg,Y.comm.dist.par)
	
	Y.comm.dist.par.g.agg<-c(Y.comm.dist.par.g.agg,Y.comm.dist.par.g)
	Y.comm.dist.par.log.n.agg<-c(Y.comm.dist.par.log.n.agg,Y.comm.dist.par.log.n)
	Y.sep.dist.par.agg<- c(Y.sep.dist.par.agg,Y.sep.dist.par)
	
#laplace(llocation = "identity", lscale = "loge", elocation = list(),
#		escale = list(), ilocation = NULL, iscale = NULL,
#		imethod = 1, zero = 2)
#fit = vglm(y ~ 1, laplace, lddat, trace = TRUE, crit = "l")
	
	
	
	r.comm.dist.par.agg<-c(r.comm.dist.par.agg,r.comm.dist.par)
	
	r.comm.dist.par.log.n.agg<-c(r.comm.dist.par.log.n.agg,r.comm.dist.par.log.n)
	r.sep.dist.par.agg<- c(r.sep.dist.par.agg,r.sep.dist.par)
	
	
	
	oos.diss.0.exp.dist.par.agg<-c(oos.diss.0.exp.dist.par.agg,oos.diss.0.exp.dist.par)
	
	
	oos.diss.A.exp.dist.par.agg<- c(oos.diss.A.exp.dist.par.agg,oos.diss.A.exp.dist.par)
	oos.diss.0.par.log.n.agg<- c(oos.diss.0.par.log.n.agg,oos.diss.0.par.log.n)
	oos.diss.A.par.log.n.agg<-c(oos.diss.A.par.log.n.agg,oos.diss.A.par.log.n)
	

	means.0<- c(means.0,mean(oos.diss.0))
	vars.0 <- c(vars.0,var(oos.diss.0))
	skews.0 <- c(skews.0,skewness(oos.diss.0))
	kurts.0 <- c(kurts.0,kurtosis(oos.diss.0))

	means.A<- c(means.A,mean(oos.diss.A))
	vars.A <- c(vars.A,var(oos.diss.A))
	skews.A <- c(skews.A,skewness(oos.diss.A))
	kurts.A <- c(kurts.A,kurtosis(oos.diss.A))

	
	
	}
#	title(plot.title)
	
#	fname<- paste(plot.title,".pdf",collapse="")
#	savePlot(filename=fname,type="pdf")
	
#	save.image(paste("params",vary.param,param.index,".RData",sep="",collapse=""))
#}

