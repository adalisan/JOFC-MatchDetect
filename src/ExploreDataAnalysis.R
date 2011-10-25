# TODO: Add comment
# 
# Author: Sancar
###############################################################################


source("./src/simulation_math_util_fn.R")
source("./src/gaussian_simulation_sim_fn.R")
source("./src/dirichlet_simulation_fn.R")
source("./src/oosMDS.R")
source("./src/smacofM.R")
source("./src/oosIM.R")
source("./src/simulateTestPlot.R")



results.dir <- "./results/"
require(MASS)
require(MCMCpack)
require(exact2x2)

par.compute <- FALSE
run.in.linux<- FALSE
compare.pom.cca<-FALSE
#c.vals<-c(0.01)
c.vals<-c(0.01)

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


	color.file<- read.csv("Cat_12.csv",header=FALSE,skip=2,as.is=TRUE)
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








#n.vals <-c(50,100,150,200,300,500)
#
n.vals<-c(150)
nmc <-  4
s<- 150

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
		p = 7,
		r = 10,
		
		q =12,
		
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
		c.val=0.01,
		#w.vals = c(0.001,0.1,0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999),
		w.vals = c(0.5,0.8,0.85,0.9,0.925,0.95,0.99,0.999),
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

param.index.count<-1
params.list <- rep(list(params),param.index.count)

#params.list[[1]]$p <- 3
#params.list[[2]]$p <- 5
#params.list[[3]]$p <- 10
#params.list[[4]]$p <- 15
#params.list[[5]]$p <- 19
#params.list[[6]]$p <- 25
vary.param<-"p"

#params.list[[1]]$r <- 1.5
#params.list[[2]]$r <- 3
#params.list[[3]]$r <- 10
#params.list[[4]]$r <- 30
#params.list[[5]]$r <- 100

#params.list[[1]]$q <- 80





FC1.agg<- array(0,dim=c(nmc,w.max.index))
FC2.agg<- array(0,dim=c(nmc,w.max.index))
FC3.agg<- array(0,dim=c(nmc,w.max.index))


for (param.index in 1:param.index.count){

params<- params.list[[param.index]]

sink(paste("log-params-",vary.param,param.index,".txt",sep="",collapse=""))

for (mc in 1:nmc){

Test.Stats<-with(params, {	
power.w.star <- 0

m<-s

power.mc= array(0,dim=c(w.max.index,len))  #power values for JOFC in this MC replicate
power.cca.mc = array(0,dim=c(len))         #power values for CCA in this MC replicate
power.pom.mc = array(0,dim=c(len))         #power values for PoM in this MC replicate
power.cca.reg.mc = array(0,dim=c(len))     #power values for reg CCA in this MC replicate


config.mismatch <-  list(frob.norm=array(0,dim=c(w.max.index))) #Frob. norm of configuration difference
#between PoM and JOFC with smallest w
min.stress.for.w.val = array(0,dim=c(w.max.index))   #minimum stress value for  smacof algorithm
pom.stress <- 0

T0.cca.reg <- array(0,dim=c(m))     #Test statistics for regularized CCA under null
TA.cca.reg <- array(0,dim=c(m))		#Test statistics for regularized CCA under alternative

T0.cca <- array(0,dim=c(m))     #Test statistics for CCA under null
TA.cca <- array(0,dim=c(m))		#Test statistics for CCA under alternative

T0.pom <- array(0,dim=c(m))    #Test statistics for PoM under null
TA.pom <- array(0,dim=c(m))    #Test statistics for JOFC under alternative

T0 <- matrix(0,w.max.index,m)   #Test statistics for JOFC under null
TA <- matrix(0,w.max.index,m)    #Test statistics for JOFC under alternative

T0.best.w <- matrix(0,2,m)    #Test statistics for JOFC (comparison of w=0.5 with optimal w*
TA.best.w <- matrix(0,2,m)

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
X1 <- xlist$X1
X2 <- xlist$X2
if (model=="gaussian")
	sigma.mc<-xlist$sigma.beta

D1 <- dist(X1)
D2 <- dist(X2)

if (verbose) print("random matched pairs generated\n")

#prescaling
if (pre.scaling) {
	s <- lm(as.vector(D1) ~ as.vector(D2) + 0)$coefficients
} else {
	s <- 1
}

#m pairs of unmatched points
if (model=="gaussian"){
	## test observations -- m pairs of matched and m pairs of unmatched
	ylist <- matched_rnorm(m, p, q, c.val, r, alpha=alpha.mc[(n+1):(n+m), ],
			sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)
	Y2A <- matched_rnorm(m, p, q, c.val, r, alpha=alpha.mc[(n+m+1):(n+m+m), ],
			sigma.alpha=sigma,old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)$X2
} else{
	ylist <- matched_rdirichlet(m, p, r, q, c.val, alpha.mc[(n+1):(n+m), ])
	Y2A <- matched_rdirichlet(m, p, r, q, c.val, alpha.mc[(n+m+1):(n+m+m), ])$X2
}
Y1 <- ylist$X1
Y20 <- ylist$X2

# Dissimilarity matrices for in-sample +out-of-sample
D10A <- as.matrix(dist(rbind(X1, Y1)))

D20 <- as.matrix(dist(rbind(X2, Y20))) * s
D2A <- as.matrix(dist(rbind(X2, Y2A))) * s
D1<-as.matrix(D1)
D2<-as.matrix(D2)
pom.config<-c()
cca.config<-c()
D2<-D2*s

if (verbose) print("PoM and CCA embedding\n")	
if (compare.pom.cca) {
	
	## ==== cca ====
	#embed in-sample measurements
	if (oos == TRUE) {
		if (c.val==0){
			if (model=="gaussian"){
				
				X1t <- smacofM(D1,ndim = p,verbose=FALSE)
				X2t <- smacofM(D2,ndim = p,verbose=FALSE)
			} else{
				X1t <- smacofM(D1,ndim = p+1,verbose=TRUE)		
				X2t <- smacofM(D2,ndim = p+1,verbose=FALSE)
			}
		} else{
			X1t <- smacofM(D=D1,ndim= pprime1,verbose=FALSE)
			X2t <- smacofM(D=D2,ndim= pprime2,verbose=FALSE)
		}
		
		xcca <- cancor(X1t, X2t)
		
		#project using projection vectors computed by CCA
		Y1t  <- (oosMDS(D10A, X1t) %*% xcca$xcoef)[, 1:d]
		Y20t <- (oosMDS(D20, X2t) %*% xcca$ycoef)[, 1:d]
		Y2At <- (oosMDS(D2A, X2t) %*% xcca$ycoef)[, 1:d]
		#cca.config<-rbind(X1t,X2t)
		
	} else {
		if (c.val==0){
			if (model=="gaussian"){
				X1t <- smacofM(D10A, ndim=p,verbose=FALSE)
				D20A <-dist(rbind(X2, Y20, Y2A))
				X2t <- smacofM(D20A, ndim=p,verbose=FALSE)
			}
			else{
				X1t <- smacofM(D10A, ndim=p+1,verbose=FALSE)
				D20A <-dist(rbind(X2, Y20, Y2A))
				X2t <- smacofM(D20A, ndim=p+1,verbose=FALSE)
				
				
			}
		} else{
			if (model=="gaussian"){
				pprime1 <- p+q
				pprime2 <- p+q
			}
			else{
				pprime1 <- p+q+2
				pprime2 <- p+q+2
				
			}
			X1t <- smacofM(D10A, ndim=pprime1,verbose=FALSE,init=cmdscale(D10A,pprime1))
			D20A <-dist(rbind(X2, Y20, Y2A))
			X2t <- smacofM(D20A, ndim=pprime2,verbose=FALSE,init=cmdscale(D20A,pprime2))
			
			
		}
		
		
		if (verbose) print("CCA embedding complete\n")
		center1 <- colMeans(X1t[1:n, ])   # column means of training obs
		center2 <- colMeans(X2t[1:n, ])
		X1t <- X1t - matrix(center1, n+m, pprime1, byrow=TRUE) # column-center training only
		X2t <- X2t - matrix(center2, n+2*m, pprime2, byrow=TRUE)
		cca <- cancor(X1t[1:n, ], X2t[1:n, ])
		Y1t <-  (X1t[(n+1):(n+m), ] %*% cca$xcoef )[, 1:d]
		Y20t <- (X2t[(n+1):(n+m), ] %*% cca$ycoef)[, 1:d]
		Y2At <- (X2t[(n+m+1):(n+2*m), ] %*% cca$ycoef)[, 1:d]
	}
	T0.cca <- rowSums((Y1t - Y20t)^2)
	TA.cca <- rowSums((Y1t - Y2At)^2)
	power.cca.mc <- get_power(T0.cca, TA.cca, size)
	
	
	if (verbose) print("CCA test statistic complete\n")
	
	##low-dimensional (regularized) CCA 
	## ==== cca ====
	#embed in-sample measurements
	if (oos == TRUE) {
		if (c.val==0){
			if (model=="gaussian"){
				
				X1t <- smacofM(D1,ndim = floor((d+p)/2),verbose=FALSE)
				X2t <- smacofM(D2,ndim = floor((d+p)/2),verbose=FALSE)
			} else{
				X1t <- smacofM(D1,ndim = floor((d+p)/2)+1,verbose=TRUE)		
				X2t <- smacofM(D2,ndim = floor((d+p)/2)+1,verbose=FALSE)
			}
		} else{
			X1t <- smacofM(D=D1,ndim= floor((d+p)/2)+1,verbose=FALSE)
			X2t <- smacofM(D=D2,ndim= floor((d+p)/2)+1,verbose=FALSE)
		}
		
		xcca <- cancor(X1t, X2t)
		
		#project using projection vectors computed by CCA
		Y1t  <- (oosMDS(D10A, X1t) %*% xcca$xcoef)[, 1:d]
		Y20t <- (oosMDS(D20, X2t) %*% xcca$ycoef)[, 1:d]
		Y2At <- (oosMDS(D2A, X2t) %*% xcca$ycoef)[, 1:d]
		#cca.config<-rbind(X1t,X2t)
		
	} else {
		if (c.val==0){
			if (model=="gaussian"){
				X1t <- smacofM(D10A, ndim=floor((d+p)/2),verbose=FALSE)
				D20A <-dist(rbind(X2, Y20, Y2A))
				X2t <- smacofM(D20A, ndim=floor((d+p)/2),verbose=FALSE)
			}
			else{
				X1t <- smacofM(D10A, ndim= floor((d+p)/2)+1,verbose=FALSE)
				D20A <-dist(rbind(X2, Y20, Y2A))
				X2t <- smacofM(D20A, ndim= floor((d+p)/2)+1,verbose=FALSE)
				
				
			}
		} else{
			if (model=="gaussian"){
				pprime1 <- p+q
				pprime2 <- p+q
			}
			else{
				pprime1 <- p+q+2
				pprime2 <- p+q+2
				
			}
			X1t <- smacofM(D10A, ndim=pprime1,verbose=FALSE,init=cmdscale(D10A,pprime1))
			D20A <-dist(rbind(X2, Y20, Y2A))
			X2t <- smacofM(D20A, ndim=pprime2,verbose=FALSE,init=cmdscale(D20A,pprime2))
			
			
		}
		
		
		if (verbose) print("CCA embedding complete\n")
		center1 <- colMeans(X1t[1:n, ])   # column means of training obs
		center2 <- colMeans(X2t[1:n, ])
		X1t <- X1t - matrix(center1, n+m, pprime1, byrow=TRUE) # column-center training only
		X2t <- X2t - matrix(center2, n+2*m, pprime2, byrow=TRUE)
		cca <- cancor(X1t[1:n, ], X2t[1:n, ])
		Y1t <-  (X1t[(n+1):(n+m), ] %*% cca$xcoef )[, 1:d]
		Y20t <- (X2t[(n+1):(n+m), ] %*% cca$ycoef)[, 1:d]
		Y2At <- (X2t[(n+m+1):(n+2*m), ] %*% cca$ycoef)[, 1:d]
	}
	T0.cca.reg <- rowSums((Y1t - Y20t)^2)
	TA.cca.reg <- rowSums((Y1t - Y2At)^2)
	power.cca.reg.mc <- get_power(T0.cca.reg, TA.cca.reg, size)
	
	
	## ==== pom = procrustes o mds ====
	if (oos == TRUE) {
		#Embed in-sample
		X1t <- smacofM(D1, ndim=d,verbose=FALSE)
		X2t <- smacofM(D2, ndim=d,verbose=FALSE)
		if (verbose) print (colMeans(X1t))
		if (verbose) print (colMeans(X2t))
		# Compute Proc from in-sample embeddings
		proc <- MCMCpack::procrustes(X2t, X1t, dilation=proc.dilation)
		# Out-of sample embed and Proc Transform dissimilarities
		Y1t  <- oosMDS(D10A, X1t)
		Y20t <- oosMDS(D20, X2t) %*% proc$R * proc$s
		Y2At <- oosMDS(D2A, X2t) %*% proc$R * proc$s
		X2tp<-X2t %*% proc$R * proc$s
		pom.config<-rbind(X1t,X2tp)
		pom.stress<- sum((as.dist(D1) - dist(X1t))^2)
		pom.stress<- pom.stress+ sum((as.dist(D2) - dist(X2tp))^2)
		if (verbose) print("PoM embedding complete\n")
	} else {
		X1t <- smacofM(D10A,ndim= d,verbose=FALSE,init=cmdscale(D10A,d))
		D20A <-dist(rbind(X2, Y20, Y2A))
		X2t <- smacofM(D20A,ndim= d,verbose=FALSE,init=cmdscale(D20A,d))
		center1 <- colMeans(X1t[1:n, ])
		center2 <- colMeans(X2t[1:n, ])
		X1t <- X1t - matrix(center1, n+m, d, byrow=TRUE) # column-center training only
		X2t <- X2t - matrix(center2, n+2*m, d, byrow=TRUE)
		proc <- MCMCpack::procrustes(X2t[1:n, ], X1t[1:n, ], dilation=proc.dilation)
		Y1t <- X1t[(n+1):(n+m), ]
		Y20t <- X2t[(n+1):(n+m), ] %*% proc$R * proc$s
		Y2At <- X2t[(n+m+1):(n+2*m), ] %*% proc$R * proc$s
	}
	
	T0.pom <- rowSums((Y1t - Y20t)^2)
	TA.pom <- rowSums((Y1t - Y2At)^2)
	power.pom.mc <- get_power(T0.pom, TA.pom, size)
	if (verbose) print("PoM test statistic complete \n")
	
	
}




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


if (oos == TRUE) {
	
	#In sample embedding
	# Form omnibus dissimilarity matrix
	M <- omnibusM(D1, D2, L)
	init.conf<-NULL
	
	if (compare.pom.cca) init.conf<- pom.config
	
	# Embed in-sample using different weight matrices (differentw values)
	X.embeds<-JOFC.Fid.Commens.Tradeoff(M,d,w.vals,separability.entries.w,init.conf=init.conf,wt.equalize=wt.equalize)
	
	Fid.Err.Term.1 <- X.embeds[[w.max.index+2]]
	Fid.Err.Term.2 <- X.embeds[[w.max.index+3]]
	Comm.Err.Term  <- X.embeds[[w.max.index+4]]
	
	Fid.Err.Sum.Term.1 <- X.embeds[[w.max.index+5]]
	Fid.Err.Sum.Term.2 <- X.embeds[[w.max.index+6]]
	Comm.Err.Sum.Term  <- X.embeds[[w.max.index+7]]
	FC.ratio  <- X.embeds[[w.max.index+8]]
	FC.ratio.2  <- X.embeds[[w.max.index+9]]
	FC.ratio.3  <- X.embeds[[w.max.index+10]]
	
	print("Fid.Err.Term.1" )
	print(Fid.Err.Term.1 )
	print("Comm.Err.Term ")
	print(Comm.Err.Term )
	min.stress.for.w.val <- X.embeds[[w.max.index+1]]
	if (verbose) print("JOFC embeddings complete\n")
	
	
	#
	# OOS Dissimilarity matrices
	#
	
	
	
	D.oos.1<-dist(Y1)
	D.oos.2.null <- dist(Y20)*s
	D.oos.2.alt <- dist(Y2A)*s
	
	#Imputing dissimilarity  entries for OOS
	if (Wchoice == "avg") {
		L.tilde.null <- (D.oos.1 + D.oos.2.null)/2
		L.tilde.alt  <- (D.oos.1 + D.oos.2.alt)/2
	} else if (Wchoice == "sqrt") {
		L.tilde.null <- sqrt((D.oos.1^2 + D.oos.2.null^2)/2)
		L.tilde.alt  <- sqrt((D.oos.1^2 + D.oos.2.alt^2)/2)
		
	} else if (Wchoice == "NA+diag(0)") {
		L.tilde.null <- matrix(NA,m,m)
		L.tilde.alt <- matrix(NA,m,m)
		diag(L.tilde.null)<- 0
		diag(L.tilde.alt)<- 0
	}
	#Form OOS omnibus matrices
	M.oos.0<- omnibusM(D.oos.1,D.oos.2.null, L.tilde.null)
	M.oos.A<- omnibusM(D.oos.1,D.oos.2.alt,  L.tilde.alt)
	
	
	
	for (l in 1:w.max.index){
		if (verbose) print("OOS embedding for JOFC for w= \n")
		if (verbose) print(w.vals[l])
		
		w.val.l <- w.vals[l]
		X <- X.embeds[[l]]
		
		
		oos.obs.flag<- c(rep(1,2*n),rep(0,2*m))
		
		#Compute Weight matrix corresponding in-sample  entries
		oos.Weight.mat.1<-w.val.to.W.mat(w.val.l,(2*n),separability.entries.w,wt.equalize)
		
		#Compute Weight matrix corresponding OOS  entries
		oos.Weight.mat.2<-w.val.to.W.mat(w.val.l,(2*m),separability.entries.w,wt.equalize)
		
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
		Y.0t<-oosIM(D=omnibus.oos.D.0,
				X=X,
				init     = "random",
				verbose  = FALSE,
				itmax    = 1000,
				eps      = 1e-8,
				W        = oos.Weight.mat,
				isWithin = oos.obs.flag,
				bwOos    = TRUE)
		
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
		Y1t<-Y.0t[1:m,]
		Y2t<-Y.0t[m+(1:m),]
		Y1t.A<-Y.At[1:m,]
		Y2At<-Y.At[m+(1:m),]
		
		if (compare.pom.cca){
		X2tp<-pom.config[n+(1:n),]
		X1t<-pom.config[(1:n),]
		
		X.0<-rbind(X1t,X2tp)
		X.a<-X[1:n,]
		X.b<-X[n+(1:n),]
		mean.a <- colMeans(X.a)
		mean.b <- colMeans(X.b)
		#	means[l,]<- c(mean.a,mean.b)
		proc.pom2JOFC <- MCMCpack::procrustes(X,X.0,dilation=FALSE,translation=TRUE)
		#proc.pom2JOFC.a <- MCMCpack::procrustes(X.a,X1t,dilation=FALSE,translation=TRUE)
		
		#proc.pom2JOFC.b <- MCMCpack::procrustes(X.b,X2tp,dilation=FALSE,translation=TRUE)
		#X.c<-rbind(X.a-mean.a,X.b-mean.b)
		#proc.pom2JOFC.a <- MCMCpack::procrustes(X.c,X.0,dilation=FALSE,translation=TRUE)
		
		config.mismatch$frob.norm[l] <- norm (proc.pom2JOFC$X.new-X.0,'F')
		#config.mismatch[l,2] <- norm (proc.pom2JOFC.a$X.new-X.0,'F')
		#config.mismatch[l,3] <- norm (proc.pom2JOFC.b$X.new-X2tp,'F')
		
		}
		#	if (verbose) print(means[l,])
		
		
		
		
		T0[l,] <- rowSums((Y1t - Y2t)^2)
		TA[l,] <- rowSums((Y1t.A - Y2At)^2)
		if (verbose) print("JOFC test statistic complete \n")
		power.mcnemar.l <- get_power(T0[l,],TA[l,],level.mcnemar)
		if (power.mcnemar.l>power.w.star){
			rival.w <- w.vals[l]
			power.w.star <- power.mcnemar.l
			w.val.rival.idx <- l
		}
		
	}
}
else {
	M0 <- omnibusM(D10A, D20, W=(D10A+D20)/2)
	MA <- omnibusM(D10A, D2A, W=(D10A+D2A)/2)
	
	X0.embeds<-JOFC.Fid.Commens.Tradeoff(M0,d,w.vals,separability.entries.w,wt.equalize=wt.equalize)
	XA.embeds<-JOFC.Fid.Commens.Tradeoff(MA,d,w.vals,separability.entries.w,wt.equalize=wt.equalize)
	for (l in 1:w.max.index){
		X0 <- X0.embeds[[l]]
		XA <- XA.embeds[[l]]
		T0[l,] <- rowSums((X0[(n+1):(n+m), ] - X0[(n+m+n+1):(n+m+n+m), ])^2)
		TA[l,] <- rowSums((XA[(n+1):(n+m), ] - XA[(n+m+n+1):(n+m+n+m), ])^2)
		#Not done yet
		#if (compare.pom.cca){
		#X.0<-rbind(X1t,X2t %*% proc$R * proc$s)
		#proc.pom2JOFC <- MCMCpack::procrustes(X,X.0,dilation=FALSE)
		#config.mismatch[l] <- norm (proc.pom2JOFC$X.new-X.0,'F')
		#}
		
	}
	
}


# Power comparison test
# In order to compare the best w^* vs w=0.5 in an unbiased way
# re-run the simulation only for w= w^* and w=0.5
# compute the contingency table using those results
#if (power.comparison.test){
## n pairs of matched points
if (model=="gaussian"){
	xlist <- matched_rnorm(n, p, q, c.val, r, alpha=alpha.mc[1:n, ],sigma.alpha=sigma,
			old.gauss.model.param=old.gauss.model.param, sigma.beta=sigma.mc)
} else{
	xlist <- matched_rdirichlet(n, p, r, q, c.val, alpha.mc[1:n, ])
}
X1 <- xlist$X1
X2 <- xlist$X2
D1 <- dist(X1)
D2 <- dist(X2)

if (verbose) print("random matched pairs generated\n")
if (pre.scaling) {
	s <- lm(as.vector(D1) ~ as.vector(D2) + 0)$coefficients
} else {
	s <- 1
}

if (model=="gaussian"){
	## test observations -- m pairs of matched and m pairs of unmatched
	ylist <- matched_rnorm(m, p, q, c.val, r, alpha=alpha.mc[(n+1):(n+m), ],sigma.alpha=sigma,
			old.gauss.model.param=old.gauss.model.param,sigma.beta=sigma.mc)
	Y2A <- matched_rnorm(m, p, q, c.val, r, alpha=alpha.mc[(n+m+1):(n+m+m), ],sigma.alpha=sigma,
			old.gauss.model.param=old.gauss.model.param,sigma.beta=sigma.mc)$X2
} else{
	ylist <- matched_rdirichlet(m, p, r, q, c.val, alpha.mc[(n+1):(n+m), ])
	Y2A <- matched_rdirichlet(m, p, r, q, c.val, alpha.mc[(n+m+1):(n+m+m), ])$X2
}
Y1 <- ylist$X1
Y20 <- ylist$X2


D10A <- as.matrix(dist(rbind(X1, Y1)))

D20 <- as.matrix(dist(rbind(X2, Y20))) * s
D2A <- as.matrix(dist(rbind(X2, Y2A))) * s
D1<-as.matrix(D1)
D2<-as.matrix(D2)

D2<-D2*s



## ==== jofc ====
if (Wchoice == "avg") {
	L <- (D1 + D2)/2
} else if (Wchoice == "sqrt") {
	L <- sqrt((D1^2 + D2^2)/2)
} else if (Wchoice == "NA+diag(0)") {
	L <- matrix(NA,n,n)
	diag(L)<- 0
}


if (oos == TRUE) {
	
	#In sample embedding
	M <- omnibusM(D1, D2, L)
	init.conf<-NULL
	
	if (compare.pom.cca) init.conf<- pom.config
	
	#
	# Use only def.w=0.5 and rival.w for w.vals
	X.embeds.compare<-JOFC.Fid.Commens.Tradeoff(M,d,c(def.w,rival.w),separability.entries.w,init.conf=init.conf,wt.equalize=wt.equalize)
	
	if (verbose) print("JOFC embeddings complete\n")
	
	
	#
	# OOS Dissimilarity matrices
	#
	
	
	
	D.oos.1<-dist(Y1)
	D.oos.2.null <- dist(Y20)
	D.oos.2.alt <- dist(Y2A)
	
	if (Wchoice == "avg") {
		L.tilde.null <- (D.oos.1 + D.oos.2.null)/2
		L.tilde.alt  <- (D.oos.1 + D.oos.2.alt)/2
	} else if (Wchoice == "sqrt") {
		L.tilde.null <- sqrt((D.oos.1^2 + D.oos.2.null^2)/2)
		L.tilde.alt  <- sqrt((D.oos.1^2 + D.oos.2.alt^2)/2)
		
	} else if (Wchoice == "NA+diag(0)") {
		L.tilde.null <- matrix(NA,m,m)
		L.tilde.alt <- matrix(NA,m,m)
		diag(L.tilde.null)<- 0
		diag(L.tilde.alt)<- 0
	}
	M.oos.0<- omnibusM(D.oos.1,D.oos.2.null, L.tilde.null)
	M.oos.A<- omnibusM(D.oos.1,D.oos.2.alt,  L.tilde.alt)
	
	
	for (l in 1:2){
		if (verbose) print(paste(rival.w))
		if (l==1){
			w.val.l <- def.w
		}
		else {
			w.val.l <- rival.w
		}
		X<-X.embeds.compare[[l]]
		
		oos.obs.flag<- c(rep(1,2*n),rep(0,2*m))
		
		oos.Weight.mat.1<-w.val.to.W.mat(w.val.l,(2*n),separability.entries.w,wt.equalize)
		
		oos.Weight.mat.2<-w.val.to.W.mat(w.val.l,(2*m),separability.entries.w,wt.equalize)
		if (!assume.matched.for.oos){
			oos.Weight.mat.2[1:m,m+(1:m)]<-0
			oos.Weight.mat.2[m+(1:m),(1:m)]<-0
		}
		if (oos.use.imputed){
			oos.Weight.mat.w <- matrix(1-w.val.l,2*n,2*m)
		} else{
			oos.Weight.mat.w <- rbind(cbind(matrix(1-w.val.l,n,m), matrix(0,n,m) ),
					cbind(matrix(0,n,m),matrix(1-w.val.l,n,m))
			)
		}
		oos.Weight.mat<-omnibusM(oos.Weight.mat.1,oos.Weight.mat.2,oos.Weight.mat.w)
		oos.Weight.mat[1:(2*n),1:(2*n)]<-0
		
		
		
		ideal.omnibus.0  <- as.matrix(dist(rbind(X1,X2,Y1,Y20)))
		ideal.omnibus.A  <- as.matrix(dist(rbind(X1,X2,Y1,Y2A)))
		omnibus.oos.D.0 <- omnibusM(M,M.oos.0,ideal.omnibus.0[1:(2*n),(2*n)+(1:(2*m))])
		omnibus.oos.D.A <- omnibusM(M,M.oos.A, ideal.omnibus.A[1:(2*n),(2*n)+(1:(2*m))])
		oos.Weight.mat[is.na(omnibus.oos.D.0)]<-0
		omnibus.oos.D.0[is.na(omnibus.oos.D.0)]<-1
		omnibus.oos.D.A[is.na(omnibus.oos.D.A)]<-1
		if (verbose) print("JOFC null omnibus OOS embedding \n")
		Y.0t<-oosIM(D=omnibus.oos.D.0,
				X=X,
				init     = "random",
				verbose  = FALSE,
				itmax    = 1000,
				eps      = 1e-8,
				W        = oos.Weight.mat,
				isWithin = oos.obs.flag,
				bwOos    = TRUE)
		
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
		Y1t<-Y.0t[1:m,]
		Y2t<-Y.0t[m+(1:m),]
		Y1t.A<-Y.At[1:m,]
		Y2At<-Y.At[m+(1:m),]
		
		
		
		T0.best.w[l,] <- rowSums((Y1t - Y2t)^2)
		TA.best.w[l,] <- rowSums((Y1t.A - Y2At)^2)
		if (verbose) print("JOFC test statistic complete \n")
		
	}
	
}
w.val.def.idx <- which(w.vals==def.w)
w.val.rival.idx<- which(w.vals==rival.w)
crit.value<-get_crit_val(T0.best.w[1,],level.mcnemar)
crit.value.2<-get_crit_val(T0.best.w[2,],level.mcnemar)
if (verbose){
	print("crit.values")
	print(crit.value)
	print(crit.value.2)
}
cont.table[1,1] <- sum(T0.best.w[1,]<=crit.value & T0.best.w[2,]<=crit.value.2) + 
		sum(TA.best.w[1,]>crit.value & TA.best.w[2,]>crit.value.2)
cont.table[1,2] <- sum(T0.best.w[1,]>crit.value & T0.best.w[2,]<=crit.value.2)  + 
		sum(TA.best.w[1,]<=crit.value & TA.best.w[2,]>crit.value.2)
cont.table[2,1] <- sum(T0.best.w[1,]<=crit.value & T0.best.w[2,]>crit.value.2)  +
		sum(TA.best.w[1,]>crit.value & TA.best.w[2,]<=crit.value.2)
cont.table[2,2] <- sum(T0.best.w[1,]>crit.value & T0.best.w[2,]>crit.value.2)   + 
		sum(TA.best.w[1,]<=crit.value & TA.best.w[2,]<=crit.value.2) 
if (verbose) print("Cont table computed \n")
if (verbose) print(cont.table)

#	}
for (w.ind in 1:w.max.index){
	power.mc[w.ind, ] <- get_power(T0[w.ind,], TA[w.ind,], size)
	#T0.agg[mc,w.ind,] <- T0[w.ind,]
	#TA.agg[mc,w.ind,] <- TA[w.ind,]
	
}

FidComm.Terms<- list(F1=Fid.Err.Term.1,F2=Fid.Err.Term.2,C=Comm.Err.Term)
FidComm.Sum.Terms <- list(F1=Fid.Err.Sum.Term.1,F2=Fid.Err.Sum.Term.2,C=Comm.Err.Sum.Term)
print(str(FidComm.Terms))
print("FC.ratio")
print(str(FC.ratio))
print("FC.ratio.2")
print(str(FC.ratio.2))
print("FC.ratio.3")
print(str(FC.ratio.3))

#list(power.mc=power.mc,power.cmp=list(cca = power.cca.mc,pom = power.pom.mc,cca.reg =power.cca.reg.mc), cont.tables=cont.table,
#		config.dist= config.mismatch, min.stress=c(min.stress.for.w.val,pom.stress),means=means,FidComm.Terms=FidComm.Terms,
#		FidComm.Sum.Terms = FidComm.Sum.Terms,F.to.C.ratio = FC.ratio, wtF.to.C.ratio=FC.ratio.2
#)
list(T0=T0,TA=TA,FC1=FC.ratio,FC2=FC.ratio.2,FC3=FC.ratio.3)
} )

T0.agg[mc,,]<-Test.Stats$T0
TA.agg[mc,,]<-Test.Stats$TA
FC1.agg[mc,]<-unlist(Test.Stats$FC1)
FC2.agg[mc,]<-unlist(Test.Stats$FC2)
FC3.agg[mc,]<-unlist(Test.Stats$FC3)

}
print("Fc1 for different mc replicates ")
print(FC1.agg)
print("Fc2 for different mc replicates ")
print(FC2.agg)
print("Fc3 for different mc replicates ")
print(FC3.agg)
sink()


lty.i.vec<- rep(1,w.max.index)

windows()
for (w.ind in 1:w.max.index){
	
	avg.T0<-c()
	for ( j in 1:nmc){
		avg.T0<-c(avg.T0,T0.agg[j,w.ind,])
	}
	T0.dens <-density( avg.T0,bw=1.2)
	par( col= colors.vec[w.ind] )
	par(lty= lty.i.vec[w.ind])
	if (w.ind==1){
	plot(T0.dens,xlim=c(0,4),ylim=c(0,0.4),main="")
	} else{
	lines(T0.dens,main="")
	}
}

							legend.txt <- params$w.vals
							if (compare.pom.cca)
								legend.txt <-c(legend.txt ,"pom","cca")
							
							legend("bottomright",legend=legend.txt,
									col=colors.vec,lty=lty.i.vec)
							plot.title<-paste("T0",vary.param,unlist(params.list[[param.index]][vary.param]))
							title(plot.title)
							fname<- paste(plot.title,".pdf",collapse="")
							savePlot(filename=fname,type="pdf")

windows()							
							
for (w.ind in 1:w.max.index) {
		
	avg.T0<-c()
	for ( j in 1:nmc){
		avg.T0<-c(avg.T0,T0.agg[j,w.ind,])
	}
	T0.dens <-density( avg.T0,bw=1.2)
	par( col= colors.vec[w.ind] )
	par(lty= 1)
	
	
	avg.TA<-c()
	for ( j in 1:nmc){
		avg.TA<-c(avg.TA,TA.agg[j,w.ind,])
	}
	TA.dens <-density( avg.TA,bw=1.2)
	
	
	if (w.ind==1){
	plot(T0.dens,xlim=c(0,12),ylim=c(0,0.5),main="")
	par(lty= 2)
	lines(TA.dens)
	} else{
	lines(T0.dens)
	par(lty= 2)
	lines(TA.dens)
	}
par(lty= 1)


}
legend.txt <- params$w.vals
							if (compare.pom.cca)
								legend.txt <-c(legend.txt ,"pom","cca")
							
							legend("bottomright",legend=legend.txt,
									col=colors.vec,lty=lty.i.vec)
							plot.title<-paste("TA",vary.param,unlist(params.list[[param.index]][vary.param]))
							title(plot.title)

							fname<- paste(plot.title,".pdf",collapse="")
							savePlot(filename=fname,type="pdf")

							

windows()
ind.selected <-c(1,6,7,8)
w.vals.selected <- params$w.vals[ind.selected]
col.selected<-c("blue","red","green","brown","gold4")
par(lwd=2)
plot(x=c(0),y=c(0),xlim=c(0,12),ylim=c(0,0.5),main="",lty=1,xlab="x",ylab="Density(x)")
for (w.it in 1:length(ind.selected)) {
	w.ind<- ind.selected[w.it]
	avg.T0<-c()
	for ( j in 1:nmc){
		avg.T0<-c(avg.T0,T0.agg[j,w.ind,])
	}
	T0.dens <-density( avg.T0,bw=1)
	
	
	par(lty= 1)
	
	
	avg.TA<-c()
	for ( j in 1:nmc){
		avg.TA<-c(avg.TA,TA.agg[j,w.ind,])
	}
	TA.dens <-density( avg.TA,bw=1)
	
	
	
	
	par( col= col.selected[w.it] )
	lines(T0.dens)
	par(lty= 2)
	lines(TA.dens)
	
	

}

legend.txt <- w.vals.selected
legend.txt.2<- c("T0","TA")
lty.vec<-c(rep(1,length(ind.selected)))
							if (compare.pom.cca)
								legend.txt <-c(legend.txt ,"pom","cca")
							
							legend("topright",legend=legend.txt,
									col=col.selected,lty=lty.vec)
							legend("topleft",legend=legend.txt.2,
									col=col.selected[1],lty=c(1,2))
							#plot.title<-paste("TA",vary.param,unlist(params.list[[param.index]][vary.param]))
							#title(plot.title)
							plot.title<-"T0_and_TA_vs_w"
							fname<- paste(plot.title,".pdf",collapse="")
							savePlot(filename=fname,type="pdf")

							



							
							
							
							
save.image(paste("params",vary.param,param.index,".RData",sep="",collapse=""))
}

