## Time-stamp: <wiki_classification_SA&DM.R zma 2010-09-23 00:06>

source("wiki.fn.R")
source("oosMDS.R")
load("wiki.RData")

##* manifold matching: Xi_0 == English, Xi_1 == French
##
## Problem:
##
##        Xi_0                                      Xi_1
## +-----------------+                     +------------------+
## |                 |        1 - 1        |                  |
## |     J.red       |<------------------->|       J.red      |
## |                 |  bw groups or obs   |                  |
## |                 |                     |                  |
## |        +--------+                     +--------+         |
## |        |        |   correspondence    |        |         |
## |        | J.blue |<===================>| J.blue |         |
## |        |        |   between groups    |        |         |
## +--------+--------+                     +--------+---------+
##                                              ^
##                                              |
##                                              +-- to be classified
##
## Method:
## 1. learn map psi: Xi_F --> Xi_E through Xi_E ~ Xi_F | J0
## 2. transformation: psi(Xi_F | J1)
## 3. train g1_E on Xi_E | J1
## 4. classification: g1_E( psi(Xi_F | J1) ) == g1_E o psi (Xi_F | J1)
##
## Two approaches for step 1:
## - p-approach:
##     - embed Xi_EJ0 Xi_FJ0 separately
##     - procrustes ==> Q
## - w-approach:
##     - embed Xi_EJ0 Xi_FJ0 via w-approach
##     - oosMDS Xi_EJ1 and Xi_FJ1 separately
##     - train g1_E
##     - classification: g1_E(Xi_FJ1 )

NDIM       <- 6
J.red      <- 0:2
J.blue     <- 3:4
N          <- length(label)
r.red      <- (1 :N)[label %in% J.red]
r.blue     <- (1 :N)[-r.red]
n.red      <- length(r.red)
n.blue     <- length(r.blue)
label.red  <- label[r.red]
label.blue <- label[r.blue]
loss <- "strain"


run.jofc <- function(D1, D2, D10A,D20,D2A,
					D.oos.1,
					D.oos.2.null ,
					D.oos.2.alt ,
					
					ideal.omnibus.0  ,
					ideal.omnibus.A ,
	
				n,m,
				d,c.val,
				model,oos,proc.dilation,
				verbose)   {
	
	
	w.max.index <- length(w.vals)
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
		M.oos.0 <- omnibusM(D.oos.1,D.oos.2.null, L.tilde.null)
		M.oos.A <- omnibusM(D.oos.1,D.oos.2.alt,  L.tilde.alt)
		
		
		
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
			Y1t<-Y.0t[1:m,]
			Y2t<-Y.0t[m+(1:m),]
			Y1t.A<-Y.At[1:m,]
			Y2At<-Y.At[m+(1:m),]
			

			
			
			
			T0[l,] <- rowSums((Y1t - Y2t)^2)
			TA[l,] <- rowSums((Y1t.A - Y2At)^2)
			}
			
	
	
}


return(list(T0=T0,TA=TA))

}

size <- seq(0, 1, 0.01)
oos.use.imputed<-FALSE


for (mc in 1: nmc) {
	
	
	left.out.samp<- sample(test.samp.size,1:N)
	sample.for.unmatched<- sample(test.samp.size,1:(N-test.samp.size))
	orig.indices<-1:N
	sample.for.unmatched.orig.index<-orig.indices[-left.out.samp][sample.for.unmatched]
	
	T0 <- matrix(0,w.max.index,m)   #Test statistics for JOFC under null
	TA <- matrix(0,w.max.index,m)    #Test statistics for JOFC under alternative
	
	D1<-TE[-left.out.samp,-left.out.samp]
 	D2<-TF[-left.out.samp,-left.out.samp]
 	D10A<- TE
	D20<-  TF
	
	D2A<-  omnibusM(D2,D2[sample.for.unmatched,sample.for.unmatched],D2[,sample.for.unmatched] )
	
	D.oos.1      <- TE[left.out.samp,left.out.samp]
	D.oos.2.null <- TF[left.out.samp,left.out.samp]
	D.oos.2.alt  <-  D2[sample.for.unmatched,sample.for.unmatched]
			
	L.in.oos.0 <- omnibusM(TE[-left.out.samp,left.out.samp],TF[-left.out.samp,left.out.samp],matrix(0,n,m))
	L.in.oos.A <- omnibusM(TE[-left.out.samp,left.out.samp],TF[-left.out.samp,sample.for.unmatched.orig.index],matrix(0,n,m))
					
	ideal.omnibus.0 <- omnibusM(omnibusM (D1,D2,Ltilde),omnibusM(D.oos.1,D.oos.null,L.tilde.oos),L.in.oos.0)
	ideal.omnibus.A <- omnibusM(omnibusM (D1,D2,Ltilde),omnibusM(D.oos.1,D.oos.null,L.tilde.oos),L.in.oos.A)
	
	n <- N -test.samp.size
	m <- test.samp.size

	d <- 21
	c.val <- 0
	
	power.w.star<- 0
				
	for (l in 1:w.val.len){
	
	w.val.l <- w.vals[l]			
	JOFC.results <- run.jofc(
				   D1, D2, D10A,D20,D2A,
					D.oos.1,
					D.oos.2.null ,
					D.oos.2.alt ,
					
					ideal.omnibus.0  ,
					ideal.omnibus.A ,
	
				n,m,
				d,c.val,
	
				model,oos,proc.dilation,w.val.l, oos.use.imputed, 
				verbose) 
	
			T0[l,] <- rowSums((Y1t - Y2t)^2)
			TA[l,] <- rowSums((Y1t.A - Y2At)^2)
			
			
			power.l <- get_power(T0[l,],TA[l,],size)
			if (verbose) print("JOFC test statistic complete \n")
			power.mcnemar.l <- get_power(T0[l,],TA[l,],level.mcnemar)
			if (power.mcnemar.l>power.w.star){
				rival.w <- w.vals[l]
				power.w.star <- power.mcnemar.l
				w.val.rival.idx <- l
		}
	}
	
	}
		
	
	
	
	
	}






##** out-of-sample embed X0.blue and X1.blue

## assume 1 - 1 correspondence between red objects in Xi_0 and Xi_1
pG <- matchManifolds(GE, GF, label, method="P", block=F, crit=loss)
pT <- matchManifolds(TE, TF, label, method="P", block=F, crit=loss)
wG <- matchManifolds(GE, GF, label, method="W", block=F, crit=loss, use.knn=T)
wT <- matchManifolds(TE, TF, label, method="W", block=F, crit=loss, use.knn=T)

## block: between Xi_0.red and Xi_1.red, there is no 1 - 1 correspondence,
## only the correspondence between classes in J.red is known.
pbG <- matchManifolds(GE, GF, label, method="P", block=T, crit=loss)
pbT <- matchManifolds(TE, TF, label, method="P", block=T, crit=loss)
wbG <- matchManifolds(GE, GF, label, method="W", block=T, crit=loss, use.knn=T)
wbT <- matchManifolds(TE, TF, label, method="W", block=T, crit=loss, use.knn=T)

save(list=c("pG", "pT", "wG", "wT", "pbG", "pbT", "wbG", "wbT",
            "n.red", "n.blue", "label.red", "label.blue", "N"),
     file="wiki_classification_SA&DM.Rdata")

##** LDA
## lda trained on X0.blue and tested on X1.blue (1 - 1)

sum(pG$class != label.blue)/n.blue
sum(pT$class != label.blue)/n.blue
sum(wG$class != label.blue)/n.blue
sum(wT$class != label.blue)/n.blue

## lda trained on X0.blue and tested on X1.blue (block)
sum(pbG$class != label.blue)/n.blue
sum(pbT$class != label.blue)/n.blue
sum(wbG$class != label.blue)/n.blue
sum(wbT$class != label.blue)/n.blue

##** combine G and T

## P-approach

## trained on X0.blue and tested on X1.blue (1 - 1)
z <- lda(cbind(pG$X0.blue, pT$X0.blue), label.blue)
classPGT <- predict(z, cbind(pG$X1.blue, pT$X1.blue))$class
sum(classPGT != label.blue)/n.blue

## trained on X0.blue and tested on X1.blue (block)
z <- lda(cbind(pbG$X0.blue, pbT$X0.blue), label.blue)
classPbGT <- predict(z, cbind(pbG$X1.blue, pbT$X1.blue))$class
sum(classPbGT != label.blue)/n.blue

## W-approach

## trained on X0.blue and tested on X1.blue (1 - 1)
z <- lda(cbind(wG$X0.blue, wT$X0.blue), label.blue)
classPGT <- predict(z, cbind(wG$X1.blue, wT$X1.blue))$class
sum(classPGT != label.blue)/n.blue      # 0.29308

## trained on X0.blue and tested on X1.blue (block)
z <- lda(cbind(wbG$X0.blue, wbT$X0.blue), label.blue)
classWbGT <- predict(z, cbind(wbG$X1.blue, wbT$X1.blue))$class
sum(classWbGT != label.blue)/n.blue
