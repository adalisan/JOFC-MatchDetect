## Time-stamp: <wiki_classification_SA&DM.R zma 2010-09-23 00:06>

source("./src/wiki.fn.R")
source("./src/oosMDS.R")
load("./src/wiki.RData")


library(vegan)
library(ROCR)
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

class.test <- sample (0:4,2,replace=FALSE)
c.0 <- sample (0:2,1)
c.1 <- sample (3:4,1)

c.0<- 0
c.1<- 1
class.test<-c(c.0,c.1)
testing.data.ind <- label%in%class.test
manifold.match.data.TE<-TE[!testing.data.ind,!testing.data.ind]
manifold.match.data.TF<-TF[!testing.data.ind,!testing.data.ind]

common.dim<-9
p<-12
matched.obs.classifier.obs.flag<-  testing.data.ind
class.labels<-label

test.N<- sum(matched.obs.classifier.obs.flag)
matched.obs.idx <- which(matched.obs.classifier.obs.flag)
W.choice<-"avg"
Wchoice<-"avg"
#Diss1<- manifold.match.data.TE
#Diss2<- manifold.match.data.TF 
#D1.all.obj<- TE
#D2.all.obj<- TF

#rm(manifold.match.data.TE)
#rm(manifold.match.data.TF)
#rm(TE)
#rm(TF)



matchManifold<-function(Diss1,Diss2,method,common.dim,D1.all.obj,D2.all.obj,
		class.labels,crit,Wchoice) {
#
# Diss1 and Diss2 are dissimilarity matrices (that are of classes not selected for discrimination)
#Joint Embed Diss1 and Diss2 using PoM ,CCA and JOFC

if (crit == "strain") {
	inMds  <- "cmdscale"
	oosMds <- "oosMDS"
} else {
	inMds  <- "smacofM"
	oosMds <- "smacofOos"
}

# Number of all paired observations
N<- dim(D1.all.obj)


#
# PoM
# Compute Procrustes to match the instance pairs belonging to 4 nontested classes
X.embed.1 <- cmdscale(Diss1,common.dim)
X.embed.2 <- cmdscale(Diss2,common.dim)

proc <- procrustes(X.embed.2,X.embed.1, scale=TRUE)
Q <- proc$rotation * proc$scale

#
# CCA
# Computation of "optimal" linear projection using CCA


X1.ld <- cmdscale(Diss1, p+1)		
X2.ld <- cmdscale(Diss2, p+1)

xcca <- cancor(X1.ld, X2.ld)


#
# JOFC
# Base embedding : Embedding of 4 nontested class instances

if (Wchoice == "avg") {
	W <- (Diss1+ Diss2)/2
} else if (Wchoice == "sqrt") {
	W <- sqrt((Diss1^2 + Diss2^2)/2)
}

M.base<-omnibusM(Diss1,Diss2,W)
X.embed.base <- cmdscale(M.base, common.dim)

return (list(pom=list(Q=Q,X.embed.1=X.embed.1,X.embed.2=X.embed.2),
		cca=list(xcca=xcca,X1.ld=X1.ld,X2.ld=X2.ld),jofc=X.embed.base))
}




EmbedTrainingSet<- function(D1.all.obj,D2.all.obj,BaseEmbed,
				matched.obs.classifier.obs.flag,Wchoice,method){

if (method=="all" || method=="pom")  	attach(BaseEmbed$pom)
if (method=="all" || method=="cca") 	attach(BaseEmbed$cca)
X.embed.base <- BaseEmbed$jofc
jofc.M.oos.embed<-NULL

X.train.1<-NULL
X.train.points.proc<- NULL
X.train.t.cca<- NULL
X.test.t.cca <- NULL
jofc.M.oos.embed <- NULL

	
#
# PoM classifier
#

if (method=="all" || method=="pom"){
X.train.1 <- match.fun(oosMds)(D1.all.obj,X.embed.1,w=!matched.obs.classifier.obs.flag)
X.test.points <- match.fun(oosMds)(D2.all.obj,X.embed.2,w=!matched.obs.classifier.obs.flag)

X.train.points.proc<-X.train.1%*%Q

}

#
# oos Embed for CCA and find lower dim. representation
#

if (method=="all" || method=="cca"){

X.train.ld.embed 		<- match.fun(oosMds)(D1.all.obj,X1.ld,w=!matched.obs.classifier.obs.flag)
X.test.points.ld.embed  <- match.fun(oosMds)(D2.all.obj,X2.ld,w=!matched.obs.classifier.obs.flag)

X.train.t.cca <- (X.train.ld.embed%*%xcca$xcoef)[,1:common.dim]
X.test.t.cca  <- (X.test.points.ld.embed%*%xcca$ycoef)[,1:common.dim]

}

if (method=="all" || method=="jofc"){

if (Wchoice == "avg") {
	W.M <- (D1.all.obj+ D2.all.obj)/2
} else if (Wchoice == "sqrt") {
	W.M <- sqrt((D1.all.obj^2 + D2.all.obj^2)/2)
}


M<- omnibusM(D1.all.obj,D2.all.obj,W.M)
print(sum(matched.obs.classifier.obs.flag))
print(dim(M))

gc()
jofc.M.oos.embed<-match.fun(oosMds)(M,X.embed.base,
		w=rep(!matched.obs.classifier.obs.flag,2))
}

detach(BaseEmbed$pom)
detach(BaseEmbed$cca)
#detach(jofc)

return(list(pom=list(X.train.pom=X.train.points.proc,X.test.pom=X.test.points),
		cca=list(X.train.cca=X.train.t.cca,X.test.cca=X.test.t.cca),
				jofc=jofc.M.oos.embed))
}












manif.match.new<- matchManifold(manifold.match.data.TE,manifold.match.data.TF,"all",common.dim,TE,TF,
		label,loss,W.choice)

if (loss == "strain") {
	inMds  <- "cmdscale"
	oosMds <- "oosMDS"
} else {
	inMds  <- "smacofM"
	oosMds <- "smacofOos"
}

Embedded.Train.Test.Set <- EmbedTrainingSet(TE,TF,manif.match.new,
				matched.obs.classifier.obs.flag,W.choice,"all")



#
# Classification: Testing step
#

#attach(Embedded.Train.Test.Set )

X.train.pom<-Embedded.Train.Test.Set$pom$X.train.pom
X.test.pom<-Embedded.Train.Test.Set$pom$X.test.pom
X.train.cca<-Embedded.Train.Test.Set$cca$X.train.cca
X.test.cca<-Embedded.Train.Test.Set$cca$X.test.cca
jofc<-Embedded.Train.Test.Set$jofc




true.pred.pom <- 0
true.pred.cca <- 0
true.pred.jofc <- 0

max.post.probs.pom<-rep(0,test.N)
max.post.probs.cca<-rep(0,test.N)
max.post.probs.jofc<-rep(0,test.N)
cont.table.pom<-matrix(0,2,2,dimnames=list(class.test,class.test))
cont.table.cca<-matrix(0,2,2,dimnames=list(class.test,class.test))
cont.table.jofc<-matrix(0,2,2,dimnames=list(class.test,class.test))

true.pred.pom.nt <- 0
true.pred.cca.nt <- 0
true.pred.jofc.nt <- 0

max.post.probs.pom.nt<-rep(0,test.N)
max.post.probs.cca.nt<-rep(0,test.N)
max.post.probs.jofc.nt<-rep(0,test.N)
cont.table.pom.nt<-matrix(0,2,2,dimnames=list(class.test,class.test))
cont.table.cca.nt<-matrix(0,2,2,dimnames=list(class.test,class.test))
cont.table.jofc.nt<-matrix(0,2,2,dimnames=list(class.test,class.test))





left.out.obs.vector<-1:test.N

for ( t in 1:length(left.out.obs.vector)){
	left.out.obs<- left.out.obs.vector[t]
	
	train.obs.flag <- left.out.obs.vector
	train.obs.flag<-train.obs.flag[train.obs.flag!=left.out.obs] 
	
	train.labels <- class.labels[matched.obs.idx[train.obs.flag]]
	test.label<-class.labels[matched.obs.idx[left.out.obs]]
	
	pom.classifier<-lda(X.train.pom[-left.out.obs,],train.labels)
	
	pred.class<-predict(pom.classifier,X.test.pom[left.out.obs,])
	max.post.probs.pom[t]<- pred.class$posterior[2]
	
	if (pred.class$class==test.label)
		true.pred.pom<-true.pred.pom+1
	if (pred.class$class==class.test[1]){
		if (pred.class$class==test.label){
			cont.table.pom[1,1]<-		cont.table.pom[1,1]+1
		} else{
			cont.table.pom[1,2]<-		cont.table.pom[1,2]+1
			
		}
		
	} else{
		if (pred.class$class==test.label){
			cont.table.pom[2,2]<-		cont.table.pom[2,2]+1
		} else{
			cont.table.pom[2,1]<-		cont.table.pom[2,1]+1
			
		}
		
	}
	
	
#
# CCA classifier
#
	

	
	cca.classifier<-lda(X.train.cca[-left.out.obs,],train.labels)
	pred.class<-predict(cca.classifier,X.test.cca[left.out.obs,])
	max.post.probs.cca[t]<- pred.class$posterior[2]
	if (pred.class$class==test.label)
		true.pred.cca<-true.pred.cca+1
	if (pred.class$class==class.test[1]){
		if (pred.class$class==test.label){
			cont.table.cca[1,1]<-		cont.table.cca[1,1]+1
		} else{
			cont.table.cca[1,2]<-		cont.table.cca[1,2]+1
			
		}
		
	} else{
		if (pred.class$class==test.label){
			cont.table.cca[2,2]<-		cont.table.cca[2,2]+1
		} else{
			cont.table.cca[2,1]<-		cont.table.cca[2,1]+1
			
		}
		
	}
	
	

	train.embed<-jofc[train.obs.flag,]
	test.embed <- jofc[test.N+left.out.obs,]
	jofc.classifier <- lda(train.embed,train.labels)
	pred.class <-predict(jofc.classifier,test.embed)
	max.post.probs.jofc[t]<- pred.class$posterior[2]
	if (pred.class$class==test.label)
		true.pred.jofc<-true.pred.jofc+1
	if (pred.class$class==class.test[1]){
		if (pred.class$class==test.label){
			cont.table.jofc[1,1]<-		cont.table.jofc[1,1]+1
		} else{
			cont.table.jofc[1,2]<-		cont.table.jofc[1,2]+1
			
		}
		
	} else{
		if (pred.class$class==test.label){
			cont.table.jofc[2,2]<-		cont.table.jofc[2,2]+1
		} else{
			cont.table.jofc[2,1]<-		cont.table.jofc[2,1]+1
			
		}
		
	}
	
	#Ideal POM classifier
      #nontransfer.pom.classifier<-lda(rbind(X.train.pom[-left.out.obs,]
	#						,X.test.pom[-left.out.obs,])
	#				,rep(train.labels,2))
	nontransfer.pom.classifier<-lda(X.test.pom[-left.out.obs,],
					train.labels)


	pred.class<-predict(nontransfer.pom.classifier,X.test.pom[left.out.obs,])
	max.post.probs.pom.nt[t]<- pred.class$posterior[2]
	
	if (pred.class$class==test.label)
		true.pred.pom.nt<-true.pred.pom.nt+1
	if (pred.class$class==class.test[1]){
		if (pred.class$class==test.label){
			cont.table.pom.nt[1,1]<-		cont.table.pom.nt[1,1]+1
		} else{
			cont.table.pom.nt[1,2]<-		cont.table.pom.nt[1,2]+1
			
		}
		
	} else{
		if (pred.class$class==test.label){
			cont.table.pom.nt[2,2]<-		cont.table.pom.nt[2,2]+1
		} else{
			cont.table.pom.nt[2,1]<-		cont.table.pom.nt[2,1]+1
			
		}
		
	}





      nontransfer.cca.classifier<-lda(X.test.cca[-left.out.obs,]
					,train.labels)

	pred.class<-predict(nontransfer.cca.classifier,X.test.cca[left.out.obs,])
	max.post.probs.cca.nt[t]<- pred.class$posterior[2]
	if (pred.class$class==test.label)
		true.pred.cca.nt<-true.pred.cca.nt+1
	if (pred.class$class==class.test[1]){
		if (pred.class$class==test.label){
			cont.table.cca.nt[1,1]<-		cont.table.cca.nt[1,1]+1
		} else{
			cont.table.cca.nt[1,2]<-		cont.table.cca.nt[1,2]+1
			
		}
		
	} else{
		if (pred.class$class==test.label){
			cont.table.cca.nt[2,2]<-		cont.table.cca.nt[2,2]+1
		} else{
			cont.table.cca.nt[2,1]<-		cont.table.cca.nt[2,1]+1
			
		}
		
	}
		train.TF.inst <- jofc[-c(left.out.obs.vector,test.N+left.out.obs),]
		nontransfer.jofc.classifier<-lda(train.TF.inst
					,train.labels)
	pred.class <-predict(nontransfer.jofc.classifier,test.embed)
		
	max.post.probs.jofc.nt[t]<- pred.class$posterior[2]
	if (pred.class$class==test.label)
		true.pred.jofc.nt<-true.pred.jofc.nt+1
	if (pred.class$class==class.test[1]){
		if (pred.class$class==test.label){
			cont.table.jofc.nt[1,1]<-		cont.table.jofc.nt[1,1]+1
		} else{
			cont.table.jofc.nt[1,2]<-		cont.table.jofc.nt[1,2]+1
			
		}
		
	} else{
		if (pred.class$class==test.label){
			cont.table.jofc.nt[2,2]<-		cont.table.jofc.nt[2,2]+1
		} else{
			cont.table.jofc.nt[2,1]<-		cont.table.jofc.nt[2,1]+1
			
		}
		
	}
		
}

ROCR.preds<-prediction(max.post.probs.pom,class.labels[matched.obs.idx])
ROCR.perf<- performance(ROCR.preds,'tpr','fpr')
windows()
plot(ROCR.perf,col="blue")

ROCR.preds<-prediction(max.post.probs.pom.nt,class.labels[matched.obs.idx])
ROCR.perf<- performance(ROCR.preds,'tpr','fpr')
plot(ROCR.perf,col="blue",lty=3,add=TRUE)

par(lty=1)




ROCR.preds<-prediction(max.post.probs.cca,class.labels[matched.obs.idx])
ROCR.perf<- performance(ROCR.preds,'tpr','fpr')


plot(ROCR.perf,col="red",add=TRUE)
ROCR.preds<-prediction(max.post.probs.cca.nt,class.labels[matched.obs.idx])
ROCR.perf<- performance(ROCR.preds,'tpr','fpr')
plot(ROCR.perf,col="red",lty=3,add=TRUE)


par(lty=1)


ROCR.preds<-prediction(max.post.probs.jofc,class.labels[matched.obs.idx])
ROCR.perf<- performance(ROCR.preds,'tpr','fpr')

plot(ROCR.perf,col="green",add=TRUE)

ROCR.preds<-prediction(max.post.probs.jofc.nt,class.labels[matched.obs.idx])
ROCR.perf<- performance(ROCR.preds,'tpr','fpr')
plot(ROCR.perf,col="green",lty=3,add=TRUE)
legend.txt<-c("pom","cca","jofc")
class.type<- c("With Domain Transfer","With Orig. Domain Data")
legend.txt <-paste(rep(legend.txt,each=2),class.type)

legend("bottomright",legend=legend.txt,
	col=rep(c("blue","red","green"),each=2),
	lty=c(1,3))


