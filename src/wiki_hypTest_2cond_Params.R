print(getwd())

debug.mode <- FALSE

source("./lib/graph_embedding_fn.R")

diss_measure <- "default" #shortest path

use.text.diss <- TRUE

if (use.text.diss){
  load("./data/wiki.RData")
  
  Diss.E<- TE
  Diss.F<- TF
} else if (diss_measure=="default"){
    load("./data/wiki.RData")
   
    Diss.E<- GE
    Diss.F<- GF
  } else {
  load("./data/wiki_adj.RData")
  cached.diss.file <- "./cache/wiki_Dice.RData"
  if (file.exists(cached.diss.file))
    load(cached.diss.file)
  else{
    GE_dice <- C_dice_weighted(en_a)
    GF_dice <- C_dice_weighted(fr_a)
    Diss.E = GE_dice
    Diss.F = GF_dice
    if ("cache" %in% dir())
      save(Diss.E,Diss.F,file="./cache/wiki_Dice.RData")
  }
}





w.vals = c(0.1,0.2,0.4,0.5,0.8,0.925,0.95,0.99,0.999)
nmc<-4
test.samp.size <- 200

par.compute <- TRUE

if (debug.mode){
  nmc <- 4
  w.vals <- c(0.1,0.8,0.9,0.99)
  Diss.E<- Diss.E[1:550,1:550]
  Diss.F<- Diss.F[1:550,1:550]
  test.samp.size <- 100
par.compute<- FALSE
}
w.val.len <- length(w.vals)

N <- dim(Diss.E)[1]






d <- 12
size <- seq(0, 1, 0.01)
#oos.use.imputed<-FALSE
oos.use.imputed<-TRUE
level.mcnemar <- 0.02
power.nmc <-array(0,dim=c(w.val.len,nmc,length(size)))

m<- test.samp.size
n<- N-(2*test.samp.size)

Wchoice<-"avg"
size <- seq(0, 1, 0.01)

#originally was false
separability.entries.w<- TRUE

wt.equalize<-FALSE
assume.matched.for.oos<-TRUE
oos.use.imputed<-FALSE
oos <- TRUE
verbose <-TRUE
pre.scaling<-FALSE


#prescaling
if (pre.scaling) {
  s <- lm(as.vector(Diss.E) ~ as.vector(Diss.F) + 0)$coefficients
} else {
  s <- 1
}

Diss.F<-Diss.F*s

