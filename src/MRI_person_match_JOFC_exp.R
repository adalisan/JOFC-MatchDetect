distC <- read.table("./data/MRI_person_match_distC.txt")
distP <- read.table("./data/MRI_person_match_distP.txt")
source("./lib/simulation_math_util_fn.R")
source("./lib/simulateTestPlot.R")
source("./src/real.data.experiment.fn.R")

rep.vec <- 1:25
size.vec <- seq(0,1,0.05)
jk.res <-list()
for (rep.i in rep.vec)
 jk.res[[rep.i]]<- run.JOFC.match.jacknife.replicate (m.i = rep.i, N = 42,test.samp.size = 1, 
                                     w.val.len = 1, Diss.E=distC, Diss.F=distP,                                                                                                                                                                                                          
                                     d=5,   oos=TRUE,   separability.entries.w=TRUE
                                     , wt.equalize=FALSE, assume.matched.for.oos=FALSE
                                     , oos.use.imputed= FALSE, w.vals=c(.8), size=size.vec,
                                     verbose=TRUE, level.critical.val = 0.05)

