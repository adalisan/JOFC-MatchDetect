
run.experiment.JOFC<-function(G,Gp,n_vals,num_iter,embed.dim,diss_measure="default"){
  
  library(optmatch)
  library(igraph)
  library(MASS)
  library(MCMCpack)
  source("./lib/graph_embedding_fn.R")
  source("./lib/simulation_math_util_fn.R")
  source("./lib/smacofM.R")
  source("./lib/oosIM.R")
  source("./lib/diffusion_distance.R")
  matched.cost<-0.01
  N<-nrow(G)
  corr.matches =matrix(0,length(n_vals),num_iter)
  
  for (n_v_i in 1:length(n_vals)){
    n_v = n_vals[n_v_i]
    for (it in 1:num_iter){
      
      #insample_logic_vec <- 1:N %in% 1:n_v
      insample_logic_vec <- 1:N %in% sample(1:N,n_v,replace=FALSE)
      print(insample_logic_vec)
      insample_logic_vec <- c(insample_logic_vec,insample_logic_vec)
      jofc.result<- JOFC.graph.custom.dist(G,Gp,in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=0.95,graph.is.directed=FALSE, vert_diss_measure=diss_measure)
      #jofc.result<-try(jofc(G,Gp, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=0.9,graph.is.directed=FALSE, oos=TRUE,use.weighted.graph=FALSE))
      #if (inherits(jofc.result,"try-error")) {
      #	print('Skipping iteration')
      #		next}
      
      jofc.res.1<-jofc.result[[1]]
      
      M.result.1<-try(solveMarriage(jofc.res.1))
      
      
      if (inherits(M.result.1,"try-error"))    {  
        print('Skipping iteration')
        next}
      
      NumofTruePairing.1<-present(M.result.1)
      corr.matches[n_v_i,it] = NumofTruePairing.1
    }
  }
  save(file.path("logs",
                 paste("JOFC_graph",Sys.time(),runif(n=1,max=100))),
       list=c(corr.matches))
  return(corr.matches)
}

bitflip_MC_rep <- function (pert,n,n_vals,embed.dim,diss_measure){
  corr.match.array.mc<-array(0,dim=c(length(n_vals),npert))
 
  for(ipert in 1:npert)
  {
    G<-ER(n,0.5)
    Y.emb<-NULL
    Gp<-bitflip(G ,pert[ipert],pert[ipert])
    corr.matches<-run.experiment.JOFC(G,Gp,n_vals,num_iter=1, embed.dim=embed.dim,diss_measure=diss_measure)
    corr.match.array.mc[,ipert] <- corr.matches
  }
  return (corr.match.array.mc)
}

bitflip_exp<-function (nmc,pert,n,n_vals,embed.dim=6)
{
  npert<-length(pert)
  corr.match.array<-array(0,dim=c(length(n_vals),nmc,npert))
  corr.match.avg<-array(0,dim=c(length(n_vals),npert))
  
  seed<-123
  set.seed(seed)
  sfInit( parallel=TRUE)
  print(sfCpus())
  corr.match.list <- sfLApply(1:nmc,bitflip_MC_rep,pert,n,n_vals,embed.dim)
  
  sfStop()
  for t in 1:length(corr.match.list)
      corr.match.array[,t,] <- corr.match.list[[t]]
     
  
  corr.match.array<-sweep(corr.match.array, 1, n-n_vals, "/")
  print(corr.match.array)
  print(dim(corr.match.array))
  colors.vec<-c( "red","blue","orange","green")
  for (i in 1:length(n_vals)) 
    for (j in 1:npert){
      corr.match.avg[i,j] <- mean(corr.match.array[i,,j])
    }
  
  #windows()
  plot(n_vals, as.vector(corr.match.avg[,1]) ,xlab="Hard seeds",
       ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[1],type="l")
  
  for(ipert in 2:npert)
  {
    lines(n_vals, as.vector(corr.match.avg[,ipert]) ,xlab="Hard seeds",ylab="Fraction of  correct matches",ylim=c(0,1),col=colors.vec[ipert])
  }  
  
  return(corr.match.avg)
}

wiki_exp <- function(num_iter,n_vals,embed.dim=13) {
  load ("./data/wiki.RData")
  corr.matches<-run.experiment.JOFC(GE,GF,n_vals,num_iter=num_iter,embed.dim,diss_measure="default")
  
}	

worm_exp <- function(num_iter,n_vals,embed.dim=6) {
  load("./data/celegansGraph.Rd")
  corr.matches<-run.experiment.JOFC(Ac,Ag,n_vals,num_iter=num_iter,embed.dim,diss_measure="default")
}  


enron_exp <- function (num_iter,n_vals_nec,embed.dim=6){
  load("./data/AAA-187As-184x184.Rbin")
  corr.matches<-run.experiment.JOFC(AAA[[130]],AAA[[131]],n_vals_vec,num_iter=num_iter,embed.dim,diss_measure="default")
  
}

