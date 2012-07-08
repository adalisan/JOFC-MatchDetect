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
	return(corr.matches)
}
			
bitflip_exp<-function (nmc,pert,n,n_vals,embed.dim=6)
{
	npert<-length(pert)
	corr.match.avg<-array(0,dim=c(length(n_vals),nmc,npert))
	
	seed<-123
set.seed(seed)
for(imc in 1:nmc)
{
	G<-ER(n,0.5)
	for(ipert in 1:npert)
	{
		Y.emb<-NULL
		Gp<-bitflip(G ,pert[ipert],pert[ipert])
		corr.matches<-run.experiment.JOFC(G,Gp,n_vals,num_iter=15,embed.dim,diss_measure="default")
		corr.match.avg[,imc,ipert] <- rowMeans(corr.matches)
	}
	
}

return(corr.match.avg)
}
		
		