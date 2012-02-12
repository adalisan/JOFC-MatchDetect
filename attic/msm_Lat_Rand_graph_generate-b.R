require(msm)

K<-3 #Num of states

num.V<- 100 #Num. of Vertices

T<-10 #Determines time from start to end for the time series of graphs
Trans.Mat<-rbind(c(0.5,0.1,0.1),
				c(0.1,0,  0.2),
				c(0.1,0.01, 0))


Traj.Ens<-sim.msm(q=Trans.Mat,maxtime=T)

gen.Graph.Pair<- function(Traj.Ens,
							inst.ensemble #The index of the trajectory  
					){
  P.0<-matrix(0,num.V,num.V)
  P.1<-matrix(0,num.V,num.V)

  P.vk<-matrix(0,num.V,K)



  for (v in 1:num.V){



    Traj.Ens<-sim.msm(q=Trans.Mat,maxtime=T)

    num.transition <- length(Traj.Ens$times)-1
    state.dur<-Traj.Ens$times[2:(num.transition+1)] -Traj.Ens$times[1:num.transition]
    for (k in 1:K){
	    times.in.k<-which(Traj.Ens$states==k)
	    total.time.in.k<-sum(state.dur[times.in.k])
	 
	    P.vk[v,k]<- total.time.in.k/T
	  }
  }

  for (i in 1:num.V)
	  for (j in i:num.V)

		  P.1[i,j] <- P.1[j,i]<-P.vk[i,]%.%P.vk[j,]
}

