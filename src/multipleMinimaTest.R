
require(gdata)
require(Hmisc)
require(rgl)
require(animation)
source("./lib/smacofM.R")
source("./lib/oosIM.R")
verbose <- FALSE
run.in.linux<-FALSE
results.dir<-"./graphs"
plot.in.3d<-TRUE
create.ani <-FALSE

fps=10

ani.options(outdir=file.path(Sys.getenv("PROJECT_DIR"),"graphs"))

meshgrid <- function(a,b) {
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
} 
raw.stress.at <-function(config){
  sum(as.dist(W)*((dist(config)-d.X)^2))
}

hessian.mat <- function (X.embed.2.norm,n) {
  hess.mat.size<-4*n^2
  hess.mat<-matrix(0,2*n,2*n)
  for (i.r in 1:n){
    for (i.d in 1:2){
      hess.idx.r<-(2*(i.r-1)+i.d)
      X.embed.2.norm.minus <- X.embed.2.norm.plus<-X.embed.2.norm
      X.embed.2.norm.plus[i.r,i.d] <- X.embed.2.norm.plus[i.r,i.d] + epsilon
      X.embed.2.norm.minus[i.r,i.d] <- X.embed.2.norm.minus[i.r,i.d] - epsilon        
      dir.deriv[w.i,i,j,]<-(raw.stress.at(X.embed.2.norm.plus)- raw.stress.at(X.embed.2.norm.minus))/(2*epsilon)
      for (j.r in i.r:n){
        for (j.d in 1:2){
          hess.idx.c<-(2*(j.r-1)+j.d)
          X.embed.2.norm.back.forw<-X.embed.2.norm.forw.forw <- X.embed.2.norm.forw.back <-  X.embed.2.norm.back.back <-X.embed.2.norm
          X.embed.2.norm.forw.forw[i.r,i.d] <- X.embed.2.norm.forw.forw[i.r,i.d] + epsilon
          X.embed.2.norm.forw.forw[j.r,j.d] <- X.embed.2.norm.forw.forw[j.r,j.d] + epsilon
          
          X.embed.2.norm.forw.back[i.r,i.d] <- X.embed.2.norm.forw.back[i.r,i.d] + epsilon
          X.embed.2.norm.forw.back[j.r,j.d] <- X.embed.2.norm.forw.back[j.r,j.d] - epsilon
          
          X.embed.2.norm.back.back[i.r,i.d] <- X.embed.2.norm.back.back[i.r,i.d] - epsilon
          X.embed.2.norm.back.back[j.r,j.d] <- X.embed.2.norm.back.back[j.r,j.d] - epsilon
          
          X.embed.2.norm.back.forw[i.r,i.d] <- X.embed.2.norm.back.forw[i.r,i.d] - epsilon
          X.embed.2.norm.back.forw[j.r,j.d] <- X.embed.2.norm.back.forw[j.r,j.d] + epsilon
          
          # Approximate entry of hessian matrix by finite difference
          hess.mat[hess.idx.r,hess.idx.c]<- (raw.stress.at(X.embed.2.norm.forw.forw) - raw.stress.at(X.embed.2.norm.forw.back)
                                             - raw.stress.at(X.embed.2.norm.back.forw) + raw.stress.at(X.embed.2.norm.back.back))/(4*epsilon*epsilon)
        }
      }
      
    }
  }
  return(hess.mat)
}

stress.plot3d <- function (time,sign.hessian.at.pt,x.coords,y.coords,
                           stress.at.loc.w,grid.seq.x,grid.seq.y,w.vals,
                           rotate.z.angle=0) {
  w.i <- min(floor(time/1.5)+1,length(w.vals))
  print('time and w.i')
  print(time)
  print(w.i)
  
  
  col.matrix = sign.hessian.at.pt[w.i,,]
  col.matrix[sign.hessian.at.pt[w.i,,]==1]="orange"
  col.matrix[sign.hessian.at.pt[w.i,,]==2]="red"
  col.matrix[sign.hessian.at.pt[w.i,,]==0]="green"
  col.matrix[sign.hessian.at.pt[w.i,,]==-1]="blue"
  col.matrix[sign.hessian.at.pt[w.i,,]==-2]="pink"
  clear3d(type="shapes")
  
  plot3d(x=x.coords, y=y.coords,
         #plot3d(unmatrix(mesh.grid.coords$x,byrow=FALSE),
         #   unmatrix(mesh.grid.coords$y,byrow=FALSE),[]
         z=stress.at.loc.w[w.i,,],col=unmatrix(col.matrix),byrow=FALSE,
         xlab="x",ylab="y",zlab="Stress" ,box=FALSE,axes=FALSE,
         xlim=c(min(x.coords),max(x.coords)),ylim=c(min(y.coords),max(y.coords)),
         zlim=c(0,max(stress.at.loc.w)))
  decorate3d(xlim=c(min(x.coords),max(x.coords)),ylim=c(min(y.coords),max(y.coords)),
             zlim=c(0,max(stress.at.loc.w)))
  surface3d(grid.seq.x,grid.seq.y,stress.at.loc.w[w.i,,])
  title3d()
  title3d(paste("Stress w=",eval(expression(w.vals[w.i]))),col="red")
  if (time==0 || w.i==1 && time<(2.0/fps)){
    new_mat=transform3d(par3d("userMatrix"),angle=rotate.z.angle,x=0,y=0,z=1)
  } else{
    new_mat=par3d("userMatrix")
  }
  return(list(userMatrix=new_mat ))
}


animate.config.w <- function () {
  
  oopt = ani.options(interval = 1, nmax = length(w.vals))
  ## use a loop to create images one by one
  for (w.i in 1:length(w.vals)) {
    par(pch=1)
    if (is.null(far.to.init.X5.for.w[[w.i]]) || nrow(far.to.init.X5.for.w[[w.i]])>0)
      plot(x=far.to.init.X5.for.w[[w.i]][,1],y=far.to.init.X5.for.w[[w.i]][,2],col="red",
           xlim=c(min(grid.seq.x),max(grid.seq.x))
           ,ylim=c(min(grid.seq.y),max(grid.seq.y)),xlab="x",ylab="y")
    else{
      
      plot.new()
      plot.window(xlim=c(min(grid.seq.x),max(grid.seq.x))
                  ,ylim=c(min(grid.seq.y),max(grid.seq.y)))
      
      axis(1)
      axis(2)
      
      title(xlab="x",ylab="y")
      
      box()
    }
    
    
    par(pch=3)
    if (is.null(close.to.init.X5.for.w[[w.i]]) || nrow(close.to.init.X5.for.w[[w.i]])>0)
      points(x=close.to.init.X5.for.w[[w.i]][,1],y=close.to.init.X5.for.w[[w.i]][,2],col="red",
             xlim=c(min(grid.seq.x),max(grid.seq.x))
             ,ylim=c(min(grid.seq.y),max(grid.seq.y)),xlab="x",ylab="y")
    par(pch=1)
    if (is.null(far.to.init.X6.for.w[[w.i]]) || (nrow(far.to.init.X6.for.w[[w.i]])>0))
      points(x=far.to.init.X6.for.w[[w.i]][,1],y=far.to.init.X6.for.w[[w.i]][,2],col="blue",
             xlim=c(min(grid.seq.x),max(grid.seq.x))
             ,ylim=c(min(grid.seq.y),max(grid.seq.y)),xlab="x",ylab="y")
    par(pch=3)
    if (is.null(close.to.init.X6.for.w[[w.i]]) || (nrow(close.to.init.X6.for.w[[w.i]])>0))
      points(x=close.to.init.X6.for.w[[w.i]][,1],y=close.to.init.X6.for.w[[w.i]][,2],col="blue",
             xlim=c(min(grid.seq.x),max(grid.seq.x))
             ,ylim=c(min(grid.seq.y),max(grid.seq.y)),xlab="x",ylab="y")
    title(paste("Final point config. of X_5 and X_6 for w=",eval(expression(w.vals[w.i]))))
    
    ani.pause() ## pause for a while ('interval')
  }
  ## restore the options
  ani.options(oopt)
  
}


animate.final.loc.stress.w <- function () {
  
  oopt = ani.options(interval = 1, nmax = length(w.vals))
  ## use a loop to create images one by one
  for (w.i in 1:length(w.vals)) {
    par(pch=1)
    stress.vals<-as.vector(stress.at.loc.w[w.i,,])
    
    
    if (((max(stress.vals))-min(stress.vals))>5E-3){
      value.intervals <- hist(stress.vals,breaks=3,
                              plot=FALSE)
      print(levels(value.intervals))
      pt.at.level.1 <- (stress.vals<=value.intervals$breaks[2])
      pt.at.level.2 <- (stress.vals>value.intervals$breaks[2]) & (stress.vals<=value.intervals$breaks[3])
      pt.at.level.3 <- (stress.vals>value.intervals$breaks[3])
      
    }else{
      #value.intervals==levels(value.intervals)[1]
      pt.at.level.1 = 1:length(stress.vals)
      pt.at.level.2 = c()
      pt.at.level.3 = c()
      
    }
    
    plot  (x=final.coords.x.5.w[w.i,,][ pt.at.level.1],
           y=final.coords.y.5.w[w.i,,][ pt.at.level.1],col="red",
           xlim=c(min(grid.seq.x),max(grid.seq.x)),
           ylim=c(min(grid.seq.y),max(grid.seq.y)),
           xlab="x",ylab="y"
    )
    points(x=final.coords.x.5.w[w.i,,][ pt.at.level.2],
           y=final.coords.y.5.w[w.i,,][ pt.at.level.2],col="purple")
    
    points  (x=final.coords.x.5.w[w.i,,][ pt.at.level.3],
             y=final.coords.y.5.w[w.i,,][ pt.at.level.3],col="blue")        
    
    par(pch=3)
    points  (x=final.coords.x.6.w[w.i,,][ pt.at.level.1],
             y=final.coords.y.6.w[w.i,,][ pt.at.level.1],col="red",
             xlim=c(min(grid.seq.x),max(grid.seq.x)),
             ylim=c(min(grid.seq.y),max(grid.seq.y))
    )
    points(x=final.coords.x.6.w[w.i,,][ pt.at.level.2],
           y=final.coords.y.6.w[w.i,,][ pt.at.level.2],col="purple")
    
    points  (x=final.coords.x.6.w[w.i,,][ pt.at.level.3],
             y=final.coords.y.6.w[w.i,,][ pt.at.level.3],col="blue")
    X[5,]<-c(1,0)
    X[6,]<-c(0,1)
    points(x=c(1,0),y=c(0,1),col="black",pch=c(1,3))
    
    legends.txt<-c("Lowest Stress","Medium","Highest","True Coords","X_5","X_6")
    legend(legend=legends.txt,x="topright",col=c("red","purple","blue","black","red","red"),
           pch=c(1,1,1,1,1,3))
    
    ## restore the options
    title(paste("Final point config. of X_5 and X_6 for w=",eval(expression(w.vals[w.i]))))
    
    title(eval(expression(w.vals[w.i])))
    ani.pause() ## pause for a while ('interval')
  }
  
  
  ani.options(oopt)
  
}



#
# Create point configuration
#


n<-7
X<-matrix(0,n,2)
X[,1]<-c(0,1,1,0,0,0,0.5)#1.01)
X[,2]<-c(0,0,1,1,0,0,0.5)#,0)
X[5,]<-c(1,0)
X[6,]<-c(0,1)
#X[8,]<-c(0,1.0)
in.sample.ind<-c(1:4,7)#,8)
d.X<- dist(X)

epsilon <-1E-3
w <- 0.99

#
# Perturb Dissimilarity matrix
#
i <- 4
j <- 5

ind <- (n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<-d.X[ind]-1.4

i<-2
j<-6
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

d.X[ind]<- d.X[ind]-1.4
#X[5,]<-c(0,1)


i<-5
j<-7
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

#d.X[ind]<- d.X[ind]+0.707

i<-5
j<-8
ind<-(n*(i-1) - i*(i-1)/2 + j-i)

#d.X[ind]<- d.X[ind]-1.4

if (verbose) print(d.X)

init.config<-X

w.vals<-c(0.1,seq(0.2,0.4,0.1),seq(0.41,0.54,0.01),seq(0.55,0.7,0.05),seq(0.75,0.85,0.01),seq(0.90,0.95,0.05),0.99)
w.vals.sp<-w.vals[c(1,13,14,15,21,26:31,35,36)]
grid.seq.x<-seq(-0.5,1.5,0.1)
grid.seq.y<-seq(-0.2,1.6,0.1)

#
# Define empty arrays for results
#
final.close.to.init.w<-rep(0, length(w.vals))

dir.deriv<-array(0,dim=c(length(w.vals),length(grid.seq.x),length(grid.seq.y),2*n))
hessian.at.pt<-array(list(c(1,1,2,1)),dim=c(length(w.vals),length(grid.seq.x),length(grid.seq.y)))
sign.hessian.at.pt<-array(0,dim=c(length(w.vals),length(grid.seq.x),length(grid.seq.y)))
stress.at.loc.w <-sign.hessian.at.pt 

min.config.stress.1.w<- rep(0, length(w.vals))
min.config.stress.2.w<- rep(0, length(w.vals))


close.to.init.X5.for.w <- list()
close.to.init.X6.for.w <- list()
far.to.init.X5.for.w   <- list()
far.to.init.X6.for.w   <- list()



final.coords.x.5.w<-array(0,dim=c(length(w.vals),length(grid.seq.x),length(grid.seq.y)))
final.coords.y.5.w<-array(0,dim=c(length(w.vals),length(grid.seq.x),length(grid.seq.y)))
final.coords.x.6.w<-array(0,dim=c(length(w.vals),length(grid.seq.x),length(grid.seq.y)))
final.coords.y.6.w<-array(0,dim=c(length(w.vals),length(grid.seq.x),length(grid.seq.y)))



for (w.i in 1:length(w.vals)){
  
  stress.at.loc<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
  
  #Set up weight vector
  w <- w.vals[w.i]
  W <- matrix(1-w,n,n)
  W[5,2] <- W[2,5]<- w
  W[6,4] <- W[4,6]<- w
  W[5,7] <- W[7,5]<- 1-w
  diag(W)<-0
  
  
  W.oos <- W
  W.oos[in.sample.ind, in.sample.ind]<-0	
  new.index.order<- c(in.sample.ind,c(5,6))
  W.oos<-W.oos[new.index.order,new.index.order]
  
  
  
  in.sample.Bool<-(1:n %in% in.sample.ind)
  X.embed.1.in <- smacofM(as.matrix(d.X)[in.sample.ind,in.sample.ind],
                          ndim    = 2,
                          W       = W[in.sample.ind,in.sample.ind],
                          init    = init.config[in.sample.Bool,]		,
                          verbose = FALSE,
                          itmax   = 1000,
                          eps     = 1e-6)
  
  X.embed.1.oos <- oosIM(D=as.matrix(d.X),
                         X=X.embed.1.in,
                         init= init.config[!in.sample.Bool,],
                         W=W.oos,
                         verbose = FALSE,
                         itmax   = 1000,
                         eps     = 1e-6,
                         bwOos = TRUE,
                         isWithin = ifelse(in.sample.Bool,1,0) )
  
  row.names(X.embed.1.oos)<-NULL
  row.names(X.embed.1.in)<-NULL
  X.embed.1.norm <- rbind(X.embed.1.in[1:4,] ,X.embed.1.oos,X.embed.1.in[5:(n-2),]) 
  if (verbose) print("Init config ")
  if (verbose) print(X.embed.1.norm)
  stress <- raw.stress.at(X.embed.1.norm)
  if (verbose) print("Init config  stress")
  if (verbose)  print(stress)
  close.to.init.1<-data.frame(x=numeric(),y=numeric())
  close.to.init.2<-data.frame(x=numeric(),y=numeric())
  far.to.init.1<-data.frame(x=numeric(),y=numeric())
  far.to.init.2<-data.frame(x=numeric(),y=numeric())
  min.stress.1<-100
  min.stress.2<-100
  
  
  
  grid.resp<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
  grid.coords<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
  
  
  mesh.grid.coords<-meshgrid(grid.seq.x,grid.seq.y)
  
  
  final.coords.x.5<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
  final.coords.y.5<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
  final.coords.x.6<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
  final.coords.y.6<-array(0,dim=c(length(grid.seq.x),length(grid.seq.y)))
  min.config.1<-matrix(0,n,2)
  min.config.2<-matrix(0,n,2)
  
  for (j in 1:length(grid.seq.y)){
    
    for (i in 1:length(grid.seq.x)){
      init.config <- X
      
      i.x<-grid.seq.x[i]
      j.y<-grid.seq.y[j]
      #i.x<-mesh.grid.coords$x[i,j] 
      #j.y<-mesh.grid.coords$y[i,j] 
      
      init.config[5,]<-c(i.x,j.y)
      init.config[6,]<-c(1-i.x,1-j.y)
      
      X.embed.2.in<-smacofM(as.matrix(d.X)[in.sample.ind, in.sample.ind],
                            ndim    = 2,
                            W       = W[in.sample.ind,in.sample.ind],
                            init    = init.config[in.sample.Bool,],
                            verbose = FALSE,
                            itmax   = 1000,
                            eps     = 1e-6)
      
      X.embed.2.oos<-oosIM(D = as.matrix(d.X),
                           X = X.embed.2.in,
                           init = init.config[!in.sample.Bool,],
                           
                           verbose =  FALSE,#(i<=(length(grid.seq.x)/2)),
                           itmax   = 1000,
                           eps     = 1e-6,
                           W = W.oos,
                           bwOos = TRUE,
                           isWithin = ifelse(in.sample.Bool,1,0))
      
      
      
      
      row.names(X.embed.2.oos)<-NULL
      row.names(X.embed.2.in)<-NULL
      
      X.embed.2.norm<- rbind(X.embed.2.in[1:4,] ,X.embed.2.oos,X.embed.2.in[5:(n-2),])		
      
      hess.mat<- hessian.mat(X.embed.2.norm,n)
      hess.mat[hess.mat==0]<-t(hess.mat)[hess.mat==0]
      
      hessian.at.pt[[w.i,i,j]]<- as.list(hess.mat)
      
      hess.mat<-hess.mat[9:12,9:12]
      hess.eig<- eigen(hess.mat,symmetric=TRUE,only.values=TRUE)
      e.vals<- hess.eig$values
      low.than.thres<- abs(e.vals)<1E-5
      e.vals[low.than.thres] <-  0 #sign(e.vals[low.than.thres])*1E-5 
      
      if (sum(e.vals<0)==0){
        if (sum((e.vals==0)>0)){
          sign.hessian.at.pt[w.i,i,j] <- 2
        } else{
          sign.hessian.at.pt[w.i,i,j] <- 1 #pos definite
        }
      } else if (sum(e.vals>0)==0) {
        if  (sum((e.vals==0)>0)){
          sign.hessian.at.pt[w.i,i,j] <- -2 #neg definite
        } else{
          sign.hessian.at.pt[w.i,i,j] <- -1 #neg definite
          
        }
        
      } else{
        sign.hessian.at.pt[w.i,i,j] <- 0  #neither pos nor neg definite  saddle point
      }
      final.coords.x.5[i,j]<-X.embed.2.norm[5,1]
      final.coords.y.5[i,j]<-X.embed.2.norm[5,2]
      final.coords.x.6[i,j]<-X.embed.2.norm[6,1]
      final.coords.y.6[i,j]<-X.embed.2.norm[6,2]
      #print("X.embed.2.norm")
      #print(X.embed.2.norm)
      stress <- sum(as.dist(W)*((dist(X.embed.2.norm)-d.X)^2))
      stress.unif.wt <- sum(((dist(X.embed.2.norm)-d.X)^2))
      stress.unif.abs <- sum(abs(dist(X.embed.2.norm)-d.X))
      #print(stress)
      #	print(stress.unif.wt)
      #	print(stress.unif.abs)		
      if (((X.embed.2.norm[6,2]>0.5)) & 
        (((X.embed.2.norm[5,1]>X.embed.2.norm[6,1])))) { 
        #& (sum(X.embed.2.norm[6,]<0.3)==1))
        
        
        grid.resp[i,j]<-1
        #	print("First Test Point")
        #	print(X.embed.2.norm[5,])
        if (stress < min.stress.1) {
          min.stress.1<-stress
          if (verbose)	print("Min stress found(real min)")
          if (verbose) print(X.embed.2.norm)
          min.config.1<-X.embed.2.norm
          if (verbose) print(min.stress.1)
          
        }
        close.to.init.1<-rbind(close.to.init.1,X.embed.2.norm[5,])
        close.to.init.2<-rbind(close.to.init.2,X.embed.2.norm[6,])
        
        
      } else{
        grid.resp[i,j]<-0
        #print("First Test Point")
        if (verbose) print(X.embed.2.norm[5,])		
        if (stress < min.stress.2){
          min.stress.2<-stress
          if (verbose) print("Min stress found(second min)")
          if (verbose) print(X.embed.2.norm)
          if (verbose) 	print(min.stress.2)
          min.config.2<-X.embed.2.norm
          
        }
        far.to.init.1<-rbind(far.to.init.1,X.embed.2.norm[5,])
        far.to.init.2<-rbind(far.to.init.2,X.embed.2.norm[6,])
        
        #} else{
        #	grid.resp[i,j] <- NA
      }
      
      if (w.i>9) {
        print("stress")
        print(stress)
        print(as.dist(W))
      }
      stress.at.loc[i,j]<- stress
      
      
      
      
    }
    
  }
  
  stress.at.loc.w[w.i,,] <-stress.at.loc
  final.coords.x.5.w[w.i,,]=final.coords.x.5
  final.coords.x.6.w[w.i,,]=final.coords.x.6
  final.coords.y.5.w[w.i,,]=final.coords.y.5
  final.coords.y.6.w[w.i,,]=final.coords.y.6
  
  #plot(x=grid.seq.x,y=grid.seq.y, col=grid.resp)
  #x.coords = [x1 x2 x3 x4 ... x1 x2 x3 x4 ...
  x.coords <- rep(grid.seq.x,length(grid.seq.y))
  #y.coords = [y1 y1 y1 y1 ... y2 y2 y2 y2 ...
  y.coords <- rep(grid.seq.y,each=length(grid.seq.x))
  #grid.resp<-grid.resp[length(grid.seq.x)A:1,]
  
  if (w %in% w.vals.sp){
    if (run.in.linux) {X11()} else {windows()}
    plot(x.coords, y.coords,
         #plot(unmatrix(mesh.grid.coords$x,byrow=FALSE),
         #     unmatrix(mesh.grid.coords$y,byrow=FALSE),
         col=ifelse(unmatrix(grid.resp,byrow=FALSE)==1,"red","black"))
  }
  print("For w.value")
  print(w)
  print("Min stress found(real min)")
  if (verbose) print(min.config.1)
  print(min.stress.1)
  print("Min stress found(second min)")
  if (verbose) print(min.config.2)
  print(min.stress.2)
  min.config.stress.1.w[w.i]<-min.stress.1
  min.config.stress.2.w[w.i]<-min.stress.2
  
  
  if(!is.vector(close.to.init.1)){
    if (run.in.linux) {X11()} else {windows()}
    
    par(pch=1)
    plot(x=close.to.init.1[,1],y=close.to.init.1[,2],col="red",
         xlim=c(min(close.to.init.1[,1],close.to.init.2[,1]),max(close.to.init.1[,1],close.to.init.2[,1]))
         ,ylim=c(min(close.to.init.1[,2],close.to.init.2[,2]),max(close.to.init.1[,2],close.to.init.2[,2])))
    
    par(pch=3)
    points(x=close.to.init.2[,1],y=close.to.init.2[,2],col="blue")
    title(paste("Final config Close to true config- w=",w,collapse=""))
    legend("topright",legend=c(expression(X[6]),expression(X[7])),col=c("red","blue"),pch=c(1,3))
    
    dev.print(paste(results.dir,"/","true-min-w",w,".png",collapse="",sep=""),device=png,width=600,height=600)
    fname<-paste(results.dir,"/","true-min-w",w,".pdf",collapse="",sep="")
    dev.copy2pdf(file=fname)
    dev.off()
    
    #select.x<- sort( sample.int(length(grid.seq.x) , 10))
    select.x<-1:length(grid.seq.x) 
    #select.y <- sort( sample.int(length(grid.seq.y) , 10))
    select.y<-1:length(grid.seq.y) 
    #The indexing (select.x,select.y) is mixed because mesh.grid function generates a matrix 
    #whose columns are for x coordinates, while for final.coords, rows are for x coordinates
    
  }
  
  if(!is.vector(far.to.init.1)){
    if (run.in.linux) {X11()} else {windows()}
    par(pch=1)
    plot(x=far.to.init.1[,1],y=far.to.init.1[,2],col="red",
         xlim=c(min(grid.seq.x),max(grid.seq.x))
         ,ylim=c(min(grid.seq.y),max(grid.seq.y)))
    par(pch=3)
    points(x=far.to.init.2[,1],y=far.to.init.2[,2],col="blue")
    title(paste("Final config Far to true config- w=",w,collapse=""))
    legend("topright",legend=c(expression(X[6]),expression(X[7])),col=c("red","blue"),pch=c(1,3))
    
    dev.print(paste(results.dir,"/","other-min-w",w,".png",collapse="",sep=""),device=png, width=600,height=600)
    fname<-paste(results.dir,"/","other-min-w",w,".pdf",collapse="",sep="")
    dev.copy2pdf(file=fname)
    
  }
  
  print("Number of final config close to initial config")
  print(dim(close.to.init.1)[1])
  final.close.to.init.w[w.i] <- dim(close.to.init.1)[1]
  
  print("Number of final config far to initial config")
  print(dim(far.to.init.1)[1])
  min.config.stress.1.w[w.i]<- min.stress.1
  min.config.stress.2.w[w.i]<- min.stress.2
  
  if (plot.in.3d){
    
    #if (!w%in%w.vals.sp) rgl.close()
    
  }
  
  if (nrow(far.to.init.1)==0) { far.to.init.1<-NULL}
  if (nrow(far.to.init.2)==0) { far.to.init.2<-NULL}
  if (nrow(close.to.init.1)==0) { close.to.init.1<-NULL}
  if (nrow(close.to.init.2)==0) { close.to.init.2<-NULL}
  
  far.to.init.X5.for.w<- c(far.to.init.X5.for.w, list(far.to.init.1)) 
  far.to.init.X6.for.w<- c(far.to.init.X6.for.w, list(far.to.init.2))
  close.to.init.X5.for.w <- c(close.to.init.X5.for.w, list(close.to.init.1))
  close.to.init.X6.for.w <- c(close.to.init.X6.for.w, list(close.to.init.2)) 
  print("length(far.to.init.X5.for.w)")
  print(length(far.to.init.X5.for.w))
  
  
  
}


min.config.stress.w.table<-rbind(min.config.stress.1.w,min.config.stress.2.w)

colnames(min.config.stress.w.table) <- w.vals
row.names(min.config.stress.w.table) <- c("Local min for real config.","Alternative local min")
#if (plot.in.3d){ rgl.close()}
#graphics.off()



if (create.ani){
  saveGIF(animate.config.w(),"config_w.gif")
  saveGIF(animate.final.loc.stress.w(),"final_loc_vals_w.gif")
  saveMovie(animate.final.loc.stress.w(),movie.name="final_loc_vals_w.mov")
  #saveVideo(animate.final.loc.stress.w(),movie.name="final_loc_vals_w.avi")
  
  open3d(windowRect=c(0,0,480,480))
  play3d(stress.plot3d,duration=65,dev = rgl.cur(),sign.hessian.at.pt=sign.hessian.at.pt,
         x.coords=x.coords,y.coords=y.coords,
         stress.at.loc=stress.at.loc.w,
         grid.seq.x=grid.seq.x,grid.seq.y=grid.seq.y,w.vals=w.vals,rotate.z.angle=pi/3)
  
  movie3d(stress.plot3d,duration=65,dev = rgl.cur(),fps=fps,sign.hessian.at.pt=sign.hessian.at.pt,
          x.coords=x.coords,y.coords=y.coords,
          stress.at.loc=stress.at.loc.w,
          grid.seq.x=grid.seq.x,grid.seq.y=grid.seq.y,w.vals=w.vals, rotate.z.angle=pi/3,
          dir=file.path(getwd(),"./cache/"),type="mov")
  
}