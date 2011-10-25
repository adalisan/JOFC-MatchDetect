
N          <- length(label)
loss <- "strain"
size<-seq(0,1,0.01)




data.subset<-sample(1:N,300)
data.subset<-sort(data.subset)
TE<- TE[data.subset,data.subset]
TF<- TF[data.subset,data.subset]

N<- 300