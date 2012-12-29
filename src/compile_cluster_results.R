abs.dir <-paste(c("./cluster/"),collapse="")
abs.dir<- c("./")	
params.list<-list()
sim.res.list<-list()
	print(abs.dir)
	mzdatafiles <- list.files(abs.dir, pattern = "*.RData+", recursive = TRUE, full.names = TRUE)
	for (filename in mzdatafiles){
		print(filename)
      load(file=filename)
       params.list<- c(params.list,list(params))
       sim.res.list<-c(sim.res.list,list(sim.res))
 	}
