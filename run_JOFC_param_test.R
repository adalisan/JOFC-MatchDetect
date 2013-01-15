args=(commandArgs(TRUE))

paramsFile=NULL
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    paramsFile= "./src/runningParams.R"

}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

if (is.null (paramsFile))
  paramsFile= "./src/runningParams.R"

paramsFile = paste("./src/runParams_par_",Sys.getenv("SGE_TASK_ID"),".R",sep="",collapse="")

library(ProjectTemplate)
library(snowfall)

print(paramsFile)

library(parallel)
debug.mode <- FALSE
run.for.Sweave <- FALSE
num.cpus <- detectCores()
load.project()

source(paramsFile)

print(params)
source("./src/JOFC-test-params.R")

save.image(file= paste("JOFC_MVN_Dir_Sim_",format(Sys.time(), "%b %d %H:%M:%S"),".RData"))

