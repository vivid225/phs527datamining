analysisID = function(id){
  source("perceptron_perallel.R")
  id <- as.numeric(id)
  percep_function(id = id)
}

dID = Sys.getenv("PBS_ARRAYID")
analysisID(id = dID)
