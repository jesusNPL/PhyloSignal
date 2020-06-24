calcAICw <- function(tabla){
  if ( ! ("geiger" %in% installed.packages())) {install.packages("geiger", dependencies = T)}
  pp <- split(tabla, seq(nrow(tabla)))
  aicValues <- list()
  for(i in 1:length(pp)){
  aicValues[[i]] <- geiger::aicw(as.numeric(pp[[i]]))
  }
  Phylo <-  as.list(paste("Phylo", 1:length(pp), sep = "_"))
  Model <- c("lambda", "BM", "noSignal")
  aicVALS <- do.call(rbind, aicValues)
  aicVALS$Phylo <- rep(Phylo, each = 3)
  aicVALS$Model <- rep(Model, each = 1)
  names(aicVALS) <- c("AICc", "dAICc", "AICcw", "Phylogeny", "Model")
  return(aicVALS)
}

#llll <- calcAICw(tabla = aicss)
