# Functions to perform test of phylogenetic signal using Pagel's lambda and Blomberg’s K

# lambdaEval function performs analysis of phylogenetic signal following https://lukejharmon.github.io/ilhabela/instruction/2015/06/02/ContinuousModels/
# Lambda is a tree transformation that stretches tip branches relative to internal branches, 
# making the tree more and more like a complete polytomy. If our estimated lambda = 0, 
# then the traits are inferred to have no phylogenetic signal. Lambda = 1 corresponds to a Brownian motion model; 
# 0 < lambda < 1 is in between.

lambdaEval <- function(tree, trait){
  
  if ( ! ("pbapply" %in% installed.packages())) {install.packages("pbapply", dependencies = T)}
  if ( ! ("geiger" %in% installed.packages())) {install.packages("geiger", dependencies = T)}
  if ( ! ("phytools" %in% installed.packages())) {install.packages("phytools", dependencies = T)}
  
  lambdaModel <- pbapply::pblapply(ape::multi2di.multiPhylo(tree, random = TRUE), geiger::fitContinuous, 
                          trait, model = "lambda")
  bmModel <- pbapply::pblapply(ape::multi2di.multiPhylo(tree, random = TRUE), geiger::fitContinuous, trait, model = "BM")
  TreeLambda0 <- pbapply::pblapply(ape::multi2di.multiPhylo(tree, random = TRUE), geiger::rescale, 
                          model = "lambda", 0.0)
  class(TreeLambda0) <- "multiPhylo"
  nosigModel <- pbapply::pblapply(TreeLambda0, geiger::fitContinuous, trait)
  
  # lambdas
  lambdas <- list()
  lambdalnL <- list()
  lambdaAICc <- list()
  lambdaK <- list()
  # Brownian motion
  bmlnL <- list()
  bmAICc <- list()
  bmK <- list()
  # No signal
  nosiglnL <- list()
  nosigAICc <- list()
  nosigK <- list()
  
  for(i in 1:length(tree)){
    # Lambda
    lambdas[[i]] <- lambdaModel[[i]]$opt$lambda
    lambdalnL[[i]] <- lambdaModel[[i]]$opt$lnL
    lambdaAICc[[i]] <- lambdaModel[[i]]$opt$aicc
    lambdaK[[i]] <- lambdaModel[[i]]$opt$k
    # Brownian motion
    bmlnL[[i]] <- bmModel[[i]]$opt$lnL
    bmAICc[[i]] <- bmModel[[i]]$opt$aicc
    bmK[[i]] <- bmModel[[i]]$opt$k
    # no signal
    nosiglnL[[i]] <- nosigModel[[i]]$opt$lnL
    nosigAICc[[i]] <- nosigModel[[i]]$opt$aicc
    nosigK[[i]] <- nosigModel[[i]]$opt$k
  }
  lambdas <- do.call(rbind, lambdas)
  lambdalnL <- do.call(rbind, lambdalnL)
  lambdaAICc <- do.call(rbind, lambdaAICc)
  lambdaK <- do.call(rbind, lambdaK)
  
  bmlnL <- do.call(rbind, bmlnL)
  bmAICc <- do.call(rbind, bmAICc)
  bmK <- do.call(rbind, bmK)
  
  nosiglnL <- do.call(rbind, nosiglnL)
  nosigAICc <- do.call(rbind, nosigAICc)
  nosigK <- do.call(rbind, nosigK)
  
  results <- data.frame(lambdas, lambdalnL, lambdaAICc, lambdaK, 
                        bmlnL, bmAICc, bmK, 
                        nosiglnL, nosigAICc, nosigK)
  names(results) <- c("lambda", "lblnL", "lbAICc", "lbK",
                      "bmlnL", "bmAICc", "bmK", 
                      "nosiglnL", "nosigAICc", "nosigK")
  return(results)
}  

#x <- lambdaEval(tree = mbCCfull100, trait = sla)

# Lambda is a tree transformation that stretches tip branches relative to internal branches, 
# making the tree more and more like a complete polytomy. If our estimated lambda = 0, 
# then the traits are inferred to have no phylogenetic signal. Lambda = 1 corresponds to a Brownian motion model; 
# 0 < lambda < 1 is in between.

lambda_model <- function(tree, trait, method = "lambda", test = TRUE, nsim = 999){
    
  if ( ! ("pbapply" %in% installed.packages())) {install.packages("pbapply", dependencies = T)}
  if ( ! ("phytools" %in% installed.packages())) {install.packages("phytools", dependencies = T)}
  
  lambdaModel <- pbapply::pblapply(tree, phytools::phylosig, trait, method, test, nsim)
  lambdas <- list()
  lambdaslogL <- list()
  lambdaslogl0 <- list()
  lambdasPval <- list()
  for(j in 1:length(tree)){
    lambdas[[j]] <- lambdaModel[[j]]$lambda
    lambdaslogL[[j]] <- lambdaModel[[j]]$logL
    lambdaslogl0[[j]] <- lambdaModel[[j]]$logL0
    lambdasPval[[j]] <- lambdaModel[[j]]$P
  }
  lambdaVals <- do.call(rbind, lambdas)
  lambdaLogL <- do.call(rbind, lambdaslogL)
  lambdalogl0 <- do.call(rbind, lambdaslogl0)
  lambdaPvals <- do.call(rbind, lambdasPval)
  
  results <- data.frame("lambda" = lambdaVals, "logL" = lambdaLogL, 
                        "logL0" = lambdalogl0, "Pval" = lambdaPvals)
  return(results)
}

#Z <- lambda_model(tree = mbCCfull100, trait = sla, method = "lambda", test = TRUE, nsim = 999)

# Blomberg’s K, compares the variance of PICs to what we would espect under a Brownian motion model. 
# K = 1 means that relatives resemble one another as much as we should expect under BM; 
# K < 1 means that there is less “phylogenetic signal” than expected under BM, while K > 1 means that there is more. 
# A significant p-value returned from phylosignal tells you that there is significant phylogenetic signal - that is, 
# close relatives are more similar than random pairs of species.

K_model <- function(tree, trait, method = "K", test = TRUE, nsim = 999){
  
  if ( ! ("pbapply" %in% installed.packages())) {install.packages("pbapply", dependencies = T)}
  if ( ! ("phytools" %in% installed.packages())) {install.packages("phytools", dependencies = T)}
  
  KModel <- pbapply::pblapply(tree, phytools::phylosig, trait, method, test, nsim)
  Ks <- list()
  nSP <- list()
  KPval <- list()
  for(j in 1:length(tree)){
    Ks[[j]] <- KModel[[j]]$K
    nSP[[j]] <- length(tree[[j]]$tip.label)
    KPval[[j]] <- KModel[[j]]$P
  }
  kVals <- do.call(rbind, Ks)
  ageRoot <- tree[[1]]$root.time
  ageStart <- floor(ageRoot)
  ageEnd <- 0
  timeVals <- rev(seq(ageEnd, ageStart, by = 1))
  nSPVals <- do.call(rbind, nSP)
  kPvals <- do.call(rbind, KPval)
  
  results <- data.frame("nSPP" = nSPVals, "K" = kVals, "KPval" = kPvals, "Time" = timeVals)
  return(results)
}

#y <- K_model(tree = mbCCfull100, trait = sla, method = "K", test = TRUE, nsim = 9999)

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
