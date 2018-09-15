# Function to perform test of phylogenetic signal using Pagel's lambda

lambdaEval <- function(tree, trait){
  library(pbapply)
  lambdaModel <- pblapply(multi2di.multiPhylo(tree, random = TRUE), geiger::fitContinuous, 
                          trait, model = "lambda")
  bmModel <- pblapply(multi2di.multiPhylo(tree, random = TRUE), geiger::fitContinuous, 
                      trait, model = "BM")
  TreeLambda0 <- pblapply(multi2di.multiPhylo(tree, random = TRUE), rescale, 
                          model = "lambda", 0)
  class(TreeLambda0) <- "multiPhylo"
  nosigModel <- pblapply(TreeLambda0, geiger::fitContinuous, trait, model = "BM")
  
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

lambda_model <- function(tree, trait, method = "lambda", test = TRUE, nsim = 999){
  lambdaModel <- pblapply(tree, phytools::phylosig, trait, method, test, nsim)
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

K_model <- function(tree, trait, method = "K", test = TRUE, nsim = 999){
  KModel <- pblapply(tree, phytools::phylosig, trait, method, test, nsim)
  Ks <- list()
  KPval <- list()
  for(j in 1:length(tree)){
    Ks[[j]] <- KModel[[j]]$K
    KPval[[j]] <- KModel[[j]]$P
  }
  kVals <- do.call(rbind, Ks)
  kPvals <- do.call(rbind, KPval)
  
  results <- data.frame("K" = kVals, "KPval" = kPvals)
  return(results)
}

#y <- K_model(tree = mbCCfull100, trait = sla, method = "K", test = TRUE, nsim = 9999)
#Z <- lambda_model(tree = mbCCfull100, trait = sla, method = "lambda", test = TRUE, nsim = 999)

calcAICw <- function(tabla){
  pp <- split(tabla, seq(nrow(tabla)))
  aicValues <- list()
  for(i in 1:100){
  aicValues[[i]] <- geiger::aicw(as.numeric(pp[[i]]))
  }
  Phylo <-  as.list(paste("Phylo", 1:100, sep = "_"))
  Model <- c("lambda", "BM", "noSignal")
  aicVALS <- do.call(rbind, aicValues)
  aicVALS$Phylo <- rep(Phylo, each = 3)
  aicVALS$Model <- rep(Model, each = 1)
  names(aicVALS) <- c("AICc", "Delta", "AICw", "Phylogeny", "Model")
  return(aicVALS)
}

#llll <- calcAICw(tabla = aicss)








