
K_test <- function(tree, trait, method, test, nsim, number_of_trees, bounds_sim){
  if ( ! ("pbapply" %in% installed.packages())) {install.packages("pbapply", dependencies = T)}
  if ( ! ("phytools" %in% installed.packages())) {install.packages("phytools", dependencies = T)}
  
  posteriors <- sample(tree, number_of_trees)
  if (length(posteriors) == 1){tree <- tree} else {tree <- posteriors}
  
  KModel <- pbapply::pblapply(tree, phytools::phylosig, trait, method, test, nsim)
  Trees <- ape::multi2di.multiPhylo(tree, random = TRUE)
  
  Ks <- list()
  KPval <- list()
  KNull <- list()
  
  for(j in 1:length(tree)){
    
    cat("Phylogenetic tree = ", j <- j, "\n")
    tr <- Trees[[j]]
    
    Ks[[j]] <- KModel[[j]]$K
    KPval[[j]] <- KModel[[j]]$P
    
    KNull[[j]] <- apply(phytools::fastBM(tr, n = nsim, sig2 = mean(ape::pic(trait, ape::multi2di(tr))^2), a = mean(trait), bounds = bounds_sim), 
                        2, phytools::phylosig, tree = tr)
    #KPvalNull[[j]] <- mean(abs(log(c(Ks[[j]], KNull[[j]]))) >= abs(log(Ks[[j]])))
  }
  
  kVals <- do.call(rbind, Ks)
  kPvals <- do.call(rbind, KPval)
  kNULL <- do.call(rbind, KNull)
  
  KPvalNull <- list()
  kNullMean <- list()
  kNullSD <- list()
  
  for(i in 1:length(KNull)){
    tmp <- KNull[[i]]
    KPvalNull[[j]] <- mean(abs(log(c(Ks[[i]], tmp))) >= abs(log(Ks[[i]])))
    kNullMean[[i]] <- mean(tmp)
    kNullSD[[i]] <- sd(tmp)
  }
  
  kPvalNulls <- do.call(rbind, KPvalNull)
  kNullMean <- do.call(rbind, kNullMean)
  kNullSD <- do.call(rbind, kNullSD)
  kSES <- ((kVals - kNullMean)/kNullSD)
  
  results <- data.frame("Kobs" = kVals, "Pvalobs" = kPvals, "kNullmean" = kNullMean, 
                        "KNullsd" = kNullSD, "PvalNull" = kPvalNulls, "kSES" = kSES)
  return(results)
}

#K_test(tree = mbCCfull100, trait = sla, method = "K", test = TRUE, nsim = 999, number_of_trees = 10, bounds_sim = c(-Inf, Inf))

