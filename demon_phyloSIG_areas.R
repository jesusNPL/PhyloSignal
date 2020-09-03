
## Function to estimate phylogenetic signal based on randomization and BM for multiple traits

phyloSIG_areas <- function(phylo, dato, area, treatment, nsim) { 
  require(phytools)
  K_val <- numeric(length = length(2:ncol(dato))) 
  P_val <- numeric(length = length(2:ncol(dato)))
  P_val_BM <- numeric(length = length(2:ncol(dato)))
  headers <- character(length = length(2:ncol(dato)))
  simRND <- list()
  simBM <- list()
  
  for(i in 2:length(names(dato))){ 
    trait <- na.omit(setNames(dato[, i], dato[, 1]))
    phylo2 <- drop.tip(phylo, setdiff(phylo$tip.label, names(trait)))
    
    # Simulation randomized tips
    tmp <- phylosig(phylo2, trait, method = "K", nsim = nsim, test = TRUE)
    
    headers[i] <- names(dato)[i]
    K_val[i] <- tmp$K
    P_val[i] <- tmp$P
    simRND[[i]] <- tmp$sim.K
    
    ##### Simulation Brownian model
    simBM[[i]] <- apply(phytools::fastBM(phylo2, nsim = nsim, sig2 = mean(ape::pic(trait, ape::multi2di(phylo2))^2), a = mean(trait)), 
                   2, phytools::phylosig, tree = phylo2)
    P_val_BM[i] <- mean(abs(log(c(tmp$K, simBM[[i]]))) >= abs(log(tmp$K)))
    cat(names(dato)[i])
    print(i)
    
  }
  simRND <- data.frame(do.call(cbind, simRND))
  names(simRND) <- names(dato)[2:ncol(dato)]
  simBM <- data.frame(do.call(cbind, simBM))
  names(simBM) <- names(dato)[2:ncol(dato)]
  
  Area <- rep(area, length(2:ncol(dato)))
  Treatment <- rep(treatment, length(2:ncol(dato)))
  Ks <- data.frame(cbind(headers, K_val, P_val, P_val_BM, Area, Treatment))
  Ks <- Ks[2:nrow(Ks), ] 
  res <- list(Ks, simRND, simBM)
  return(res)
}
