demon_K_test_ME <- function(tree, trait, method, test, nsim, ME, bounds_sim){
  if ( ! ("phytools" %in% installed.packages())) {install.packages("phytools", dependencies = TRUE)}

  Ks <- phytools::phylosig(tree = tree, x = trait, method = "K", test = TRUE, se = ME, nsim)
  K_obs <- Ks$K
  K_pval <- Ks$P
  K_sim <- Ks$sim.K
  K_sim_mean <- mean(K_sim)
  K_sim_SD <- sd(K_sim)
  K_sim_SES <- ((K_obs - mean(K_sim))/sd(K_sim))
  
  K_null <- apply(phytools::fastBM(tree, n = nsim, 
                                   sig2 = mean(ape::pic(trait, ape::multi2di(tree))^2), 
                                   a = mean(trait), bounds = bounds_sim), 2, 
                  phytools::phylosig, tree = tree, se = ME)
  
  K_null_vec <- numeric(length = nsim)
  for(i in 1:length(K_null)) {
    K_null_vec[i] <- K_null[[i]]$K
  }
  K_null_pval <- mean(abs(log(c(K_obs, K_null_vec))) >= abs(log(K_obs)))
  K_null_mean <- mean(K_null_vec)
  K_null_SD <- sd(K_null_vec)
  K_null_SES <- ((K_obs - mean(K_null_vec))/sd(K_null_vec))
  
  tabla <- data.frame("K_obs" = K_obs, "K_pval" = K_pval, "K_sim_mean" = K_sim_mean, "K_sim_SD" = K_sim_SD, "K_sim_SES" = K_sim_SES, 
                      "K_BM_pval" = K_null_pval, "K_BM_mean" = K_null_mean, "K_BM_SD" = K_null_SD, "K_BM_SES" = K_null_SES)
  return(tabla)
}

#demon_K_test_ME(tree = phylo_sel, 
 #               trait = setNames(trait_mean_log[, 2], trait_mean_log[, 1]),  
  #              ME = setNames(trait_SE_log[, 2], trait_SE_log[, 1]), 
   #             method = "K", test = TRUE, nsim = 1000, 
    #            bounds_sim = c(-Inf, Inf))


