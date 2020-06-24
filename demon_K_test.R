demon_K_test <- function(tree, trait, method, test, nsim, bounds_sim){
  if ( ! ("phytools" %in% installed.packages())) {install.packages("phytools", dependencies = TRUE)}

  Ks <- phytools::phylosig(tree, trait, method, test, nsim)
  K_obs <- Ks$K
  K_pval <- Ks$P
  K_sim <- Ks$sim.K
  K_sim_mean <- mean(K_sim)
  K_sim_SD <- sd(K_sim)
  K_sim_SES <- ((K_obs - mean(K_sim))/sd(K_sim))
  
  K_null <- apply(phytools::fastBM(tree, n = nsim, 
                                   sig2 = mean(ape::pic(trait, ape::multi2di(tree))^2), 
                                   a = mean(trait), bounds = bounds_sim), 2, 
                  phytools::phylosig, tree = tree)
  
  K_null_pval <- mean(abs(log(c(K_obs, K_null))) >= abs(log(K_obs)))
  K_null_mean <- mean(K_null)
  K_null_SD <- sd(K_null)
  K_null_SES <- ((K_obs - mean(K_null))/sd(K_null))
  
  tabla <- data.frame("K_obs" = K_obs, "K_pval" = K_pval, "K_sim_mean" = K_sim_mean, "K_sim_SD" = K_sim_SD, "K_sim_SES" = K_sim_SES, 
                      "K_null_pval" = K_null_pval, "K_null_mean" = K_null_mean, "K_null_SD" = K_null_SD, "K_null_SES" = K_null_SES)
  return(tabla)
}
