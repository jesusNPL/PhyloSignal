lambdaEval <- function(tree, trait){
  options(scipen = 999)
  if ( ! ("geiger" %in% installed.packages())) {install.packages("geiger", dependencies = T)}
  if ( ! ("phytools" %in% installed.packages())) {install.packages("phytools", dependencies = T)}
  
  tree0 <- geiger::rescale(tree, model = "lambda", 0.0)
  
  lambdaModel <- geiger::fitContinuous(tree, trait, model = "lambda")
  bmModel <- geiger::fitContinuous(tree, trait, model = "BM")
  nosigModel <- geiger::fitContinuous(tree0, trait)
  
  # Lambda
  lambda_VAL <- lambdaModel$opt$lambda
  lambda_lnL <- lambdaModel$opt$lnL
  lambda_AICc <- lambdaModel$opt$aicc
  lambda_Npar <- lambdaModel$opt$k
  
  lambda <- data.frame(lambda_VAL, lambda_lnL, lambda_AICc, lambda_Npar)
  colnames(lambda) <- c("Value", "lnL", "AICc", "Npar")
  
  # Brownian motion
  BM_VAL <- NA
  BM_lnL <- bmModel$opt$lnL
  BM_AICc <- bmModel$opt$aicc
  BM_Npar <- bmModel$opt$k
  
  BM <- data.frame(BM_VAL, BM_lnL, BM_AICc, BM_Npar)
  colnames(BM) <- c("Value", "lnL", "AICc", "Npar")
  
  # no signal
  nosig_VAL <- NA
  nosig_lnL <- nosigModel$opt$lnL
  nosig_AICc <- nosigModel$opt$aicc
  nosig_Npar <- nosigModel$opt$k
  
  nosig <- data.frame(nosig_VAL, nosig_lnL, nosig_AICc, nosig_Npar)
  colnames(nosig) <- c("Value", "lnL", "AICc", "Npar")
  
  models <- rbind.data.frame(lambda, BM, nosig)
  rownames(models) <- c("Lambda", "BM", "noSignal")
 
   ### Model comparison
  comparison <- as.numeric(models[, 3])
  names(comparison) <- c("Lambda", "BM", "noSignal")
  modelComparison <- aicw(comparison)
  
  results <- list(models, modelComparison)
  return(results)
}  



