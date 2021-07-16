

make_K <- function(phylo, trait){
  n <- length(trait)  # trait is X parameter in Blomberg et al. n is parameter n.
  vcv <- ape::vcv.phylo(phylo)  # parameter V
  inverse.vcv <- solve(vcv)  # parameter V^-1
  root.value <- sum(inverse.vcv %*% trait)/sum(inverse.vcv)  # parameter a hat
  # phylogenetic corrected trait mean
  MSEo <- (t(trait - root.value) %*% (trait - root.value))/(n - 1)
  # the mean squared error of the trait
  MSE <- (t(trait - root.value) %*% inverse.vcv %*% (trait - root.value))/(n - 1)
  # the mean squared error of the trait given the VCV
  ratio.obs <- MSEo/MSE  # observed ratio of mean squared errors
  ratio.expected <- (sum(diag(vcv)) - (n/sum(inverse.vcv)))/(n - 1)
  k <- as.numeric(ratio.obs/ratio.expected)  # standardized and can be compared across studies
  
  names(k) <- "Blomberg's_K"
  return(k)
}

#tr <- ape::read.tree("Downloads/comparative.phylo.txt")
#traits <- read.table("Downloads/comparative.traits.txt", 
 #                    sep = "\t", header = T, row.names = 1)

#make_K(phylo = tr, trait = setNames(traits$trait1, rownames(traits)))
