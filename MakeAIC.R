#The general form for calculating AIC:
#  AIC = -2*ln(likelihood) + 2*K  where ln is the natural logarithm, (likelihood) is the 
# value of the likelihood and K is the number of parameters in the model.

# Corrected AIC
#  AICc = -2*ln(likelihood) + 2*K + (2*K*(K+1))/(n-K-1) variables are same as in AIC

### Function to calculate AIC
# Input a vector of likelihood and a vector of the number of parameters of each likelihood
# Output a vector of AIC in the same order
AICcalc <- function(x, y){
  mod <- (x) # likelihood values
  parm <- (y) # number of paramenters
  #aic <- ((-2*log(x)) + (2*y)) # calculate AIC
  aic <- ((-2*(x)) + (2*y)) # calculate AIC with loglikelihood values
  return(aic)
}

AICcClac <- function(x, y, n){
  mod <- (x) # likelihood values
  parm <- (y) # number of paramenters
  n <- 26 # samples in the dataset
  #aicc <- -2*log(x) + 2*y + ((2*y*(y+1))/(n-y-1)) # calculate AICc
  aicc <- -2*(x) + 2*y + ((2*y*(y+1))/(n-y-1)) # calculate AICc
    return(aicc)
}

### Function to calculate Delta AIC
# Input is a vector of AIC values for various models (from the same dataset)
# Output is a vector of deltas in the same order
aic.delta <- function(x){
  l <- min(x) # calculate the lowest AIC score
  #d <- (x-l) # get differences
  d <- x-min(x)
  return(d)
}


### Function to calculate AIC weights (following Posada and Buckley 2004)
# Input is a vector of AIC values for various models (from the same dataset)
# Output is a vector of weights in the same order
make.aic.weight<-function(x){
  l<-min(x) # calculate the lowest AIC score
  d <- exp(-0.5*(x-l)) # get differences and exponentiate
  wg <- d / sum(d)
  return(wg)
}

# Function to calculate AICW using the package qpcR

makeWeights <- function(x){
  W <- list(NA)
  for(i in 1:length(x)){
    print(i)
    dat <- x[[i]]
    w <- akaike.weights(dat[,8])
    W[[i]] <- w
  }
  W <- as.data.frame(do.call("cbind", W))
  return(W)
}

