
library(tidyr)
library(dplyr)
library(phytools)

### Load data
trait_data <- read.csv("Data/Total all gardens per code_TS_CS_complete_mod_VT_over_all_mean.csv")
phylo_DATA <- ape::read.tree("Data/tree_all_sp")

##### Data preparation ####
data = subset(trait_data, select = -c(Garden.x) )

# Raw data
trait_data <- data %>%
  select(sp, VA_mean, VD_mm, VF, S, Dmh_um, Ik_um, 
         Ktheo, Fp_mod, Ap_mod, Dpit, N_vas_trach)

# Log data
trait_data_log <- trait_data
trait_data_log[, 2:12] <- apply(trait_data_log[, 2:12], 1, FUN = log)

### Data aggregation
# raw data mean
trait_mean <- aggregate(. ~ sp, trait_data, FUN = mean) # aggregating to mean value per species

# log data mean
trait_mean_log <- aggregate(. ~ sp, trait_data_log, FUN = mean)

### Phylogenetic data
phylo_sel <- drop.tip(phylo_DATA, setdiff(phylo_DATA$tip.label, trait_mean$sp)) # drop species not in the trait table

### Calculating measurement error
sp_count <- trait_data %>% 
  group_by(sp) %>% 
  count()

nInd <- mean(sp_count$n) # there are 4.6 individuals per species 

trait_SD <- aggregate(. ~ sp, trait_data, FUN = sd)
trait_SD_log <- aggregate(. ~ sp, trait_data_log, FUN = sd)

mean(trait_SD$VA_mean)
mean(trait_SD_log$VA_mean)

seSize <- numeric(length = 11)
seSize_log <- numeric(length = 11)

for(i in 1:length(seSize)) {
  
  tm1 <- trait_SD[, 2:12][, i]
  tm2 <- trait_SD_log[, 2:12][, i]
  
  seSize[i] <- round(mean(tm1)/sqrt(nInd), 4) # Measurement error for raw data
  seSize_log[i] <- round(mean(tm2)/sqrt(nInd), 4) # Measurement error for log data
}

names(seSize) <- names(trait_SD)[2:12]
names(seSize_log) <- names(trait_SD_log)[2:12]

### Standard error per species per trait - We will use these estimates to calculate the PS

SE <- function(x) {
  res <- sd(x)/sqrt(sum(!is.na(x)))
}

trait_SE <- aggregate(. ~ sp, trait_data, FUN = SE)
trait_SE_log <- aggregate(. ~ sp, trait_data_log, FUN = SE)

##### Phylogenetic signal #####
source("Functions/demon_K_test.R")
source("Functions/demon_K_test_ME.R")

### Without measurement error
ps_list <- list()
ps_list_log <- list() 

### With measurement error
ps_list_ME <- list()
ps_list_ME_log <- list() 

for(i in 2:ncol(trait_mean)) { 
  
  print(names(seSize)[i])
  
  # Without measurement error
  tmp <- demon_K_test(tree = phylo_sel, 
                               trait = setNames(trait_mean[, i], trait_mean[, 1]), 
                               method = "K", test = TRUE, nsim = 1000, 
                               bounds_sim = c(-Inf, Inf))
  
  tmp_log <- demon_K_test(tree = phylo_sel, 
                      trait = setNames(trait_mean_log[, i], trait_mean_log[, 1]), 
                      method = "K", test = TRUE, nsim = 1000, 
                      bounds_sim = c(-Inf, Inf))

  tmp$Trait <- names(trait_mean)[i]
  tmp_log$Trait <- names(trait_mean_log)[i]
  
  ps_list[[i]] <- tmp
  ps_list_log[[i]] <- tmp_log
  
  print(paste0("PS estimated for ", names(trait_mean)[i], " without Measurement Error")) 
  
  # With measurement error
  
  tmp_ME <- demon_K_test_ME(tree = phylo_sel, trait = setNames(trait_mean[, i], trait_mean[, 1]), 
                            ME = setNames(trait_SE[, i], trait_SE[, 1]), method = "K", 
                            test = TRUE, nsim = 1000, bounds_sim = c(-Inf, Inf))
  
  tmp_ME_log <- demon_K_test_ME(tree = phylo_sel, trait = setNames(trait_mean_log[, i], trait_mean_log[, 1]), 
                            ME = setNames(trait_SE_log[, i], trait_SE_log[, 1]), method = "K", 
                            test = TRUE, nsim = 1000, bounds_sim = c(-Inf, Inf))
  
  tmp_ME$Trait <- names(trait_mean)[i]
  tmp_ME_log$Trait <- names(trait_mean_log)[i]
  
  ps_list_ME[[i]] <- tmp_ME
  ps_list_ME_log[[i]] <- tmp_ME_log
  
  print(paste0("PS estimated for ", names(trait_mean)[i], " with Measurement Error")) 
  
}

# These ARE the final tables with the phylogenetic signal for all traits

### Without measurement error
# raw
ps_all <- do.call(rbind, ps_list)
ps_all
# log
ps_all_log <- do.call(rbind, ps_list_log)
ps_all_log

write.csv(ps_all, file = "Results_PS/PS_traits_raw.csv")
write.csv(ps_all_log, file = "Results_PS/PS_traits_log.csv")

### With measurement error
# raw
ps_all_ME <- do.call(rbind, ps_list_ME)
ps_all_ME
# log
ps_all_ME_log <- do.call(rbind, ps_list_ME_log)
ps_all_ME_log

write.csv(ps_all, file = "Results_PS/PS_traits_ME_raw.csv")
write.csv(ps_all_log, file = "Results_PS/PS_traits_ME_log.csv")
