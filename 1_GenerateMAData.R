# Generate Binary Event Meta-Analysis Data

library(future.apply)
plan(multisession, workers = 80)
source("0_DGMs.R")  # Data Generating Methods

# Function to Replicate Data ---------------------------------------------------
SaveDataSets <- function(theta, k, tau2, p_ic_init, min_n, max_n, num_reps,
                         dir) {
  
  pRandom_data_sets <- replicate(num_reps, 
                                 pRandomGenerateMAData(theta = theta,
                                                       k = k,
                                                       tau2 = tau2,
                                                       p_ic_init = p_ic_init,
                                                       min_n = min_n,
                                                       max_n = max_n))
  pCFixed_data_sets <- replicate(num_reps, 
                                 pCfixedGenerateMAData(theta = theta,
                                                       k = k,
                                                       tau2 = tau2,
                                                       p_ic_init = p_ic_init,
                                                       min_n = min_n,
                                                       max_n = max_n))
  pRandomB_data_sets <- replicate(num_reps, 
                                  pRandomBGenerateMAData(theta = theta,
                                                         k = k,
                                                         tau2 = tau2,
                                                         p_ic_init = p_ic_init,
                                                         min_n = min_n,
                                                         max_n = max_n))
  
  # Save as RData files
  save(pRandom_data_sets,
       file = paste0(dir, "/pRandom,Theta=", theta,
                     ",Tau2=", tau2,
                     ",P_ic=", p_ic_init,
                     ",MinN=", min_n,
                     ",MaxN=", max_n, ".RData"))
  save(pCFixed_data_sets,
       file = paste0(dir, "/pCFixed,Theta=", theta,
                     ",Tau2=", tau2,
                     ",P_ic=", p_ic_init,
                     ",MinN=", min_n,
                     ",MaxN=", max_n, ".RData"))
  save(pRandomB_data_sets,
       file = paste0(dir, "/pRandomB,Theta=", theta,
                     ",Tau2=", tau2,
                     ",P_ic=", p_ic_init,
                     ",MinN=", min_n,
                     ",MaxN=", max_n, ".RData"))
}

# Generate DGM combinations with expand.grid and set parameter values ----------
theta_sim_vals <- seq(-1.5, 1.5, 0.5)
tau2_sim_vals <- c(0, 0.2, 0.4, 0.8)
p_ic_init_small <- c(0.01, 0.05, 0.1)
p_ic_init_large <- c(0.01, 0.05)

# small sample sizes
dgm_param_combs_small <- expand.grid("theta" = theta_sim_vals,
                                     "tau2" = tau2_sim_vals, 
                                     "p_ic_init" = p_ic_init_small)
dgm_param_combs_small$min_n <- 10
dgm_param_combs_small$max_n <- 100

# large sample sizes
dgm_param_combs_large <- expand.grid("theta" = theta_sim_vals,
                                     "tau2" = tau2_sim_vals, 
                                     "p_ic_init" = p_ic_init_large)
dgm_param_combs_large$min_n <- 100
dgm_param_combs_large$max_n <- 500

# combine small and large into one data frame
dgm_param_combs <- rbind(dgm_param_combs_large, dgm_param_combs_small)

num_reps <- 2000
k <- 10
seed <- 1234

# Parallelize Code -------------------------------------------------------------
rdsdir <- "./RDSDataSets"
if (!dir.exists(rdsdir))
  dir.create(rdsdir)
##
future_mapply("SaveDataSets",
              theta = dgm_param_combs$theta,
              tau2 = dgm_param_combs$tau2,
              p_ic_init = dgm_param_combs$p_ic_init,
              min_n = dgm_param_combs$min_n,
              max_n = dgm_param_combs$max_n,
              MoreArgs = list(k = k,
                              num_reps = num_reps),
              future.seed = seed,
              dir = rdsdir)
