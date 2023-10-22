# Apply the GLMM method to Previously Generated Data Sets

library(future.apply)
plan(multisession, workers = 2)
library(metafor)
library(BiasedUrn)

AllpCFixedMAResults <- function(theta, tau2, p_ic_init, min_n, max_n,
                                excludeDBZ, num_reps, scenario_num,
                                datadir, resdir) {
  
  # Create list to store results -----------------------------------------------
  #   Each item in the list contains the results for the ith data set.
  #   Each item has 12 MA output statistics (columns) for the 14 data set 
  #   versions (rows)
  single_ma_data_set_results <- matrix(ncol = 12, nrow = 1)
  row.names(single_ma_data_set_results) <- c("GLMM")
  colnames(single_ma_data_set_results) <- c("TE.fixed", "upper.fixed",
                                            "lower.fixed", "pval.fixed", 
                                            "TE.random", "upper.random", 
                                            "lower.random", "pval.random",
                                            "tau2", "num.single.zeros",
                                            "num.double.zeros",
                                            "glmm.did.not.work")
  ma_results <- rep(list(single_ma_data_set_results), num_reps)
  
  # Read in data ---------------------------------------------------------------
  load(paste0(
    datadir, "/pCFixed,Theta=", theta, ",Tau2=", tau2, ",P_ic=", p_ic_init, 
              ",MinN=", min_n, ",MaxN=", max_n, ".RData"))
  
  
  # Run data set through CC and MA methods -------------------------------------
  for (j in 1:num_reps) {
    
    single_ma_data_set <- as.data.frame(pCFixed_data_sets[, j])
    
    # Count the number of single- and double-zero studies in the original data -
    TRT_no_event <- single_ma_data_set$TRT_n - single_ma_data_set$TRT_event
    CTRL_no_event <- single_ma_data_set$CTRL_n - single_ma_data_set$CTRL_event
    
    # Double-zeros ---------------------------
    # double-zero events
    double_zeros_loc <- (single_ma_data_set$TRT_event + 
                           single_ma_data_set$CTRL_event == 0) |
      # double-zero non-events
      (TRT_no_event + CTRL_no_event == 0) |
      # double-zero with no TRT event and no CTRL non-event
      (single_ma_data_set$TRT_event == 0 & CTRL_no_event == 0) |
      # double-zero with no CTRL event and no TRT non-event
      (single_ma_data_set$CTRL_event == 0 & TRT_no_event == 0)
    ma_results[[j]][, "num.double.zeros"] <- sum(double_zeros_loc)
    
    # Single-zeros ---------------------------
    # zero events
    single_zeros_loc <- ((single_ma_data_set$TRT_event == 0 | 
                            single_ma_data_set$CTRL_event == 0) |
                           # zero non-events
                           (TRT_no_event == 0 | 
                              CTRL_no_event == 0)) &
      # not a double-zero
      !double_zeros_loc
    ma_results[[j]][, "num.single.zeros"] <- sum(single_zeros_loc)
    
    # If TRUE, exclude double-zero studies
    if (excludeDBZ) { 
      single_ma_data_set <- single_ma_data_set[!double_zeros_loc, ]
    }
    
    # Run data set(s) through GLMM method --------------------------------------
    if (tau2 == 0) {
      
      ma_glmm_res <- rma.glmm(measure = "OR", 
                              ai = single_ma_data_set$TRT_event, 
                              n1i =single_ma_data_set$TRT_n, 
                              ci = single_ma_data_set$CTRL_event, 
                              n2i = single_ma_data_set$CTRL_n, 
                              model = "CM.AL",
                              method = "EE")
      ma_results[[j]][, c(1:9, 12)] <- c("TE.fixed" = ma_glmm_res$beta, 
                                         "upper.fixed" = ma_glmm_res$ci.ub, 
                                         "lower.fixed" = ma_glmm_res$ci.lb, 
                                         "pval.fixed" = ma_glmm_res$pval,
                                         "TE.random" = NA, 
                                         "upper.random" = NA, 
                                         "lower.random" = NA, 
                                         "pval.random" = NA, 
                                         "tau2" = NA,
                                         "glmm.did.not.work" = 0)
    } else {
      
      ma_results[[j]][, c(1:9, 12)] <- tryCatch({
        ma_glmm_res <- rma.glmm(measure = "OR", 
                                ai = single_ma_data_set$TRT_event, 
                                n1i =single_ma_data_set$TRT_n, 
                                ci = single_ma_data_set$CTRL_event, 
                                n2i = single_ma_data_set$CTRL_n, 
                                model = "CM.AL",
                                control = list(optCtrl = list(iter.max = 1000)))
        c("TE.fixed" = NA, 
          "upper.fixed" = NA, 
          "lower.fixed" = NA, 
          "pval.fixed" = NA,
          "TE.random" = ma_glmm_res$beta, 
          "upper.random" = ma_glmm_res$ci.ub, 
          "lower.random" = ma_glmm_res$ci.lb, 
          "pval.random" = ma_glmm_res$pval, 
          "tau2" = ma_glmm_res$tau2,
          "glmm.did.not.work" = 0)
        
      }, 
      error = function(msg) {
        ma_glmm_res <- rma.glmm(measure = "OR", 
                                ai = single_ma_data_set$TRT_event, 
                                n1i =single_ma_data_set$TRT_n, 
                                ci = single_ma_data_set$CTRL_event, 
                                n2i = single_ma_data_set$CTRL_n, 
                                model = "CM.AL",
                                method = "EE")
        c("TE.fixed" = NA, 
          "upper.fixed" = NA, 
          "lower.fixed" = NA, 
          "pval.fixed" = NA,
          "TE.random" = ma_glmm_res$beta, 
          "upper.random" = ma_glmm_res$ci.ub, 
          "lower.random" = ma_glmm_res$ci.lb, 
          "pval.random" = ma_glmm_res$pval, 
          "tau2" = ma_glmm_res$tau2,
          "glmm.did.not.work" = 1)
      })
    }
  }
  
  # save results ---------------------------------------------------------------
  saveRDS(ma_results, paste0(resdir, "/pCFixedResults_", scenario_num, "_GLMM.rds"))
}

# Generate MA combinations with expand.grid ------------------------------------
theta_sim_vals <- seq(-1.5, 1.5, 0.5)
tau2_sim_vals <- c(0, 0.2, 0.4, 0.8)
p_ic_init_small <- c(0.01, 0.05, 0.1)
p_ic_init_large <- c(0.01, 0.05)
excludeDBZ <- c(TRUE, FALSE)

# small sample sizes
param_combs_small <- expand.grid("theta" = theta_sim_vals,
                                 "tau2" = tau2_sim_vals, 
                                 "p_ic_init" = p_ic_init_small,
                                 "excludeDBZ" = excludeDBZ)
param_combs_small$min_n <- 10
param_combs_small$max_n <- 100

# large sample sizes
param_combs_large <- expand.grid("theta" = theta_sim_vals,
                                 "tau2" = tau2_sim_vals, 
                                 "p_ic_init" = p_ic_init_large,
                                 "excludeDBZ" = excludeDBZ)
param_combs_large$min_n <- 100
param_combs_large$max_n <- 500

# combine small and large into one data frame
param_combs <- rbind(param_combs_large, param_combs_small)

num_reps <- 2
seed <- 1234

# Parallelize Code -------------------------------------------------------------
rdsdir <- "./RDSDataSets"
resdir <- "./Results_GLMM"
if (!dir.exists(resdir))
  dir.create(resdir)

future_mapply("AllpCFixedMAResults",
              scenario_num = row.names(param_combs),
              theta = param_combs$theta,
              tau2 = param_combs$tau2,
              min_n = param_combs$min_n,
              max_n = param_combs$max_n,
              p_ic_init = param_combs$p_ic_init,
              excludeDBZ = param_combs$excludeDBZ,
              MoreArgs = list(num_reps = num_reps),
              future.seed = seed,
              datadir = rdsdir,
              resdir = resdir)
