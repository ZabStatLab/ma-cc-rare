# Apply Various Continuity Corrections, Heterogeneity Estimators, and Pooling
# Methods to Previously Generated Data Sets

require(devtools)
install_version("meta", 
                version = "4.16",
                repos = "http://cran.us.r-project.org")

library(future.apply)
plan(multisession, workers = 2)
source("0_ApplyAllCCs.R")
source("0_MAResults.R")

AllpCFixedMAResults <- function(theta, tau2, p_ic_init, min_n, max_n,
                                pooling_method, tau2_estimator, excludeDBZ,
                                num_reps, scenario_num,
                                datadir, resdir) {
  
  # Create list to store results -----------------------------------------------
  #   Each item in the list contains the results for the ith data set.
  #   Each item has 12 MA output statistics (columns) for the 14 data set 
  #   versions (rows)
  single_ma_data_set_results <- matrix(ncol = 12, nrow = 14)
  row.names(single_ma_data_set_results) <- c("original_data_set",
                                             "cac_cc_1", "cac_cc_0.01", 
                                             "czc_cc_1", "czc_cc_0.01", 
                                             "cas_cc_1", "cas_cc_0.01",
                                             "ta_cc_1", "ta_cc_0.01", 
                                             "e_cc_1_fixed", "e_cc_1_random",
                                             "e_cc_0.01_fixed", 
                                             "e_cc_0.01_random",
                                             "exclude_zeros_cc")
  colnames(single_ma_data_set_results) <- c("TE.fixed", "upper.fixed",
                                            "lower.fixed", "pval.fixed", 
                                            "TE.random", "upper.random", 
                                            "lower.random", "pval.random",
                                            "tau2", "num.single.zeros",
                                            "num.double.zeros", 
                                            "e.cc.approx.used")
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
    
    # Apply CCs ----------------------------------------------------------------
    # cc_results will be a list of 14 versions of single_ma_data_set followed 
    # by 4 items noting the number of times the approximate solution was used 
    # for the Empirical CC method
    cc_results <- ApplyAllCCs(single_ma_data_set = single_ma_data_set, 
                              pooling_method = pooling_method, 
                              tau2_estimator = tau2_estimator,
                              excludeDBZ = excludeDBZ)
    
    # Note the number of studies where the approximate solution was used
    ma_results[[j]]["e_cc_1_fixed", "e.cc.approx.used"] <- 
      cc_results$e_cc_1_fixed_approx
    ma_results[[j]]["e_cc_1_random", "e.cc.approx.used"] <- 
      cc_results$e_cc_1_random_approx
    ma_results[[j]]["e_cc_0.01_fixed", "e.cc.approx.used"] <- 
      cc_results$e_cc_0.01_fixed_approx
    ma_results[[j]]["e_cc_0.01_random", "e.cc.approx.used"] <- 
      cc_results$e_cc_0.01_random_approx
    
    # Run data set(s) through methods ------------------------------------------
    # For MH, original data set can be analyzed without a CC applied
    if (pooling_method == "MH") {
      
      for (i in 1:nrow(ma_results[[j]])) {
        single_ma_data_set_cc <- cc_results[[i]] 
        ma_results[[j]][i, 1:9] <- MAResults(single_ma_data_set = 
                                               single_ma_data_set_cc, 
                                             single_ma_data_set_uncorrected = 
                                               single_ma_data_set, 
                                             pooling_method = pooling_method, 
                                             tau2_estimator = tau2_estimator)
      }
      
    } else {  # For other pooling methods, skip analyzing original data
      
      for (i in 1:(nrow(ma_results[[j]]) - 1)) {  # -1 to skip original data
        # +1 to skip original data
        single_ma_data_set_cc <- cc_results[[i + 1]]  
        ma_results[[j]][i + 1, 1:9] <- MAResults(single_ma_data_set = 
                                                   single_ma_data_set_cc, 
                                                 single_ma_data_set_uncorrected = 
                                                   single_ma_data_set, 
                                                 pooling_method = 
                                                   pooling_method, 
                                                 tau2_estimator = 
                                                   tau2_estimator)
      }
    }
  }
  
  # save results ---------------------------------------------------------------
  saveRDS(ma_results, paste0(resdir, "/pCFixedResults_", scenario_num, ".rds"))
}

# Generate MA combinations with expand.grid ------------------------------------
theta_sim_vals <- seq(-1.5, 1.5, 0.5)
tau2_sim_vals <- c(0, 0.2, 0.4, 0.8)
p_ic_init_small <- c(0.01, 0.05, 0.1)
p_ic_init_large <- c(0.01, 0.05)
pooling_method <- c("MH", "Inverse", "SSW", "SA")
tau2_estimator <- c("DL", "SJ", "PM", "IPM") 
excludeDBZ <- c(TRUE, FALSE)

# small sample sizes
param_combs_small <- expand.grid("theta" = theta_sim_vals,
                                 "tau2" = tau2_sim_vals, 
                                 "p_ic_init" = p_ic_init_small,
                                 "pooling_method" = pooling_method,
                                 "tau2_estimator" = tau2_estimator,
                                 "excludeDBZ" = excludeDBZ)
param_combs_small$min_n <- 10
param_combs_small$max_n <- 100

# large sample sizes
param_combs_large <- expand.grid("theta" = theta_sim_vals,
                                 "tau2" = tau2_sim_vals, 
                                 "p_ic_init" = p_ic_init_large,
                                 "pooling_method" = pooling_method,
                                 "tau2_estimator" = tau2_estimator,
                                 "excludeDBZ" = excludeDBZ)
param_combs_large$min_n <- 100
param_combs_large$max_n <- 500

# combine small and large into one data frame
param_combs <- rbind(param_combs_large, param_combs_small)

num_reps <- 2
seed <- 1234

# Parallelize Code -------------------------------------------------------------
rdsdir <- "./RDSDataSets"
resdir <- "./Results"
if (!dir.exists(resdir))
  dir.create(resdir)

future_mapply("AllpCFixedMAResults",
              scenario_num = row.names(param_combs),
              theta = param_combs$theta,
              tau2 = param_combs$tau2,
              min_n = param_combs$min_n,
              max_n = param_combs$max_n,
              p_ic_init = param_combs$p_ic_init,
              pooling_method = param_combs$pooling_method,
              tau2_estimator = param_combs$tau2_estimator,
              excludeDBZ = param_combs$excludeDBZ,
              MoreArgs = list(num_reps = num_reps),
              future.seed = seed,
              datadir = rdsdir,
              resdir = resdir)
