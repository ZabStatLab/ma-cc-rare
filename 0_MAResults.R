# Run Individual MA Data Set Through metabin/metafor, Separated by Heterogeneity 
# Variance Estimator

source("0_IPM.R")
library(meta)
library(metafor)

MAResults <- function(single_ma_data_set, single_ma_data_set_uncorrected,
                      pooling_method, tau2_estimator) {
  
  # Use metafor for the SA pooling method --------------------------------------
  if (pooling_method == "SA") {
    
    # Fixed-effect results
    ma_results_FE <- metafor::rma.uni(ai = single_ma_data_set$TRT_event,
                                      n1i = single_ma_data_set$TRT_n,
                                      ci = single_ma_data_set$CTRL_event,
                                      n2i = single_ma_data_set$CTRL_n,
                                      measure = "OR",
                                      add = 0,
                                      method = "FE",
                                      weighted = FALSE)  # SA pooling method
    
    # Random-effects results
    if (tau2_estimator == "IPM") {  # Pre-calculate tau2 estimate with IPM code
      
      IPM_est <- IPM(single_ma_data_set_uncorrected = 
                       single_ma_data_set_uncorrected,
                     maxit = 100)
      
      ma_results_RE <- metafor::rma.uni(ai = single_ma_data_set$TRT_event,
                                        n1i = single_ma_data_set$TRT_n,
                                        ci = single_ma_data_set$CTRL_event,
                                        n2i = single_ma_data_set$CTRL_n,
                                        measure = "OR",
                                        add = 0,
                                        tau2 = IPM_est,
                                        weighted = FALSE)  # SA pooling method
      
    } else {  # Use built-in tau2 estimators
      
      ma_results_RE <- metafor::rma.uni(ai = single_ma_data_set$TRT_event,
                                        n1i = single_ma_data_set$TRT_n,
                                        ci = single_ma_data_set$CTRL_event,
                                        n2i = single_ma_data_set$CTRL_n,
                                        measure = "OR",
                                        add = 0,
                                        method = tau2_estimator,
                                        weighted = FALSE)  # SA pooling method
    }
    
    # Combine results
    ma_results <- list("TE.fixed" = as.numeric(ma_results_FE$beta), 
                       "upper.fixed" = ma_results_FE$ci.ub, 
                       "lower.fixed" = ma_results_FE$ci.lb, 
                       "pval.fixed" = ma_results_FE$pval, 
                       "TE.random" = as.numeric(ma_results_RE$beta), 
                       "upper.random" = ma_results_RE$ci.ub, 
                       "lower.random" = ma_results_RE$ci.lb, 
                       "pval.random" = ma_results_RE$pval, 
                       "tau2" = ma_results_RE$tau2)

  } else {  # Use metabin for the other pooling methods ------------------------
    
    if (tau2_estimator == "IPM") {  # Pre-calculate tau2 estimate with IPM code
      
      IPM_est <- IPM(single_ma_data_set_uncorrected = 
                       single_ma_data_set_uncorrected,
                     maxit = 100)
      
      ma_results <- meta::metabin(event.e = single_ma_data_set$TRT_event, 
                                  n.e = single_ma_data_set$TRT_n,
                                  event.c = single_ma_data_set$CTRL_event, 
                                  n.c = single_ma_data_set$CTRL_n,
                                  method = pooling_method,
                                  tau.preset = sqrt(IPM_est),
                                  sm = "OR",
                                  incr = 0)
      
    } else {  # Use built-in tau2 estimators
      
      ma_results <- meta::metabin(event.e = single_ma_data_set$TRT_event, 
                                  n.e = single_ma_data_set$TRT_n,
                                  event.c = single_ma_data_set$CTRL_event, 
                                  n.c = single_ma_data_set$CTRL_n,
                                  method = pooling_method,
                                  method.tau = tau2_estimator,
                                  sm = "OR",
                                  incr = 0)
    }
  }
  
  ma_results_stats <- c(ma_results$TE.fixed,  # fixed LOR
                        ma_results$upper.fixed,  # fixed upper CI
                        ma_results$lower.fixed,  # fixed lower CI
                        ma_results$pval.fixed,  # fixed pval
                        ma_results$TE.random,  # random LOR
                        ma_results$upper.random,  # random upper CI
                        ma_results$lower.random,  # random lower CI
                        ma_results$pval.random,  # random pval
                        ma_results$tau2)  # tau2 estimate
  return(ma_results_stats)
}
