# Code to Implement Various Continuity Corrections to Individual MA Data Sets

# 1. Constant to All Cells CC (CAC_CC)
# 2. Constant to Zero Cells CC (CZC_CC)
# 3. Constant to All Studies CC (CAS_CC)
# 4. Treatment Arm CC (TA_CC)
# 5. Empirical CC (E_CC)
# 6. Exclude all single- and double-zero studies (Exclude_Zeros, Exlucde_SZ_DZ)

source("0_MAResults.R")  # needed for Empirical CC
library(nleqslv)  # needed for Empirical CC

################################################################################
# 1. Constant to All Cells CC (CAC-CC)
################################################################################
# Adds a small value to all cells in both treatment arms in each of the
# studies that contain any zero observed events

CAC_CC <- function(single_ma_data_set, cc) {
  
  TRT_no_event <- single_ma_data_set$TRT_n - single_ma_data_set$TRT_event
  CTRL_no_event <- single_ma_data_set$CTRL_n - single_ma_data_set$CTRL_event
  
  for (i in 1:length(single_ma_data_set$TRT_event)) {
    # if there are zero observed events (or non-events) in either treatment arm  
    if (sum(single_ma_data_set[i, ] == 0) > 0 | 
        TRT_no_event[i] == 0 | 
        CTRL_no_event[i] == 0) {
      row_i <- single_ma_data_set[i, ]
      single_ma_data_set[i, "TRT_event"] <- row_i$TRT_event + cc
      single_ma_data_set[i, "CTRL_event"] <- row_i$CTRL_event + cc
      single_ma_data_set[i, "TRT_n"] <- row_i$TRT_n + 2 * cc
      single_ma_data_set[i, "CTRL_n"] <- row_i$CTRL_n + 2 * cc
    }
  }
  return(single_ma_data_set)
}

################################################################################
# 2. Constant to Zero Cells CC (CZC_CC)
################################################################################
# Adds a small value only to the cells containing zero observed events in each 
# study 

CZC_CC <- function(single_ma_data_set, cc){
  
  # correct for zero observed events
  single_ma_data_set$TRT_event <- ifelse(single_ma_data_set$TRT_event == 0, 
                                         cc,
                                         single_ma_data_set$TRT_event)
  single_ma_data_set$CTRL_event <- ifelse(single_ma_data_set$CTRL_event == 0,
                                          cc, 
                                          single_ma_data_set$CTRL_event)
  
  # correct for zero observed non-events
  TRT_no_event <- single_ma_data_set$TRT_n - single_ma_data_set$TRT_event
  CTRL_no_event <- single_ma_data_set$CTRL_n - single_ma_data_set$CTRL_event
  single_ma_data_set$TRT_n <- ifelse(TRT_no_event == 0, 
                                     single_ma_data_set$TRT_n + cc,
                                     single_ma_data_set$TRT_n)
  single_ma_data_set$CTRL_n <- ifelse(CTRL_no_event == 0,
                                      single_ma_data_set$CTRL_n + cc, 
                                      single_ma_data_set$CTRL_n)
  return(single_ma_data_set)
}

################################################################################
# 3. Constant to All Studies CC (CAS_CC)
################################################################################
# Adds a small value to all cells in all studies in the meta-analysis, 
# regardless if any study contains zero observed events 

CAS_CC <- function(single_ma_data_set, cc) {
  
  single_ma_data_set$TRT_event <- single_ma_data_set$TRT_event + cc
  single_ma_data_set$CTRL_event <- single_ma_data_set$CTRL_event + cc
  single_ma_data_set$TRT_n <- single_ma_data_set$TRT_n + (2 * cc)
  single_ma_data_set$CTRL_n <- single_ma_data_set$CTRL_n + (2 * cc)
  
  return(single_ma_data_set)
}

################################################################################
# 4. Treatment Arm CC (TA_CC)
################################################################################
# Developed by Sweeting et al, 2004
# Adds small, study- and treatment-specific values to all cells in studies 
# containing zero observed events. 
# Note: n_t * R = n_C, so no need to compute R (the sample size imbalance) 

TA_CC <- function(single_ma_data_set, prop_constraint) {
  
  # Calculate a vector of study-specific corrections in the TRT arm (TCC, k_T)
  # and CTRL arm (CCC, k_C)
  TCC <- single_ma_data_set$TRT_n / 
    (single_ma_data_set$TRT_n + single_ma_data_set$CTRL_n) * prop_constraint
  CCC <- single_ma_data_set$CTRL_n / 
    (single_ma_data_set$TRT_n + single_ma_data_set$CTRL_n) * prop_constraint
  
  # If a study has no zero events (or zero non-events), then set the correction 
  # to zero
  TRT_no_event <- single_ma_data_set$TRT_n - single_ma_data_set$TRT_event
  CTRL_no_event <- single_ma_data_set$CTRL_n - single_ma_data_set$CTRL_event
  study_w_no_events <- single_ma_data_set$TRT_event != 0 & 
    single_ma_data_set$CTRL_event != 0 &
    TRT_no_event != 0 & 
    CTRL_no_event != 0
  TCC[study_w_no_events] <- 0
  CCC[study_w_no_events] <- 0
  
  # If ccs are less than zero, throw error
  if (sum(TCC < 0) > 0  | sum(CCC < 0) > 0) {
    stop("Treatment Arm CC Corrections are Negative")
  } 
  
  # Apply the corrections
  single_ma_data_set$TRT_event <- single_ma_data_set$TRT_event + TCC
  single_ma_data_set$TRT_n <- single_ma_data_set$TRT_n + (2 * TCC)
  single_ma_data_set$CTRL_event <- single_ma_data_set$CTRL_event + CCC
  single_ma_data_set$CTRL_n <- single_ma_data_set$CTRL_n + (2 * CCC)
  
  return(single_ma_data_set)
}

################################################################################
# 5. Empirical CC (E_CC)
################################################################################
# Developed by Sweeting et al, 2004
# Adds small, study- and treatment-specific values to all cells in studies 
# containing zero observed events. Prior odds ratios are computed that influence
# the correction value.

# Called by E_CC: solve nonlinear system of two equations ----------------------
GetEmpiricalCCs <- function(n_it, n_ic, omega_hat, prop_constraint) {
  
  # System of two equations:
  SysTwoEq <- function(x) {
    c(x[1] * n_ic - omega_hat * x[2] * n_it - 
        x[1] * x[2] * (omega_hat  - 1) - 0,
      x[1] + x[2] - prop_constraint)
  }
  
  # Jacobian matrix of two equations:
  SysTwoEqJac <- function(x) {
    n <- length(x)
    Df <- matrix(numeric(n * n), n, n)
    # partial derivative wrt x1 of equation 1:
    Df[1, 1] <- n_ic - x[2] * (omega_hat  - 1)  
    # partial derivative wrt x2 of equation 1:
    Df[1, 2] <- omega_hat * n_it - x[1] * (omega_hat  - 1)  
    # partial derivative wrt x1 of equation 2:
    Df[2, 1] <- 1  
    # partial derivative wrt x2 of equation 2:
    Df[2, 2] <- 1  
    return(Df)
  }
  
  sol <- nleqslv(c(1, 1), SysTwoEq, jac = SysTwoEqJac)
  
  # if exact solution converged and corrections are both positive, use exact
  if (sol$termcd == 1 & sol$x[1] >= 0 & sol$x[2] >= 0) {  
    used_approx <- FALSE
    TCC <- sol$x[1]
    CCC <- sol$x[2]
  } else {  # if exact solution fails, use approximate solution
    used_approx <- TRUE
    TCC <- (n_it * omega_hat) / (n_it * omega_hat + n_ic) * prop_constraint
    CCC <- n_ic / (n_it * omega_hat + n_ic) * prop_constraint
  }
  return(list(c(TCC, CCC), used_approx))
}


E_CC <- function(single_ma_data_set, prop_constraint, tau2_estimator,
                 pooling_method) {
  
  # Get prior log odds ratios --------------------------------------------------
  # remove studies with zero observed events and zero observed non-events 
  TRT_no_event <- single_ma_data_set$TRT_n - single_ma_data_set$TRT_event
  CTRL_no_event <- single_ma_data_set$CTRL_n - single_ma_data_set$CTRL_event
  data_no_zeros <- single_ma_data_set[!(single_ma_data_set$TRT_event == 0 | 
                                          single_ma_data_set$CTRL_event == 0 |
                                          TRT_no_event == 0 |
                                          CTRL_no_event == 0), ]
  
  ma_results <- MAResults(single_ma_data_set = data_no_zeros, 
                          single_ma_data_set_uncorrected = data_no_zeros, 
                          pooling_method = pooling_method, 
                          tau2_estimator = tau2_estimator)
  omega_hat_fixed <- exp(ma_results[[1]])  # prior fixed-effect OR
  omega_hat_random <- exp(ma_results[[5]])  # prior random-effects OR
  
  # solve nonlinear system of two equations to get corrections -----------------
  ecc_fixed <- rep(list(c(0, 0)), length(single_ma_data_set$TRT_event))
  ecc_random <- rep(list(c(0, 0)), length(single_ma_data_set$TRT_event))
  ecc_fixed_approx <- rep(list(c(0, 0)), length(single_ma_data_set$TRT_event))
  ecc_random_approx <- rep(list(c(0, 0)), length(single_ma_data_set$TRT_event))
  
  for (i in 1:length(single_ma_data_set$TRT_event)) {
    e_fixed <- GetEmpiricalCCs(n_it = single_ma_data_set$TRT_n[i], 
                               n_ic = single_ma_data_set$CTRL_n[i], 
                               omega_hat = omega_hat_fixed, 
                               prop_constraint = prop_constraint)
    e_random <- GetEmpiricalCCs(n_it = single_ma_data_set$TRT_n[i], 
                                n_ic = single_ma_data_set$CTRL_n[i], 
                                omega_hat = omega_hat_random, 
                                prop_constraint = prop_constraint)
    ecc_fixed[[i]] <- e_fixed[[1]]
    ecc_random[[i]] <- e_random[[1]]
    ecc_fixed_approx[[i]] <- e_fixed[[2]]
    ecc_random_approx[[i]] <- e_random[[2]]
  }
  
  # Apply the corrections ------------------------------------------------------
  single_ma_data_set_fixed <- single_ma_data_set
  single_ma_data_set_random <- single_ma_data_set
  
  TCC_fixed <- sapply(ecc_fixed, function(y) y[1])
  CCC_fixed <- sapply(ecc_fixed, function(y) y[2])
  TCC_random <- sapply(ecc_random, function(y) y[1])
  CCC_random <- sapply(ecc_random, function(y) y[2])
  
  # If a study has no zero events (or zero non-events), then set the correction
  # to zero
  study_w_no_events <- single_ma_data_set$TRT_event != 0 & 
    single_ma_data_set$CTRL_event != 0 &
    TRT_no_event != 0 & 
    CTRL_no_event != 0
  TCC_fixed[study_w_no_events] <- 0
  CCC_fixed[study_w_no_events] <- 0
  TCC_random[study_w_no_events] <- 0
  CCC_random[study_w_no_events] <- 0
  
  # If ccs are less than zero, throw error
  if (sum(TCC_fixed < 0) > 0  | sum(CCC_fixed < 0) > 0 |
      sum(TCC_random < 0) > 0  | sum(CCC_random < 0) > 0) {
    stop("Empirical CC Corrections are Negative")
  } 
  
  single_ma_data_set_fixed$TRT_event <- single_ma_data_set$TRT_event + TCC_fixed
  single_ma_data_set_fixed$TRT_n <- single_ma_data_set$TRT_n + (2 * TCC_fixed)
  single_ma_data_set_fixed$CTRL_event <- single_ma_data_set$CTRL_event + 
    CCC_fixed
  single_ma_data_set_fixed$CTRL_n <- single_ma_data_set$CTRL_n + (2 * CCC_fixed)
  
  single_ma_data_set_random$TRT_event <- single_ma_data_set$TRT_event + 
    TCC_random
  single_ma_data_set_random$TRT_n <- single_ma_data_set$TRT_n + (2 * TCC_random)
  single_ma_data_set_random$CTRL_event <- single_ma_data_set$CTRL_event + 
    CCC_random
  single_ma_data_set_random$CTRL_n <- single_ma_data_set$CTRL_n + 
    (2 * CCC_random)
  
  return(list(single_ma_data_set_fixed, single_ma_data_set_random,
              sum(unlist(ecc_fixed_approx)), sum(unlist(ecc_random_approx))))
}

################################################################################
# 6. Exclude all single- and double-zero studies (Exclude_Zeros)
################################################################################
# Omits all studies containing single- and double-zero studies.
# This CC can only be used when 2 or more studies in the meta-analysis contain 
# at least one event in both study arms.

Exclude_Zeros <- function(single_ma_data_set) {
  
  TRT_no_event <- single_ma_data_set$TRT_n - single_ma_data_set$TRT_event
  CTRL_no_event <- single_ma_data_set$CTRL_n - single_ma_data_set$CTRL_event
  
  # Find double-zero studies ---------------------------------------------------
  # double-zero events
  double_zeros_loc <- (single_ma_data_set$TRT_event + 
                         single_ma_data_set$CTRL_event == 0) |
    # double-zero non-events
    (TRT_no_event + CTRL_no_event == 0) |
    # double-zero with no TRT event and no CTRL non-event
    (single_ma_data_set$TRT_event == 0 & CTRL_no_event == 0) |
    # double-zero with no CTRL event and no TRT non-event
    (single_ma_data_set$CTRL_event == 0 & TRT_no_event == 0)
  
  # Find single-zero studies ---------------------------------------------------
  # zero events
  single_zeros_loc <- ((single_ma_data_set$TRT_event == 0 | 
                          single_ma_data_set$CTRL_event == 0) |
                         # zero non-events
                         (TRT_no_event == 0 | 
                            CTRL_no_event == 0)) &
    # not a double-zero
    !double_zeros_loc
  
  # Exclude single- and double-zero studies ------------------------------------
  single_ma_data_set <- single_ma_data_set[!single_zeros_loc & 
                                             !double_zeros_loc, ]
  
  return(single_ma_data_set)
}
