# Apply Various Continuity Corrections, Heterogeneity Estimators, and Pooling
# Methods to the Anti-TNF Data Set

library(tidyverse)
library(readxl)
library(meta)

source("0_ApplyAllCCs.R")
source("0_MAResults.R")
set.seed(1234)

# read in data -----------------------------------------------------------------
single_ma_data_set <- read_excel("Anti-TNF.xlsx") %>% 
  rename(TRT_event = 'TRT e',
         TRT_n = 'TRT t',
         CTRL_event = 'CTRL e',
         CTRL_n = 'CTRL t')
single_ma_data_set <- as.data.frame(single_ma_data_set)

# Choose parameter values ------------------------------------------------------
pooling_method <- "SSW"  # "MH", "Inverse", "SSW", "SA"
tau2_estimator <- "SJ"  # "DL", "SJ", "PM", "IPM"
excludeDBZ <- TRUE  # TRUE, FALSE

# Create list to store results -------------------------------------------------
#   Each item in the list contains the results for the ith data set.
#   Each item has 12 MA output statistics (columns) for the 14 data set 
#   versions (rows)
ma_results <- matrix(ncol = 12, nrow = 14)
row.names(ma_results) <- c("original_data_set",
                           "cac_cc_1", "cac_cc_0.01", 
                           "czc_cc_1", "czc_cc_0.01", 
                           "cas_cc_1", "cas_cc_0.01",
                           "ta_cc_1", "ta_cc_0.01", 
                           "e_cc_1_fixed", "e_cc_1_random",
                           "e_cc_0.01_fixed", 
                           "e_cc_0.01_random",
                           "exclude_zeros_cc")
colnames(ma_results) <- c("TE.fixed", "upper.fixed",
                          "lower.fixed", "pval.fixed", 
                          "TE.random", "upper.random", 
                          "lower.random", "pval.random",
                          "tau2", "num.single.zeros",
                          "num.double.zeros", 
                          "e.cc.approx.used")

# Count the number of single- and double-zero studies in the original data -----
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
ma_results[, "num.double.zeros"] <- sum(double_zeros_loc)

# Single-zeros ---------------------------
# zero events
single_zeros_loc <- ((single_ma_data_set$TRT_event == 0 | 
                        single_ma_data_set$CTRL_event == 0) |
                       # zero non-events
                       (TRT_no_event == 0 | 
                          CTRL_no_event == 0)) &
  # not a double-zero
  !double_zeros_loc
ma_results[, "num.single.zeros"] <- sum(single_zeros_loc)


# Run data set through CC methods ----------------------------------------------

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

# Note the number of studies where the approximate solution of E_CC was used
ma_results["e_cc_1_fixed", "e.cc.approx.used"] <- 
  cc_results$e_cc_1_fixed_approx
ma_results["e_cc_1_random", "e.cc.approx.used"] <- 
  cc_results$e_cc_1_random_approx
ma_results["e_cc_0.01_fixed", "e.cc.approx.used"] <- 
  cc_results$e_cc_0.01_fixed_approx
ma_results["e_cc_0.01_random", "e.cc.approx.used"] <- 
  cc_results$e_cc_0.01_random_approx

# Run data set through MA methods ----------------------------------------------
# For MH, original data set can be analyzed without a CC applied
if (pooling_method == "MH") {
  
  for (i in 1:nrow(ma_results)) {
    single_ma_data_set_cc <- cc_results[[i]] 
    ma_results[i, 1:9] <- MAResults(single_ma_data_set = 
                                      single_ma_data_set_cc, 
                                    single_ma_data_set_uncorrected = 
                                      single_ma_data_set, 
                                    pooling_method = pooling_method, 
                                    tau2_estimator = tau2_estimator)
  }
  
} else {  # For other pooling methods, skip analyzing original data
  
  for (i in 1:(nrow(ma_results) - 1)) {  # -1 to skip original data
    # +1 to skip original data
    single_ma_data_set_cc <- cc_results[[i + 1]]  
    ma_results[i + 1, 1:9] <- MAResults(single_ma_data_set = 
                                          single_ma_data_set_cc, 
                                        single_ma_data_set_uncorrected = 
                                          single_ma_data_set, 
                                        pooling_method = 
                                          pooling_method, 
                                        tau2_estimator = 
                                          tau2_estimator)
  }
}

ma_results