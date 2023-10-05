# Apply Various Continuity Corrections to Individual MA Data Sets

# 1. Constant to All Cells CC (CAC_CC)
# 2. Constant to Zero Cells CC (CZC_CC)
# 3. Constant to All Studies CC (CAS_CC)
# 4. Treatment Arm CC (TA_CC)
# 5. Empirical CC (E_CC)
# 6. Exclude all single- and double-zero studies (Exclude_Zeros, Exclude_SZ_DZ)

source("0_CCCode.R")  # Contains the CAC_CC, CZC_CC, CAS_CC, TA_CC, E_CC, and 
# the Exclude_Zeros functions

ApplyAllCCs <- function(single_ma_data_set, pooling_method, tau2_estimator,
                        excludeDBZ) {
  
  # 1. Constant to All Cells CC (CAC_CC) ---------------------------------------
  cac_cc_1 <- CAC_CC(single_ma_data_set, cc = 0.5)
  cac_cc_0.01 <- CAC_CC(single_ma_data_set, cc = 0.005)
  
  # 2. Constant to Zero Cells CC (CZC_CC) --------------------------------------
  czc_cc_1 <- CZC_CC(single_ma_data_set, cc = 0.5)
  czc_cc_0.01 <- CZC_CC(single_ma_data_set, cc = 0.005)
  
  # 3. Constant to All Studies CC (CAS_CC) -------------------------------------
  cas_cc_1 <- CAS_CC(single_ma_data_set, cc = 0.5)
  cas_cc_0.01 <- CAS_CC(single_ma_data_set, cc = 0.005)
  
  # 4. Treatment Arm CC (TA_CC) ------------------------------------------------
  ta_cc_1 <- TA_CC(single_ma_data_set, prop_constraint = 1)
  ta_cc_0.01 <- TA_CC(single_ma_data_set, prop_constraint = 0.01)
  
  # 5. Empirical CC (E_CC) -----------------------------------------------------
  # Returns a list of two data sets: one based on a fixed-effect prior log OR, 
  #   one based on a random-effects prior log OR
  e_cc_1 <- E_CC(single_ma_data_set, 
                 prop_constraint = 1, 
                 tau2_estimator,
                 pooling_method)
  e_cc_0.01 <- E_CC(single_ma_data_set, 
                    prop_constraint = 0.01, 
                    tau2_estimator,
                    pooling_method)
  
  # 6. Exclude all single- and double-zero studies (Exclude_Zeros) -------------
  exclude_zeros_cc <- Exclude_Zeros(single_ma_data_set)
  
  return(list("original_data_set" = single_ma_data_set,
              "cac_cc_1" = cac_cc_1, 
              "cac_cc_0.01" = cac_cc_0.01, 
              "czc_cc_1" = czc_cc_1, 
              "czc_cc_0.01" = czc_cc_0.01, 
              "cas_cc_1" = cas_cc_1,
              "cas_cc_0.01" = cas_cc_0.01, 
              "ta_cc_1" = ta_cc_1, 
              "ta_cc_0.01" = ta_cc_0.01, 
              "e_cc_1_fixed" = e_cc_1[[1]],
              "e_cc_1_random" = e_cc_1[[2]],
              "e_cc_0.01_fixed" = e_cc_0.01[[1]],
              "e_cc_0.01_random" = e_cc_0.01[[2]],
              "exclude_zeros_cc" = exclude_zeros_cc,
              # Return number of studies where approximate solution is used
              "e_cc_1_fixed_approx" = e_cc_1[[3]],
              "e_cc_1_random_approx" = e_cc_1[[4]],
              "e_cc_0.01_fixed_approx" = e_cc_0.01[[3]],
              "e_cc_0.01_random_approx" = e_cc_0.01[[4]]))
}