# Apply Various Continuity Corrections, Heterogeneity Estimators, and Pooling
# Methods to the Anti-TNF Data Set

library(tidyverse)
library(readxl)
library(meta)
library(metafor)

source("0_ApplyAllCCs.R")
source("0_IPM.R")

settings.meta(digits = 2)

# read in data -----------------------------------------------------------------
anti_tnf <- read_excel("Anti-TNF.xlsx") %>% 
  rename(TRT_event = 'TRT e',
         TRT_n = 'TRT t',
         CTRL_event = 'CTRL e',
         CTRL_n = 'CTRL t')
anti_tnf <- as.data.frame(anti_tnf)

################################################################################
#                    Results reported in publication
################################################################################
# apply all correction methods to the data
cc_results <- ApplyAllCCs(single_ma_data_set = anti_tnf, 
                          pooling_method = "SA", 
                          tau2_estimator = "DL",
                          excludeDBZ = TRUE)

# Section 5 --------------------------------------------------------------------

# common-effect framework, MH pooling method, TA-CC:
m_MH <- metabin(TRT_event, 
                TRT_n,
                CTRL_event,
                CTRL_n,
                sm = "OR",
                method = "MH",
                data = anti_tnf,
                incr = "TACC",
                random = FALSE)
m_MH  

# common-effect framework, MH pooling method, MH-No-CC:
m_MH_noCC <- update(m_MH, incr = 0)
m_MH_noCC

# common-effect framework, SA pooling method, E-CC:
m_SA_ECC <- rma.uni(ai = TRT_event,
                    n1i = TRT_n,
                    ci = CTRL_event,
                    n2i = CTRL_n,
                    measure = "OR",
                    data = cc_results$e_cc_1_fixed,
                    add = 0,
                    method = "FE",
                    weighted = FALSE)
c(exp(m_SA_ECC$beta), exp(m_SA_ECC$ci.lb), exp(m_SA_ECC$ci.ub), m_SA_ECC$pval)

# common-effect framework, SSW pooling method, CAC-CC with k=0.005:
m_SSW_CACCC <- update(m_MH, 
                      method = "SSW",
                      incr = 0.005, 
                      allincr = FALSE, 
                      addincr = FALSE)
m_SSW_CACCC

# random-effects framework, SSW pooling method, TA-CC, SJ heterogeneity est.: 
m_SSW_TACC <- metabin(TRT_event, 
                      TRT_n,
                      CTRL_event,
                      CTRL_n,
                      sm = "OR",
                      method = "SSW",
                      method.tau = "SJ",
                      data = cc_results$ta_cc_1,
                      incr = 0,
                      common = FALSE)
m_SSW_TACC

# Table 7 ----------------------------------------------------------------------
m_DL <- metabin(TRT_event, 
                TRT_n,
                CTRL_event,
                CTRL_n,
                sm = "OR",
                method = "MH",
                method.tau = "DL",
                data = anti_tnf,
                incr = "TACC",
                common = FALSE)

m_SJ <- update(m_DL, method.tau = "SJ")
m_PM <- update(m_DL, method.tau = "PM")
m_IPM <- update(m_DL, tau.preset = sqrt(IPM(anti_tnf)))

m_IPM
m_PM
m_DL
m_SJ
