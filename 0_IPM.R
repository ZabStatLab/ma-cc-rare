# Code to Implement the Improved Paule-Mandel (IPM) Heterogeneity Variance
# Estimator

# Optimization function for IPM tau2 estimator ---------------------------------
QGen <- function(tau2, mu, theta_i, theta_SA, n_it, n_ic, k) {
  var_i_IPM <- ((exp(-mu - theta_SA + (tau2 / 2)) + 2 +
                   exp(mu + theta_SA + (tau2 / 2))) / 
                  (n_it + 1)) + ((exp(-mu) + 2 + exp(mu)) / (n_ic + 1))
  ss_weights <- 1 / (var_i_IPM + tau2)
  theta_IPM <- sum(theta_i * ss_weights) / sum(ss_weights)
  sum(ss_weights * (theta_i - theta_IPM)^2) - (k - 1)
}

# Get IPM tau2 estimate -------------------------------------------------------
IPM <- function(single_ma_data_set_uncorrected, maxit = 100, signiftau2 = 6) {
  
  # Add 0.5 to all cells of uncorrected data
  cc <- 0.5
  a <- single_ma_data_set_uncorrected$TRT_event + cc
  b <- single_ma_data_set_uncorrected$TRT_n - 
    single_ma_data_set_uncorrected$TRT_event + cc
  c <- single_ma_data_set_uncorrected$CTRL_event + cc
  d <- single_ma_data_set_uncorrected$CTRL_n - 
    single_ma_data_set_uncorrected$CTRL_event + cc
  n_it <- single_ma_data_set_uncorrected$TRT_n + (2 * cc)
  n_ic <- single_ma_data_set_uncorrected$CTRL_n + (2 * cc)
  k <- length(single_ma_data_set_uncorrected$TRT_event)
  
  mu_i <- log(c / d)  # log odds in control group
  mu <- mean(mu_i)  # average log odds in control group
  # Note that mu is not defined in Bhaumik et al. Using mu = mu_i resulted in a
  # low estimate compared to what Bhaumik et al report from their example data.
  # Using mu = mean(mu_i) resulted in an estimate much closer to what they
  # report.
  theta_i <- log((a * d) / (b * c))  # study-specific log odds ratios
  theta_SA <- sum(theta_i) / k  # combined log odds ratio (unweighted average)
  tau2 <- 0  # initialize

  # IPM variance (Equation 14):
  var_i_IPM <- ((exp(-mu - theta_SA + (tau2 / 2)) + 2 +
                   exp(mu + theta_SA + (tau2 / 2))) / 
                  (n_it + 1)) + ((exp(-mu) + 2 + exp(mu)) / (n_ic + 1))
  # in text below Equation 14 (ss = shared-strength):
  ss_weights <- 1 / (var_i_IPM + tau2)
  theta_IPM <- sum(theta_i * ss_weights) / sum(ss_weights)  # Equation 7
  Q_gen <- sum(ss_weights * (theta_i - theta_IPM)^2)  # part of Equation 13
  
  if (Q_gen < (k - 1)) {
    
    tau2 <- 0
    
  } else {
    
    tau2 <- uniroot(QGen, 
                    interval = c(0, 100), 
                    mu = mu, 
                    theta_i = theta_i, 
                    theta_SA = theta_SA,
                    n_it = n_it, 
                    n_ic = n_ic, 
                    k = k,
                    maxiter = 100)$root
    
  }
  
  return(tau2)
}