
##########################################
# n: number of patients
# p: number of clinical endpoints
# k: number of factors
# TRT_prob: probability to randomly generate treatment, where 1 is active and 0 is placebo  
# time: Time vector
# A_prob: probability to randomly sample from (1,0) for A matrix
# E0_null: E0 population level NULL parameter 
# E50_null: E50 population level NULL parameter 
# Emax_null: Emax population level NULL parameter 
# beta_E0:  treatment effect in E0
# beta_E50: treatment effect in E50
# beta_Emax: treatment effect in Emax
# theta_RI_sd: random intercept theta_RI sd
# theta_resid_sd: random residual sd
# y_resid: y residual sd

# 
# seed = it
# n = n_subj
# p = n_var
# k = n_fac
# TRT_prob = TRT_prob
# endpoint_proportions = endpoint_proportions
# time = weeks_measure
# cor_ar1_lag = cor_ar1_lag
# E0_null = E0_null
# E50_null = E50_null
# Emax_null = Emax_null
# beta_E0 = beta_E0
# beta_E50 = beta_E50
# beta_Emax = beta_Emax
# theta_RI_sd = theta_RI_sd
# theta_resid_sd = theta_resid_sd
# y_resid_sd = y_resid_sd

simulate_data <- function(seed = 1,
                      n = 200,
                      p = 10,
                      k = 5,
                      TRT_prob = 0.5,
                      endpoint_proportions = c(0.2, 0.3, 0.2, 0.1, 0.2),
                      time = c(0, 4, 8, 12, 16),
                      cor_ar1_lag = 0.2,
                      E0_null = 0,
                      E50_null = 7,
                      Emax_null = 3,
                      beta_E0 = 0.25,
                      beta_E50 = -0.9,
                      beta_Emax = 1.5,
                      theta_RI_sd = 1,
                      theta_resid_sd = 1,
                      y_resid_sd = 1) {

  require(dplyr)
  set.seed(seed)
  ## treatment
  fixed_proportion_set <- c(rep(0, floor(TRT_prob * n)), rep(1, n - floor(TRT_prob * n)))
  trt <- sample(fixed_proportion_set, n, replace = FALSE)
  trt_char <- ifelse(trt == 0, "Placebo", "Active")
  # loading matrix A
  # Constructing block-wise A, where 
  # each entry is either 0 or 1, and
  # each endpoint is 1 ("active") in exactly one factor
  if (sum(endpoint_proportions) != 1) {
    stop("'endpoint_proportions' must sum to one")
  }
  if (length(endpoint_proportions) != k) {
    stop("length of 'endpoint_proportions' must equal k")
  }
  
  
  intervals <- c(0, cumsum(endpoint_proportions))
  A_loadings <- matrix(0, nrow = k, ncol = p)
  for (fac_curr in 1:n_fac) {
    which_fill_in <- which((1:p / p) > intervals[fac_curr] & (1:p / p) <= intervals[fac_curr + 1])
    A_loadings[fac_curr, which_fill_in] <- 1
  }

  rownames(A_loadings) <- paste0("K", c(1:k))
  colnames(A_loadings) <- paste0("P", c(1:p))
  
  # dataset with patient ID and treatment
  dbase <- data.frame(ID = c(1:n),
                      TRT = trt_char)
  

  ## CAR1
  cs1AR1 <- nlme::corAR1(cor_ar1_lag, form = ~ week)
  d_demo <- data.frame(week = time)
  cs1AR1. <- nlme::Initialize(cs1AR1, data = d_demo)
  CAR1_theta <- as.matrix(nlme::corMatrix(cs1AR1.))
  
    
  #simulate Z matrix
  Z_matrix <- NULL
  # K_th <- 1
  for (K_th in 1:k){
    
    # E0 is equal to E0_null + trt effect if active arm and plus random intercept
    Z_E0 <- E0_null * (1 + t(beta_E0 %*% t(as.matrix(trt)))) + rnorm(n, mean = 0, sd = theta_RI_sd)
    
    # E50 and Emax does not have random error
    Z_E50 <- E50_null * (1 + t(beta_E50 %*% t(as.matrix(trt))))  
    Z_Emax <- Emax_null * (1 + t(beta_Emax %*% t(as.matrix(trt))))  
    
    Z_kth <- NULL
    for (t_th in unique(time)){
      Z_nt <- Z_E0 + (Z_Emax * (t_th)) / (Z_E50 + t_th) #+ rnorm(n, mean = 0, sd = theta_resid_sd)
      Z_nt_tmp <- data.frame(Z = Z_nt,
                             Time = t_th,
                             K = paste0("Z", K_th),
                             ID = 1:n,
                             TRT = trt_char,
                             stringsAsFactors = FALSE)
      Z_kth <- rbind(Z_kth, Z_nt_tmp)
    }
    
    Z_kth2 <- Z_kth %>% tidyr::spread(Time, Z) 
    # test <- (Z_kth2 %>% select(!c(K, ID, TRT)) %>% as.matrix())
    # test - test[, 1]
    # Z_kth3 <- (Z_kth2 %>% select(!c(K, ID, TRT)) %>% as.matrix()) %*% CAR1_theta
    Z_kth3 <- (Z_kth2 %>% select(!c(K, ID, TRT)) %>% as.matrix())
    if (theta_resid_sd > 0) {
      Z_kth3 <- Z_kth3 + mvnfast::rmvn(n = n, mu = rep(0, length(time)), sigma = CAR1_theta * theta_resid_sd^2)
    }
    colnames(Z_kth3) <- time
    
    Z_kth_final <- reshape2::melt(cbind(Z_kth2 %>% select(c(K, ID, TRT)), Z_kth3), id = c("K", "ID", "TRT")) %>% 
      dplyr::rename(Time = variable, Z = value)
    Z_matrix <- rbind(Z_matrix, Z_kth_final)
  }
  
  Z <- Z_matrix %>% tidyr::spread(K, Z)
  
  # simulate Y matrix
  Z_scores <- Z %>% select(!c("ID", "Time", "TRT"))
  
  Y_resid <- matrix(data = rnorm(n = n * p * length(time), mean = 0, sd = y_resid_sd), 
                    nrow = n * length(time), 
                    ncol = p)
  
  Y_scores <- as.matrix(Z_scores) %*% (A_loadings) + Y_resid
  
  # A_loadings
  # image(cor(Y_scores, Z_scores))
  colnames(Y_scores) <- paste0("Y", 1:p)
  Y <- cbind(Z %>% select(c("ID", "Time", "TRT")), Y_scores)
  Y$Time <- as.numeric(as.character(Y$Time))
  
  # Calc expected values of Z and Y for treatment and placebo groups
  E_Z <- matrix(NA, length(time), 2, dimnames = list(time, paste0("trt_", 0:1)))
  E_Y <- array(NA, dim = c(length(time), p, 2), dimnames = list(time, paste0("Y", 1:p), paste0("trt_", 0:1)))
  for (trt_curr in 0:1) {
    trt_nam <- paste0("trt_", trt_curr)
    Z_E0 <- E0_null * (1 + beta_E0 * trt_curr)
    Z_E50 <- E50_null * (1 + beta_E50 * trt_curr)  
    Z_Emax <- Emax_null * (1 + beta_Emax * trt_curr)
    E_Z[, trt_nam] <- Z_E0 + (Z_Emax * (time)) / (Z_E50 + time)
    mean_scores <- NULL
    for (kc in 1:k) {
      mean_scores <- cbind(mean_scores, E_Z[, trt_nam])
    }
    E_Y[, , trt_nam] <- mean_scores  %*% A_loadings
  }

  
  output <- list(A_loadings = A_loadings, 
                 Z = Z, 
                 Y = Y, 
                 E_Z = E_Z,
                 E_Y = E_Y,
                 time = time,
                 dbase = dbase,
                 parameters = c(seed = seed, 
                                n = n, 
                                p = p,
                                k = k,
                                TRT_prob = TRT_prob,
                                E0_null = E0_null,
                                E50_null = E50_null,
                                Emax_null = Emax_null,
                                beta_E0 = beta_E0,
                                beta_E50 = beta_E50,
                                beta_Emax = beta_Emax))
  
  return(output)

}
