
simulate_data <- function(control,
                           seed = 1,
                            zero_ATE = FALSE) {

  require(dplyr)
  set.seed(seed)

  fixed_proportion_set <- c(rep(0, floor(control$TRT_prob * control$n_subj)), 
                            rep(1, control$n_subj - floor(control$TRT_prob * control$n_subj)))
  trt <- sample(fixed_proportion_set, control$n_subj, replace = FALSE)
  trt_char <- ifelse(trt == 0, "Placebo", "Active")
  # loading matrix A
  # Constructing block-wise A, where 
  # each entry is either 0 or 1, and
  # each endpoint is 1 ("active") in exactly one factor
  if (sum(control$endpoint_proportions) != 1) {
    stop("'control$endpoint_proportions' must sum to one")
  }
  if (length(control$endpoint_proportions) != control$n_fac) {
    stop("length of 'control$endpoint_proportions' must equal control$n_fac")
  }
  
  # Create loadings matrix or use loadings from FUTURE-2 
  if (control$use_empirical_loadings) {
    A_loadings <- readRDS("data/empirical_loadings_FUTURE-2.RDS")
  } else {
    intervals <- c(0, cumsum(control$endpoint_proportions))
    A_loadings <- matrix(0, nrow = control$n_fac, ncol = control$n_var)
    for (fac_curr in 1:control$n_fac) {
      which_fill_in <- which((1:control$n_var / control$n_var) > intervals[fac_curr] & 
                               (1:control$n_var / control$n_var) <= intervals[fac_curr + 1])
      A_loadings[fac_curr, which_fill_in] <- 1
    }
  }

  rownames(A_loadings) <- paste0("K", c(1:control$n_fac))
  colnames(A_loadings) <- paste0("P", c(1:control$n_var))
  
  # dataset with patient ID and treatment
  dbase <- data.frame(ID = c(1:control$n_subj),
                      TRT = trt_char)
  
  ## CAR1
  cs1AR1 <- nlme::corAR1(control$cor_ar1_lag, form = ~ week)
  d_demo <- data.frame(week = control$weeks_measure)
  cs1AR1. <- nlme::Initialize(cs1AR1, data = d_demo)
  CAR1_theta <- as.matrix(nlme::corMatrix(cs1AR1.))

  Z_matrix <- NULL
  # K_th <- 1
  for (K_th in 1:control$n_fac){
    size_of_final_week_effect <- ifelse(trt_char == "Placebo", 
                                     control$size_of_final_week_effect_placebo,
                                     ifelse(zero_ATE, control$size_of_final_week_effect_placebo, control$size_of_final_week_effect_active))
    random_intercept <- rnorm(control$n_subj, mean = 0, sd = control$theta_RI_sd)
    Z_kth <- NULL
    for (t_th in unique(control$weeks_measure)){
      Z_nt <- random_intercept + emax_curve(t_eval = t_th, 
                                               time_of_final_week_effect = control$time_pt_eval, 
                                               size_of_final_week_effect = size_of_final_week_effect,
                                               time_of_half_max_effect = control$time_of_half_max_effect
                                                )
      Z_nt_tmp <- data.frame(Z = Z_nt,
                             Time = t_th,
                             K = paste0("Z", K_th),
                             ID = 1:control$n_subj,
                             TRT = trt_char,
                             stringsAsFactors = FALSE)
      Z_kth <- rbind(Z_kth, Z_nt_tmp)
    }
    Z_kth2 <- Z_kth %>% tidyr::spread(Time, Z) 
    Z_kth3 <- (Z_kth2 %>% select(!c(K, ID, TRT)) %>% as.matrix())
    if (control$theta_resid_sd > 0) {
      Z_kth3 <- Z_kth3 + mvnfast::rmvn(n = control$n_subj, 
                                       mu = rep(0, length(control$weeks_measure)), 
                                       sigma = CAR1_theta * control$theta_resid_sd^2)
    }
    colnames(Z_kth3) <- control$weeks_measure
    Z_kth_final <- reshape2::melt(cbind(Z_kth2 %>% select(c(K, ID, TRT)), Z_kth3), id = c("K", "ID", "TRT")) %>% 
      dplyr::rename(Time = variable, Z = value)
    Z_matrix <- rbind(Z_matrix, Z_kth_final)
  }
  # Calc expected values of Z and Y for treatment and placebo groups
  E_Z <- matrix(NA, length(control$weeks_measure), 2, dimnames = list(control$weeks_measure, paste0("trt_", 0:1)))
  E_Y <- array(NA, dim = c(length(control$weeks_measure), control$n_var, 2), 
               dimnames = list(control$weeks_measure, paste0("Y", 1:control$n_var), paste0("trt_", 0:1)))
  for (trt_curr in 0:1) {
    trt_nam <- paste0("trt_", trt_curr)
    E_Z[, trt_nam] <- emax_curve(t_eval = control$weeks_measure, 
                  time_of_final_week_effect = control$time_pt_eval, 
                  size_of_final_week_effect = ifelse(trt_curr == 0 | zero_ATE, control$size_of_final_week_effect_placebo,
                                                     control$size_of_final_week_effect_active),
                  time_of_half_max_effect = control$time_of_half_max_effect
            )
    mean_scores <- NULL
    for (kc in 1:control$n_fac) {
      mean_scores <- cbind(mean_scores, E_Z[, trt_nam])
    }
    E_Y[, , trt_nam] <- mean_scores %*% A_loadings
  }

  Z <- Z_matrix %>% tidyr::spread(K, Z)
  
  # simulate Y matrix
  Z_scores <- Z %>% select(!c("ID", "Time", "TRT"))
  Y_resid <- matrix(data = rnorm(n = control$n_subj * control$n_var * length(control$weeks_measure), 
                                 mean = 0, 
                                 sd = control$y_resid_sd), 
                    nrow = control$n_subj * length(control$weeks_measure), 
                    ncol = control$n_var)
  Y_scores <- as.matrix(Z_scores) %*% (A_loadings) + Y_resid
  colnames(Y_scores) <- paste0("Y", 1:control$n_var)
  Y <- cbind(Z %>% select(c("ID", "Time", "TRT")), Y_scores)
  Y$Time <- as.numeric(as.character(Y$Time))
  
  output <- list(A_loadings = A_loadings, 
                 Z = Z, 
                 Y = Y, 
                 E_Z = E_Z,
                 E_Y = E_Y,
                 time = control$weeks_measure,
                 dbase = dbase,
                 seed = seed, 
                 zero_ATE = zero_ATE,
                 control = control)
  return(output)
}
  
#   
#                                 n = control$n_subj, 
#                                 p = control$n_var,
#                                 k = control$n_fac,
#                                 TRT_prob = control$TRT_prob,
#                                 E0_null = control$E0_null,
#                                 E50_null = control$E50_null,
#                                 Emax_null = control$Emax_null,
#                                 beta_E0 = control$beta_E0,
#                                 beta_E50 = control$beta_E50,
#                                 beta_Emax = control$beta_Emax))
#   
#   return(output)
# 
# }
#                       n = 200,
#                       p = 10,
#                       k = 5,
#                       TRT_prob = 0.5,
#                       endpoint_proportions = c(0.2, 0.3, 0.2, 0.1, 0.2),
#                       time = c(0, 4, 8, 12, 16),
#                       cor_ar1_lag = 0.2,
#                       theta_RI_sd = 1,
#                       ATE = 1,
#                       E0_null = 0,
#                       E50_null = 7,
#                       Emax_null = 3,
#                       beta_E0 = 0.25,
#                       beta_E50 = -0.9,
#                       beta_Emax = 1.5,
#                       sim_type = c("orig_emax", "simplified_emax")[2]
# ) {
# 
# 
# if (sim_type == "orig_emax") {
#   Z_matrix <- NULL
#   # K_th <- 1
#   for (K_th in 1:control$n_fac){
#     # E0 is equal to E0_null + trt effect if active arm and plus random intercept
#     Z_E0 <- control$E0_null * (1 + t(control$beta_E0 %*% t(as.matrix(trt)))) + rnorm(control$n_subj, mean = 0, sd = control$theta_RI_sd)
#     # E50 and Emax does not have random error
#     Z_E50 <- control$E50_null * (1 + t(control$beta_E50 %*% t(as.matrix(trt))))  
#     Z_Emax <- control$Emax_null * (1 + t(control$beta_Emax %*% t(as.matrix(trt))))  
#     Z_kth <- NULL
#     for (t_th in unique(control$weeks_measure)){
#       Z_nt <- Z_E0 + (Z_Emax * (t_th)) / (Z_E50 + t_th) #+ rnorm(n, mean = 0, sd = theta_resid_sd)
#       Z_nt_tmp <- data.frame(Z = Z_nt,
#                              Time = t_th,
#                              K = paste0("Z", K_th),
#                              ID = 1:control$n_subj,
#                              TRT = trt_char,
#                              stringsAsFactors = FALSE)
#       Z_kth <- rbind(Z_kth, Z_nt_tmp)
#     }
#     Z_kth2 <- Z_kth %>% tidyr::spread(Time, Z) 
#     Z_kth3 <- (Z_kth2 %>% select(!c(K, ID, TRT)) %>% as.matrix())
#     if (control$theta_resid_sd > 0) {
#       Z_kth3 <- Z_kth3 + mvnfast::rmvn(n = control$n_subj, mu = rep(0, length(control$weeks_measure)), sigma = CAR1_theta * control$theta_resid_sd^2)
#     }
#     colnames(Z_kth3) <- control$weeks_measure
#     Z_kth_final <- reshape2::melt(cbind(Z_kth2 %>% select(c(K, ID, TRT)), Z_kth3), id = c("K", "ID", "TRT")) %>% 
#       dplyr::rename(Time = variable, Z = value)
#     Z_matrix <- rbind(Z_matrix, Z_kth_final)
#   }
#   # Calc expected values of Z and Y for treatment and placebo groups
#   E_Z <- matrix(NA, length(control$weeks_measure), 2, dimnames = list(control$weeks_measure, paste0("trt_", 0:1)))
#   E_Y <- array(NA, dim = c(length(control$weeks_measure), control$n_var, 2), dimnames = list(control$weeks_measure, paste0("Y", 1:control$n_var), paste0("trt_", 0:1)))
#   for (trt_curr in 0:1) {
#     trt_nam <- paste0("trt_", trt_curr)
#     Z_E0 <- control$E0_null * (1 + control$beta_E0 * trt_curr)
#     Z_E50 <- control$E50_null * (1 + control$beta_E50 * trt_curr)  
#     Z_Emax <- control$Emax_null * (1 + control$beta_Emax * trt_curr)
#     E_Z[, trt_nam] <- Z_E0 + (Z_Emax * (control$weeks_measure)) / (Z_E50 + control$weeks_measure)
#     mean_scores <- NULL
#     for (kc in 1:control$n_fac) {
#       mean_scores <- cbind(mean_scores, E_Z[, trt_nam])
#     }
#     E_Y[, , trt_nam] <- mean_scores  %*% A_loadings
#   }
# }