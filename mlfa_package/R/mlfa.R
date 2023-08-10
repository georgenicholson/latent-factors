# make non-random
# set.seed(1)


#' @title
#' Multivariate Longitudinal Factor Analysis.
#' 
#' @description 
#' Fit \code{gamm} models marginally for each factor (stage 3) that allow to predict longitudinal latent trajectories.
#' 
#' @details
#' This function performs stage 3 of the proposed method. This corresponds to the
#' functionality of fitting \code{gamm} models of a certain type for all factors, using the package \emph{mgcv} package.
#'
#' This function returns an \code{mlfa} object that has various attributes. The most important result attribute is:
#' \itemize{
#'   \item \code{mlfa$gamm_objects}: A list of \code{gamm} objects, where each \code{gamm} object has two main attributes: \code{gamm_object$gam} and \code{gamm_object$lme},
#'                      corresponding to the two fitted sub-models of the \code{gamm(...)} function call.
#' }
#'
#' Furthermore, there are several helper functions for an mlfa object that allow the following functionality:
#' \itemize{
#'   \item \code{plot.mlfa(...)}: This function predicts and plots the latent trajectories of all fitted gamm models, implicitly only using the fixed components.
#'                   The latent trajectories are plotting the fit over time, inferring time from what was passed to the \code{mfd} call in stage 1 and 2.
#'   \item \code{summary.mlfa(...)}: The function concatenates the summaries of the gam and lme models fitted in all gamm models (see \code{$gamm_objects}).
#' }
#' 
#' Further extensions:
#' \itemize{
#'   \item parallelize the fitting of the per-factor \code{gamm} models
#'   \item possible: have mfd$data_out / data as input instead of mfd -> stage 3 decoupled
#' }
#'
#' @param mfd \code{mfd} object; A fitted mfd object (result object of stage 1 and 2).
#' @param type int; The type of gamm model to fit. Default: 1 (type 1).
#' @param df int; The number of degrees of freedom of the fitted splines. The higher \code{df}, the wigglier the fitted spline.
#' @param subject_data tibble; A tibble with a column \code{mfd$subject} and various other columns that contain subject specific information.
#'                             This information can be used in the model fitting together with the results estimated in stage 1 and 2 (e.g. the scores \code{z}).
#'                             Default: \code{NULL} (not merging any data).
#' @param type1.strata character vector; The strata variable(s) of the type 1 \code{gamm} model.
#' @param type1.group string; The group variable of the type 1 \code{gamm} model. Usually \code{mfd$subject}.
#' @param type2.strata character vector; The strata variable(s) of the type 2 \code{gamm} model.
#' @param type2.group string; The group variable of the type 2 \code{gamm} model. Usually \code{mfd$subject}.
#' @return \code{mlfa} object; an object of S3 class \code{mlfa}.
#'               # Primary (necessary inputs)


# 
# mfd <- mfd_object0
# type <- 2
# df <- 5
# subject_data = d_subject_sub
# type1.strata = NULL
# type1.group = NULL
# type2.strata = "TRT01P"
# type2.group = "USUBJID"
# type3.strata = "newTRT"
# type3.group = "USUBJID"
# 

mlfa <- function(mfd, type = 1, df,
                 # Secondary inputs
                 subject_data = NULL,
                 arm_label = "newTRT",
                 placebo_label = "Placebo",
                 # type 1
                 type1.strata = NULL, type1.group = NULL,
                 # type 2
                 type2.strata = NULL, type2.group = NULL,
                 # type 3
                 type3.strata = NULL, type3.group = NULL) {
  
  mlfa <- structure(list(), class = "mlfa")
  
  # must assign NULL to attribute first
  attr(mlfa, "mfd") <- NULL
  mlfa$mfd <- mfd
  attr(mlfa, "data") <- NULL
  attr(mlfa, "gamm_objects") <- NULL
  mlfa$gamm_objects <- list()
  attr(mlfa, "fac_dim") <- NULL
  mlfa$fac_dim <- mfd$fac_dim
  attr(mlfa, "type") <- NULL
  mlfa$type <- type
  
  # storing 2 attributes from mfd here for convenience
  attr(mlfa, "time") <- NULL
  mlfa$time <- mfd$time
  attr(mlfa, "subject") <- NULL
  mlfa$subject <- mfd$subject
  
  # type 1
  if (type == 1) {
    attr(mlfa, "type1.strata") <- NULL
    attr(mlfa, "type1.group") <- NULL
  } else if (type == 2) {
    attr(mlfa, "type2.strata") <- NULL
    attr(mlfa, "type2.group") <- NULL
  } else if (type == 3) {
    attr(mlfa, "type3.strata") <- NULL
    attr(mlfa, "type3.group") <- NULL
  }
  
  if (!is.null(subject_data)) {
    # merge subject_data with data_out of mfd object
    # only merge columns of subject_data which are not already in mfd$data_out
    not_merged_cols <- setdiff(colnames(mfd$data_out), c(mfd$time, mfd$subject, colnames(mfd$A), colnames(mfd$Z), colnames(mfd$S)))  # but time and subject must be kept
    subject_data_temp <- subject_data %>% dplyr::select(-dplyr::one_of(not_merged_cols))
    mlfa$data <- dplyr::inner_join(mfd$data_out, subject_data_temp, by=c(mfd$subject))
  } else {
    mlfa$data <- mfd$data_out
  }
  
  # Fitting of stage 3 models
  # i=1
  print("Fit latent factor trajectories (stage 3) started...")
  
  
  
  
  for (i in 1:mlfa$mfd$fac_dim) {
    print(paste0("-> Fitting factor ", as.character(i), " of ", as.character(mlfa$mfd$fac_dim), " started..."))
    if (type == 1) {
      mlfa$type1.strata <- type1.strata
      mlfa$type1.group <- type1.group
      mlfa$gamm_objects[[i]] <- mlfa_type1(data = mlfa$data, df = df, factor = i, time = mlfa$time, strata = mlfa$type1.strata, group = mlfa$type1.group)
    }
    if (type == 2) {
      mlfa$type2.strata <- type2.strata
      mlfa$type2.group <- type2.group
      mlfa$gamm_objects[[i]] <- mlfa_type2(data = mlfa$data,
                                           df = df,
                                           factor = i,
                                           time = mlfa$time,
                                           strata = mlfa$type2.strata,
                                           group = mlfa$type2.group)
    }
    if (type == 3) {
      mlfa$type3.strata <- type3.strata
      mlfa$type3.group <- type3.group
      mlfa$gamm_objects[[i]] <- mlfa_type3(mlfa$data,
                                           df = df,
                                           factor = i,
                                           time = mlfa$time,
                                           strata = mlfa$type3.strata,
                                           group = mlfa$type3.group,
                                           arm_label = arm_label,
                                           placebo_label = placebo_label,
                                           m_smooth = 2)
    }
    print(paste0("-> Fitting factor ", as.character(i), " of ", as.character(mfd$fac_dim), " done."))
  }
  print("Fit latent factor trajectories (stage 3) done.")
  
  return(mlfa)
}


#' @title
#' Fit a \emph{Gaussian Additive Mixture Model} (type 2).
#' 
#' @description
#' This helper function wraps a common longitudinal model, especially for biomedical data, described in details.
#' 
#' @details
#' The underlying mgcv call is as follows:
#'
#' \code{
#' gamm(   
#'   formula = z ~ s(time, k = df) + strata[1] + s(time, k = df, by=factor(strata[1], ordered=TRUE)) + strata[2] + ... ,
#'   random = list(group=~1),
#'   correlation=corCAR1(form=~time|group),
#'   weights=varFixed(~s^2),
#'   data = data
#' )
#' }
#'
#' Possible variables are: 
#' \itemize{
#'    \item time: the timesteps, e.g. in weeks
#'    \item strata: the trial variable
#'    \item group: the subject variable
#' }
#' 
#' Possible extensions: 
#' \itemize{
#'    \item extend strata to character vector -> currently causes issue with corCAR
#' }
#'
#' @param data tibble; A tibble data frame where each row contains an observation
#'                        of multiple variables from a certain series (e.g.
#'                        a subject) and at a certain timestep (e.g. week). The columns
#'                        contain the series, the timestep and arbitrarily many
#'                        other variables.
#' @param df integer; The degrees of freedom of the spline bases.
#' @param factor integer; The factor (column of \code{Z} and \code{S}) chosen that shall be regressed on.
#' @param time string: The column name in \code{data} of the longitudinal variable chosen.
#' @param strata string: The column name in \code{data} that specify the factor for which an intercept and spline is estimated.  
#' @param group string: The column name in \code{data} that specify the factor for which a random intercept is estimated.
#' @return gamm_object \code{gamm} object; a \code{gamm} object as reutrn by the \code{gamm(...)} call.
#' 
#' 
#' 


mlfa_type2 <- function(data, df, factor = 1, time, strata, group) {
  formula <- paste0(paste0("z", as.character(factor), " ~ ", "s(", time, ", k = ", as.character(df), ")"))
  for (i in 1:length(strata)) {
    formula <- paste0(formula, " + ", 
                      strata[i], " + ", 
                      "s(", time, ", k = ", as.character(df), ", by = factor(", strata[i], ", ordered = TRUE))")
  }
  formula <- stats::as.formula(formula)
  formula_varFixed <- stats::as.formula(paste0("~ ", "s", as.character(factor), "^2"))
  formula_corCar1 <- stats::as.formula(paste0("~ ", time, " | ", group))
  random_list <- list()
  random_list[[group]] <- ~1
  
  gamm_object <- mgcv::gamm(   # !!(sym(paste("z", as.character(PC_sel), sep = "")))
    formula = formula,
    random = random_list, # sym(group)
    correlation=nlme::corCAR1(form=formula_corCar1),
    weights=nlme::varFixed(formula_varFixed),
    data = data
  )
  
  # stored for plotting
  gamm_object$time <- time
  gamm_object$strata <- strata
  gamm_object$group <- group
  gamm_object$factor <- factor
  
  return(gamm_object)
}


# 
# 
# mlfa$type3.strata <- type3.strata
# mlfa$type3.group <- type3.group
# i=1
# data = mlfa$data
# df = df
# factor = i
# time = mlfa$time
# strata = mlfa$type3.strata
# group = mlfa$type3.group
# arm_label <- "newTRT"
# placebo_label <- "Placebo"
# m_smooth = 2

# data$TRT01P <- relevel(C(factor(data$TRT01P), contr = "contr.treat"), ref = "placebo")


mlfa_type3 <- function(data,
                       df = 5,
                       factor = 1,
                       time = "AVISITN",
                       strata = "SEX",
                       group = "USUBJID",
                       arm_label = "newTRT",
                       placebo_label = "Placebo",
                       m_smooth = NA) {
  data[, placebo_label] <- 1
  form_no_inter <- as.formula(paste("~ -1 +", arm_label))
  arm_hot <- model.matrix(object = form_no_inter, 
                          model.frame(formula = form_no_inter, data = data))
  arm_hot <- arm_hot[, colnames(arm_hot) != paste0(arm_label, placebo_label)]
  colnames(arm_hot) <- gsub(" ", "", colnames(arm_hot))
  data[, colnames(arm_hot)] <- arm_hot
  formula <- paste0("z", as.character(factor), " ~ ", 
                    "s(", time, ", k = ", df, ", m = ", m_smooth, ")")
  for (nam_curr in colnames(arm_hot)) {
    formula <- paste0(formula, " + ", 
                      "s(", time, ", k = ", df, ", by = ", nam_curr, ", m = ", m_smooth,  ", pc = 0)")
  }
  formula
  
  formula <- stats::as.formula(formula)
  formula_varFixed <- stats::as.formula(paste0("~ ", "s", as.character(factor), "^2"))
  formula_corCar1 <- stats::as.formula(paste0("~ ", time, " | ", group))
  random_list <- list()
  random_list[[group]] <- ~1
  # gamm_object <- mgcv::gamm(formula = as.formula(formula),
  #                           random = list(USUBJID = ~ 1), # sym(group)
  #                           correlation = nlme::corCAR1(form = formula_corCar1),
  #                           weights=nlme::varFixed(formula_varFixed),
  #                           data = data )
  
  
  
  
  # time_fac <- factor(AVISITN)
  data$time_fac <- factor(data$AVISITN)
  t_unique <- unique(data$AVISITN)
  lme_out <- nlme::lme(fixed = z1 ~ -1 + time_fac + arm_hot:time_fac,
                       random = list(USUBJID = ~ 1), # sym(group)
                       # correlation = nlme::corCAR1(form = formula_corCar1),
                       weights = nlme::varFixed(~s1^2),
                       data = data)
  
  summ_curr <- summary(lme_out)
  mn_fixed <- summ_curr$coefficients$fixed
  cov_fixed <- summ_curr$varFix
  se_fixed <- sqrt(diag(cov_fixed))
  coef_plac <- paste0("time_fac", t_unique)
  coef_treat <- paste0(coef_plac, ":arm_hot", colnames(arm_hot))
  
  plot(mn_fixed[coef_plac], ty = "l")
  lines(mn_fixed[coef_plac] + 2 * se_fixed[coef_plac], lty = 2)
  lines(mn_fixed[coef_plac] - 2 * se_fixed[coef_plac], lty = 2)
  
  d_subject_sub$newTRT <- gsub(" ", "", d_subject_sub$newTRT)
  
  trt_unique <- unique(d_subject_sub$newTRT)
  par(mfrow = c(2, 2))  
  for (trt_curr in trt_unique) {
    if (trt_curr == placebo_label) {
      coef_names <- paste0("time_fac", t_unique)
    } else {
      coef_names <- paste0("time_fac", t_unique, ":arm_hotnewTRT", trt_curr)
    }
    plot(mn_fixed[coef_names], ty = "l")
    lines(mn_fixed[coef_names] + 2 * se_fixed[coef_names], lty = 2)
    lines(mn_fixed[coef_names] - 2 * se_fixed[coef_names], lty = 2)
  }
  
  
  
  plot(mn_fixed[nam_curr], ty = "l")
  lines(mn_fixed[coef_plac] + 2 * se_fixed[coef_plac], lty = 2)
  lines(mn_fixed[coef_plac] - 2 * se_fixed[coef_plac], lty = 2)
  
  
  
  plot(mn_fixed[coef_plac])
  coef_treat %in% names(mn_fixed)
  
  corFix <- summ_curr$corFixed
  corFix <- summ_curr$modelStruct$
    
    summary(lme_out)$coefficients
  coef_names <- paste0("time_fac", t_unique)
  plot(time_fac)
  str(coef(lme_out))
  # stored for plotting
  gamm_object$time <- time
  gamm_object$strata <- strata
  gamm_object$group <- group
  gamm_object$factor <- factor
  
  
  
  
  return(gamm_object)
}



#' @title
#' Fit a \emph{Gaussian Additive Mixture Model} (type 1).
#' 
#' @description
#' This helper function wraps a common longitudinal model, especially for biomedical data, described in details.
#' 
#' @details
#' The underlying \code{mgcv} call is as follows:
#' 
#' \code{
#' gamm(  
#'   formula = z ~ s(time, k = df) + strata[1] + s(time, k = df, by=factor(strata[1], ordered=TRUE)) + strata[2] + ... ,
#'   random = list(group=~1),
#'   data = data
#' )
#' }
#'
#' Possible variables are: 
#' \itemize{
#'    \item \code{time}: the timesteps, e.g. in weeks
#'    \item \code{strata}: the trial variable
#'    \item \code{group}: the subject variable
#' }
#'
#' @param data tibble; A tibble data frame where each row contains an observation
#'                        of multiple variables from a certain series (e.g.
#'                        a subject) and at a certain timestep (e.g. week). The columns
#'                        contain the series, the timestep and arbitrarily many
#'                        other variables.
#' @param df integer; The degrees of freedom of the spline bases.
#' @param factor integer; The factor (column of \code{Z} and \code{S}) chosen that shall be regressed on.
#' @param time string: The column name in \code{data} of the longitudinal variable chosen.
#' @param strata character vector: The column names in \code{data} that specify the factor for which an intercept and spline is separately estimated.
#' @param group string: The column name in \code{data} that specify the factor for which a random intercept is estimated.
#' @return gamm_object \code{gamm} object; a gamm object as reutrn by the \code{gamm(...)} call.
mlfa_type1 <- function(data, df, factor = 1, time, strata, group) {
  formula <- paste0(paste0("z", as.character(factor), " ~ ", "s(", time, ", k = ", as.character(df), ")"))
  for (i in 1:length(strata)) {
    formula <- paste0(formula, " + ", strata[i], " + ", "s(", time, ", k = ", as.character(df), ", by = factor(", strata[i], ", ordered = TRUE))")
  }
  formula <- stats::as.formula(formula)
  
  random_list <- list()
  random_list[[group]] <- ~1
  
  gamm_object <- mgcv::gamm(   # !!(sym(paste("z", as.character(PC_sel), sep = "")))
    formula = formula,
    random = random_list, # sym(group)
    data = data
  )
  
  # stored for plotting
  gamm_object$time <- time
  gamm_object$strata <- strata
  gamm_object$group <- group
  gamm_object$factor <- factor
  
  return(gamm_object)
}




#' @title
#' Plot (and predict) the latent trajectories of stage 3.
#' 
#' @description
#' This function first predicts and then plots the longitudinal trajectories per factor.
#'
#' @param x mlfa S3 object
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par). See \code{?plot} for further details.
#' @param save_dir string; the directory to save the plots in.
plot.mlfa <- function(x, ..., save_dir = NULL) {
  plot_list <- list()
  for (i in 1:x$fac_dim) {
    print(paste0("Plotting factor ", as.character(i), " ..."))
    gamm_object <- x$gamm_objects[[i]]
    data_predict <- predict_traj(gamm_object = gamm_object, data = x$data)
    
    if (x$type == 1) {
      p <- plot_traj_type1(gamm_object = gamm_object, data = data_predict)
    }
    if (x$type == 2) {
      p <- plot_traj_type2(gamm_object = gamm_object, data = data_predict)
    }
    plot_list[[i]] <- p
  }
  print("Plotting done.")
  return(plot_list)
}





#' @title 
#' Predict the longitudinal trajectories.
#' 
#' @description
#' Predict the longitudinal trajectory of a \code{gamm} model (single factor), using the fitted fixed-effects only.
#'
#' @param gamm_object gamm; A \code{gamm} object, returned as an attribute of the mlfa function call.
#' @param data tibble; The data frame used as input to \code{mlfa} function call.
#' @param conf_SEs int; the number of Standard Errors (SEs) displayed as confidence bands around mean. By default 2 SEs.
#' @param step int; The step size over which a prediction should be made for the \code{time} variable.
#' @return tibble; A tibble that contains a grid expansion of the fixed component variables (part of formula in \code{gamm} call),
#'                 and their corresponding prediction (fit (\code{fit} attribute), upper (\code{ub} attribute) and lower (\code{lb} attribute) confidence bonds).
#'                 Note that only the fixed variables, not the variables of the random effects, are used for prediction.
predict_traj <- function(gamm_object, data, conf_SEs = 2, step = 0.1) {
  expand_list <- list()
  # special treatment for "time" variable to make "smooth" prediction.
  expand_list[[gamm_object$time]] <- seq(data %>% dplyr::select(gamm_object$time) %>% min(), data %>% dplyr::select(gamm_object$time) %>% max(), by=step)
  # all.vars(m1$gam$formula) returns all variables in the formula
  # exclude the first one which is 'z'
  for (formula_var in all.vars(gamm_object$gam$formula)[2:length(all.vars(gamm_object$gam$formula))]) {
    # time was added above already
    if (!(formula_var == gamm_object$time)) {
      expand_list[[formula_var]] <- data %>% dplyr::select(formula_var) %>% unique() %>% dplyr::pull()
    }
  }
  data_predict <- expand.grid(expand_list)
  
  predict_output <- stats::predict(gamm_object$gam, newdata = data_predict, se.fit = TRUE)
  # predicted mean
  data_predict$fit <- predict_output$fit
  # lower confidence band
  data_predict$lb <- predict_output$fit - conf_SEs*predict_output$se.fit
  # upper confidence band
  data_predict$ub <- predict_output$fit + conf_SEs*predict_output$se.fit
  
  data_predict <- dplyr::as_tibble(data_predict)
  
  return(data_predict)
}


#' @title 
#' Latent trajectory (type 1)
#' 
#' @description
#' Plot the latent trajectory over time (type 1).
#'
#' @param gamm_object gamm object; The \code{gamm} object used to predict the data.
#' @param data tibble; A tibble containing the predicted data, as outputted by \code{predict_traj(...)}.
#' @param fill boolean; whether to fill the standard error bars with color or not.
#' @return ggplot object.
plot_traj_type1 <- function(gamm_object, data, fill = TRUE) {
  # make a "joint" strata factor for plotting
  group <- gamm_object$strata[1]
  if (length(gamm_object$strata) > 1) {
    for (i in 2:length(gamm_object$strata)) {
      group <- paste0(group, "_", gamm_object$strata[i])
    }
  }
  # adjust the data for plotting with joint factor
  joint_factor <- data %>% select(gamm_object$strata[1]) %>% pull()
  if (length(gamm_object$strata) > 1) {
    for (i in 2:length(gamm_object$strata)) {
      joint_factor <- paste0(joint_factor, "_", data %>% select(gamm_object$strata[i]) %>% pull())
    }
  }
  new_column <- as_tibble(factor(joint_factor))
  colnames(new_column) <- group
  data <- dplyr::bind_cols(data, new_column)
  
  p <- ggplot2::ggplot(mapping = ggplot2::aes_string(x=gamm_object$time, y = "fit", group = group, fill = group), data=data)
  
  factor <- all.vars(gamm_object$gam$formula)[1]
  factor <- factor[2:length(factor)]  # just the number
  if (fill) {
    p<- p+  scale_color_brewer(palette="Paired")+
      geom_line(size = 2,mapping = ggplot2::aes_string(y = "fit", colour = group), data = data) +
      geom_ribbon(mapping = ggplot2::aes_string(ymin = "lb", ymax = "ub", colour = group, fill = group),data = data, alpha = 0.3, linetype=0)
  }
  else {
    p<- p+  scale_color_brewer(palette="Paired")+
      geom_line(size = 2,mapping = ggplot2::aes_string(y = "fit", colour = group), data = data) +
      geom_ribbon(mapping = ggplot2::aes_string(ymin = "lb", ymax = "ub", colour = group, fill = NA),data = data, alpha = 0.3, linetype= "dashed")
  }
  
  p<-p+xlab("Time")+ylab("Z Scores")+
    theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) +
    theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) +
    theme(legend.position = "bottom", legend.text=element_text(size=18)) +
    theme(legend.key.width = unit(2, "cm")) + 
    theme(strip.text.y = element_text(size = 18), strip.text.x = element_text(size = 18))+ 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ggtitle(paste0("Latent Factor ", as.character(gamm_object$factor), " over ", gamm_object$time))+   
    theme_bw() 
  
  
  ############
  methods::show(p)
  return(p)
}





#' @title
#' Latent trajectory (type 2)
#' 
#' Plot the latent trajectory over time (type 2).
#'
#' @param gamm_object gamm object; The \code{gamm} object used to predict the data.
#' @param data tibble; A tibble containing the predicted data, as outputted by \code{predict_traj(...)}.
#' @param save_path string; The path where to save the plotted file.
#' @return ggplot object.
plot_traj_type2 <- function(gamm_object, data, save_path = NULL) {
  # exactly same plot as for type 1
  p <- plot_traj_type1(gamm_object, data, save_path)
  return(p)
}


#' @title
#' Summary of mlfa object.
#' 
#' @description
#' Print a summary of an mlfa object by concatenating the summaries of the \code{gam} and \code{lme} objects in its attributes.
#' 
#' @param object mlfa; An \code{mlfa} object.
#' @param ... additional arguments affecting the summary produced.
summary.mlfa <- function(object, ...) {
  summaries <- list()
  for (i in 1:length(object$gamm_objects)) {
    gamm_object <- object$gamm_objects[[i]]
    summaries[[2*i-1]] <- summary(gamm_object$gam)
    summaries[[2*i]] <- summary(gamm_object$lme)
    print(paste0("Summaries of gamm object ", as.character(i), ":"))
    print(paste0("gam object (model ", as.character(i), "):"))
    print(summaries[[2*i-1]])
    print(paste0("lme object (model ", as.character(i), "):"))
    print(summaries[[2*i]])
  }
  return(summaries)
}









# --------------------------------------------------
# --------------------------------------------------
# --------------------------------------------------




#
# # fit.longitudinal.model(Z1 ~ age + sex + trt01p + spline(AVISITN) # G.'s model #, random=~USUBJID, cor=AR1(), weights = object$S, data = longitudinal.modeldata)
# # gmm <- gamm(z ~ treat + s(vis, k=10) + s(vis, k=10, by=factor(treat, ordered=TRUE)), random = list(subj=~1), correlation=corAR1(form=~1|subj), weights=1/z.se, data=dat)
#
# # parameters
# z_fac <- !!sym("z_fac_1")  # interpret as symbol for input in formulas
# spline_df <- 50
#
# # radial basis function
# # TODO can this be used? https://cran.r-project.org/web/packages/kergp/kergp.pdf
# # covRadial ?
# # corGaus?
#
# formula <- as.formula(cat(z_fac), " ~ ", "TODO")
#
# mlfa <- gamm(formula = z_fac ~
#             subject_ID   # subject-specific intercept
#             + s(visit_ID, bs = "tp")  # bs specifies type of smoothing, here: "tp": thin plate regression spline
#             + s(subject_ID, visit_ID, bs="re")
#              correlation = corGaus(form = ~ series_step), random = list(u=~1|subj_ID), data = mfd$data_out)
#
# s(labVar1, time)
#
#
# gmm <- gamm(z ~ treat + s(vis, k=10) + s(vis, k=10, by=factor(treat, ordered=TRUE)), random = list(subj=~1), correlation=corAR1(form=~vis | subj), weights=1/z.se, data=dat)
#
# Questions:
# - how to assess model quality?
# - nlme
#
#
# # TODO questions to Matthias:
# # - how to implement radial basis function kernel?
# -
#
#
# # TODO separate s for every
#
# # epsilon?
#
# # TODO subset by subject ???
# # TODO using groupedData package?
# # TODO how to use spline coefficients
#
#
#
#
