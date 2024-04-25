# make non-random
# set.seed(1)

# dependencies =============================================================
# TODO is this required here?
library(tidyverse)
library(pcaMethods)
library(itertools)  # for zip(...) like for looping
# previously required dependencies (kept for now) -------
# library(gsubfn)  # need 0.7-0 or later; for assigning multiple values to list
# library(MASS)  # ginv is imported from MASS (see NAMESPACE file and call with MASS:: below)
# library(rapportools) # for is.boolean(...) function
library("gplots")


#' @title
#' Multivariate Factor Decomposition.
#'
#' @description
#' Estimate factor loadings and socres (stage 1 and 2) in preparation of longitudinal analysis.
#' 
#' @details
#' This function performs stage 1 and 2 (and preparatory steps for stage 3) of the proposed method. This corresponds to the following functionality:
#' Preprocessing: 
#' \itemize{
#'   \item normalization of the input data (per measurement standardize -> all measurements on same scale)
#'   \item removing rows with NAs only
#' }
#' Further useful functions for preprocessing in this package are: \code{as_ranks(...)}, \code{sampling(...)}.
#' Stage 1: 
#' \itemize{
#'   \item Estimating a loadings matrix \code{A} via Probabilistic Principal Component Analysis (PPCA) -> see \code{estimate_loadings(...)} function for details
#'   \item Inducing sparsity into A by varimax/promax rotation
#'   \item Inducing sparsity into A by a threshold based method.
#' }
#' Stage 2:
#' \itemize{
#'   \item Estimating the scores matrix \code{Z} and standard error \code{S} with linear regression,
#'         corresponding to the transformed loadings matrix \code{A} as outputted in stage 1 -> see \code{estimate_scores(...)} function for details
#' }
#' This function returns an \code{mfd} object that has various attributes. The most important result attributes are:
#' \itemize{
#'   \item \code{mfd$A} -> the loadings matrix, dimensions: (PxK)
#'   \item \code{mfd$Z} -> the scores matrix, dimensions: (NxK)
#'   \item \code{mfd$S} -> the standard error matrix, dimensions: (NxK)
#'   \item \code{mfd$data_out} -> a tibble that contains
#'       \itemize{
#'           \item the columns of \code{mfd$Z}
#'           \item the columns of \code{mfd$S}
#'           \item the columns \code{time} and \code{subject} in \code{data}
#'      }
#'   \item \code{mfd$ppca_object} -> a ppca object as fitted in stage 1.
#' }
#'
#' Thus, given a data matrix \code{Y} with dimensions (NxP), which are \code{vars} columns in \code{data} after the steps detailed in Preprocessing above, the output of mfd then approximately satisfies: \code{Y = mfd$Z \%*\% t(mfd$A)}
#'
#' Furthermore, there are several helper functions for an mfd object created with this function that allow the following functionality:
#' \itemize{
#' \item \code{plot.mfd(mfd_object)}: this helper function plots
#'   \itemize{
#'      \item the cumulative R^2 score over the numer of principal components (PCs) with decreasing eigenvalue chosen
#'      \item a heatmap of the loadings matrix \code{mfd$A}
#'      \item the trajectories of the columns of \code{mfd$Z} over \code{time}, grouped by \code{subject}
#'   }
#' \item \code{summary.mfd(mfd_object)}: this helper function prints further diagnostic information (\code{summary} call) of the \code{ppca} object in stage 1.
#' }
#'
#' Differences to G.'s code:
#' \itemize{
#'     \item no second sparsity removal of variables
#' }
#'
#' Further extensions: 
#' \itemize{
#'     \item interactive mode where n_PCs is chosen based on user input when seeing the variance explained over number of PCs plot
#' }
#'
#' @param data tibble; (preprocessed) tibble to perform stages 1 and 2 on.
#'                        The tibble expects the columns \code{vars}, "time", and "subject" to be present.
#'                        Each row represents the observations from a certain subject at a specific timestep.
#'                        \code{data} requires all values in the \code{vars} columns to be numeric.
#'                        You can convert columns to a numeric representation with the helper function \code{convert_to_numeric(...)}.
#'                        \code{convert_to_numeric(...)} currently supports columns of type binary, unordered categorical, unordered categorical with baseline, 
#'                        and ordered categorical.
#' @param vars string vector; variables (columns) to select of data for operations.
#' @param time string; variable (column) in data that represents the longitudinal dimension. The sequence identifier is passed in parameter subject.
#'                     Note that this variable is only stored in the mfd object and enforces that the column is present, but is not directly used during stage 1 and 2.
#' @param subject string; variable (column) in data that represents a subject/ID/sequence that has multiple entries in column time
#'                        Note that this variable is only stored in the mfd object and enforces that the column is present, but is not directly used during stage 1 and 2.
#' @param loadings.n_PCs int; number of Principal Components selected (if loadings.n_PCs=="max" select the maximum possible number of PCs, usually len(vars)-1).
#' @param loadings.sparse_abs_thres positive double in [0.0; 1.0]; sparsity inducing threshold, all values in \code{A} with a value smaller than \code{sparse_thres} are set to 0  (\code{A} was normalized to have values in [0;1] before).
#' @param scores.type how to compute the scores posterior: 'lin_reg' which does a simple least squares to estimate the scores.
#'                                                         'a_posteriori' which computes a posterior mean and variance.
#' @param normalize.type string; Normalize data for loadings and scores estimation before calculation of PPCA. Options include \code{"per_group_center_per_meas_standard"}, \code{"per_meas_standard"}, \code{"per_scale"}, \code{"per_largest_abs"}, \code{"per_center_baseline_scale"}, \code{"per_center_meas"}, OR \code{"none"}.
#'             \code{"per_group_center_per_meas_standard"}: removes subject specific variability and centers the entire measurement;
#'             \code{"per_meas_standard"}: normalize by subtracting the mean and dividing by SD;
#'             \code{"per_center_meas"}: normalize by subtracting the column mean(s);
#'             \code{"per_scale"}: normalize with dividing by SD;
#'             \code{"per_largest_abs"}: normalize with dividing by largest absolute value column wise;
#'             \code{"per_center_baseline_scale"}: normalize with subtracting the baseline column mean(s) and then dividing by SD;
#'             \code{"none"}: without any normalization.
#' @param standardize_A_type string; standardization for \code{A} after PPCA step for \code{loadings.sparse_abs_thres}=0 or >0 and after inducing sparsity if \code{loadings.sparse_abs_thres}>0. 
#'                           For options of standardization type, \code{type = "largest_abs"} is to standardize by dividing by largest absolute value of column;  \code{type = "orthonormal"} is to obtain the orthonormal/singular vector decomposition matrix 
#'                        so that colSums(mat^2)=I; \code{type = "L2"} is to divide by the square root of column wise summation after square so that colSums(mat^2)=I; \code{type = "none"} is to keep original matrix without any standardization. The default standardization approach is "L2".  
#' @param Time_base string or numeric; if Value in data that represents the baseline in the longitudinal dimension. 
#' @param seed \code{NULL} or int; make the function random (\code{NULL}) or deterministic (int)
#' @param loadings.rotation \code{"varimax"} or \code{"promax"}; the type of sparcity inducing rotation applied to the loadings matrix after PPCA
#' @return \code{mfd} S3 object; an object of an S3 class \code{mfd}, with attributes as discussed above and accessible via \code{$}.


# 
# data = sub_data$time
# vars = var_names
# loadings.n_PCs = n_fac
# seed = 1
# time = "Time"
# subject = "ID"
# normalize.type = c("per_meas_standard", "standardize", "per_center_baseline_scale")[1]
# Time_base = 0
# standardize_A_type = c("none", "L2", "largest_abs")[3]
# loadings.sparse_abs_thres = 0
# loadings.rotation = c("varimax", "promax")[1]
# scores.type = 'a_posteriori'
# trace = FALSE

mfd <- function(# primary (necessary) inputs
          data,
          vars,
          time,
          subject,
          # secondary inputs
          loadings.n_PCs = "max",
          loadings.sparse_abs_thres = 0,
          scores.type='a_posteriori',
          normalize.type="per_center_baseline_scale",
          Time_base=0,
          standardize_A_type = "L2",
          # tertiary inputs
          seed = 1,
          loadings.rotation = "promax",
          trace = TRUE) {
          
  # checks ---
  if (!tibble::is_tibble(data)) {
    stop("data is not a tibble.")
  }
  # at least two variables must be present in group to reasonably estimate principal component of them
  if (ncol(data) < 2) {
    stop("Less than two columns present in data -> cannot reasonable estimate PCs.")
  }
  if (!all(vars %in% colnames(data))) {
    stop("vars must be columns of data.")
  }
  col_types <- sapply(data %>% select(vars), class)
  for (i in length(vars)) {
    col <- vars[i]
    col_type <- col_types[i]
    if (!"numeric" %in% col_type) {
      stop(paste0(col, " is not numeric."))
    }
  }
  if (is.null(time) || !time %in% colnames(data)) {
    stop("time must be column of data.")
  }
  if (is.null(subject) || !subject %in% colnames(data)) {
    stop("subject must be column of data.")
  }
  if ((!(loadings.n_PCs == "max") && (!(loadings.n_PCs == as.integer(loadings.n_PCs)) || !(loadings.n_PCs >= 2) || !(loadings.n_PCs < length(vars)))) ) {
    stop("loadings.n_PCs must be either 'max' or a positive integer >= 2 and < length of vars.")
  }
  if (!(is.null(seed)) && (!(seed == as.integer(seed)) || (!seed >= 1))) {
    stop("seed must be NULL or positive interger.")
  }
  if (!loadings.rotation %in% c("promax", "varimax")) {
    stop("loadings.rotation must be either 'varimax' or 'promax'.")
  }
  
  
  # (checks end) ---
  
  # create object and object attributes ===
  # ===
  # create object of class multivariate longitudinal factor analysis
  # this object will carry the results of all operations and will also be returned from this function
  mfd <- structure(list(), class = "mfd")
  
  # initialize all attributes as NULL (filled later, if existent) ---
  # documenation on attributes with S3: http://adv-r.had.co.nz/Data-structures.html#attributes
  attr(mfd, "data_out") <- NULL
  
  attr(mfd, "A") <- NULL
  attr(mfd, "Z") <- NULL
  attr(mfd, "S") <- NULL
  attr(mfd, "fac_dim") <- NULL
  attr(mfd, "ppca_object") <- NULL
  
  attr(mfd, "time") <- NULL
  mfd$time <- time
  attr(mfd, "subject") <- NULL
  mfd$subject <- subject
  
  # main procedural logic of multivariate longitudinal factor analysis ===
  # ===
  # set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # assign data as variable within function
  data <- data
  
  # normalize data for loadings and scores estimation ---
  data_scores<-data_loadings <- normalize(data = data, 
                                          columns = vars,
                                          subject=subject, 
                                          type = normalize.type,
                                          Time = time,
                                          Time_base = Time_base)
  
  # select columns and transform to matrix ---
  Y_scores<-Y_loadings <- select_columns(data = data_loadings, columns = vars)
  Y_scores<-Y_loadings <- as.matrix(Y_loadings)
  
  mfd$mean_location_shift <- colMeans(data[, vars]) - colMeans(Y_loadings[, vars])
  mfd$sd_scale_factor <- apply(data[, vars], 2, sd, na.rm = T) / apply(Y_loadings[, vars], 2, sd, na.rm = T)
  
  
  # remove rows that contain NAs only in Y_loadings/Y_scores, also remove in data_loadings/data_scores
  list[Y_loadings, Y_scores, data_loadings, data_scores] <- 
    remove_all_NA_rows(Y_loadings, Y_scores, data_loadings, data_scores)
  if (trace) {
    print("Preprocessing done.")
  }
  # TODO potentially change estimate_loadings, estimate_scores, too -> one instead of two parameters
  if (trace) {
    print("Estimate loadings (stage 1) started...")
  }
  list[A, ppca_object] <- estimate_loadings(Y_in = Y_loadings, 
                                            type = "n_PCs", 
                                            var_explained = NULL,
                                            n_PCs = loadings.n_PCs,
                                            rotation = loadings.rotation, 
                                            conv_thres = 1e-4,
                                            ppca_seed = if (is.null(seed)) NA else seed)
  A_unnormalized <- A
  A <- normalize_loadings(mat = A,
                          type = standardize_A_type )
  
  if (abs(loadings.sparse_abs_thres)>0){ 
    A <- sparse_heur_thres(A, thres = loadings.sparse_abs_thres)
    A <- normalize_loadings(mat = A,type =standardize_A_type)
  }  
  
  
  if (trace) {
    print("Estimate loadings (stage 1) done.")
    print("Estimate scores (stage 2) started...")
  }
  estimate_scores_out <- estimate_scores(Y_in = Y_scores, A_in = A, type = scores.type)
  list[Z, S] <- estimate_scores_out[c("Z", "S")]  # type is "a_posteriori" OR "lin_reg"
  sigma_square_ml <- estimate_scores_out$sigma_square_ml
  if (trace) {
    print("Estimate scores (stage 2) done.")
  }  
  # renaming of rows and columns
  rownames(A) <- vars
  colnames(A) <- paste("a", seq(1:dim(A)[2]), sep = "")
  # each column of Z and S corresponds to one principal component (Z and S dimensions: (NxK))
  # name columns according to which factor (principal component) they belong to --> will be columns in resulting tibble
  colnames(Z) <- paste("z", seq(1:dim(Z)[2]), sep = "")
  colnames(S) <- paste("s", seq(1:dim(Z)[2]), sep = "")
  
  mfd$A <- A
  mfd$A_unnormalized <- A_unnormalized
  mfd$Z <- Z
  mfd$S <- S
  mfd$sigma_square_ml <- sigma_square_ml
  mfd$ppca_object <- ppca_object
  
  mfd$fac_dim <- dim(A)[2]
  
  # Following steps prepare the input data for stage 3, given outputs of stage 1 and 2.
  # The preparation includes the followings steps:
  # #
  if (sum(is.infinite(A)) > 0) {
    stop("A contains inf/-inf (or NA) values.")
  }
  
  # convert matrices to tibble for merging
  Z <- dplyr::as_tibble(Z)
  S <- dplyr::as_tibble(S)
  
  # merge Z and S with original data
  # data_scores has the normalization used for stage 3
  data_out <- data_scores %>% dplyr::select(time, subject)  # only use time and subject column
  data_out <- data_out %>% dplyr::bind_cols(list(Z, S))
  
  # assign data_out to mfd object
  mfd$data_out <- data_out
  
  return(mfd)
}





#' @title
#' Plot various diagnostics of stage 1 and 2.
#' 
#' @description
#' This wrapper function calls several plot functions that produce diagnostic and result plots of stage 1 and 2,
#'  given the corresponding \code{mfd} object.
#' 
#' @details
#' Three types of plots are computed:
#' \itemize{
#'   \item Cummulated R^2 over chosen principal components
#'   \item Loadings heatmap
#'   \item Longitudinal \code{z} trajectories (one for each factor).
#' }
#'
#' Called through generic \code{plot(...)} function.
#' Documentation on generic functions: http://adv-r.had.co.nz/S3.html#undefined
#'
#' @param x mfd; the result object of the mfd function.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par). See \code{?plot} for further details.
#' @param save_dir string; the directory where to save the plots in.
#' @param colour_by string = NULL; see \code{plot_z_trajectory(...)} for details.
#' @param subject_data tibble = NULL; A tibble with a column \code{mfd$subject} and various other columns that contain subject specific information.
#'                             This information can be used e.g. for styles (colouring) of various plots.
#'                             Default: \code{NULL} (not merging any data).
plot.mfd <- function(x, ..., save_dir = NULL, colour_by = NULL, subject_data = NULL) {
  plot_mfd_ppca(x, save_path = if (!is.null(save_dir)) file.path(save_dir, "CumRsquare_over_PCs.pdf") else NULL)
  plot_loadings_heatmap(x, save_path = if (!is.null(save_dir)) file.path(save_dir, "loadings_heatmap.pdf") else NULL)
  for (i in 1:x$fac_dim) {
    plot_z_trajectory(x, save_path = if (!is.null(save_dir)) file.path(save_dir, paste0("z", as.character(i), "_over_", as.character(x$time), ".pdf")) else NULL, factor = i, time = x$time, subject = x$subject,
                      colour_by = colour_by, subject_data = subject_data, selected_subjects = NULL)
  }
}



#' @title
#' Plot cummulative R^2 of \code{ppca} fit.
#' 
#' @description
#' Plot cummulative R^2 of \code{ppca} fit over (number) of principal components (stage 1 result).
#'
#' @param mfd_object mfd; the result object of \code{mfd}.
#' @param font.size number; the font size displayed in the plot.
#' @param save_path string; the path where to save the file to. Must end with \code{".pdf"}.
#' @param label.size number; the font size of the label above each bar.
#' @param label.vjust number; the vertical height adjustment for label displayed in each bar.
#' @param ylab string; the y axis name displayed.
#' @param labels vector; the label name for each lading.
#' @return ggplot object.
#' @color bar plot color
plot_mfd_ppca <- function(mfd_object,color="#9999CC",font.size=18, save_path = NULL,label.size=7,label.vjust=-0.25, ylab="Cumulated R^2",labels=c("Domain1","Domain2","Domain3")) {
  data <- dplyr::as_tibble(data.frame(PCs = base::colnames(mfd_object$A), R2cum = mfd_object$ppca_object@R2cum))
  
  if (!is.null(labels)){
    data$PCs <- factor(data$PCs, levels=paste0("a", seq(1:dim(mfd_object$A)[2])),labels = labels)  # levels for sorting
  } else{
    data$PCs <- factor(data$PCs, levels=paste0("a", seq(1:dim(mfd_object$A)[2])))  # levels for sorting
  }
  
  p <- ggplot2::ggplot(data=data, ggplot2::aes(x=PCs, y=R2cum)) +
    ggplot2::geom_bar(stat="identity", fill = color) + ylim(c(0,1))+
    ggplot2::geom_text(aes(label=round(R2cum,2)), position=position_dodge(width=0.9), vjust=-label.vjust,size=label.size)+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.y = ggplot2::element_text(11)) 
  
  p <- p + ggplot2::ggtitle("Variance explained by Principal Components") + 
    ggplot2::xlab("Principal Components (sorted)") + ggplot2::ylab(ylab)+
    theme_bw() + 
    theme(axis.text.x = element_text(size=font.size), axis.text.y = element_text(size=font.size)) +
    theme(axis.title.x = element_text(size=font.size), axis.title.y = element_text(size=font.size)) +
    theme(legend.position = "bottom", legend.text=element_text(size=font.size)) +
    theme(legend.key.width = unit(2, "cm")) + 
    theme(strip.text.y = element_text(size = font.size), strip.text.x = element_text(size = font.size))
  
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, dpi=2000, device = cairo_pdf)
  }
  # methods::show(p)
  return(p)
}



#' @title
#' Plot heatmap of loading matrix \code{A} or latent trajectory weight matrix from \code{calculate_latent_weights}.
#' 
#' @description
#' Plot a heatmap of the coefficient matrix \code{A} (variables by (sparse) principal components) or latent trajectory weight matrix transformed from matrix \code{A}.
#'
#' @param input_coefficients matrix; the result \code{A} matrix from \code{mfd} or the matrix from \code{calculate_latent_weights}.
#' @param save_path string; the path where to save the file to. Must end with \code{".pdf"}.
#' @param axis_y_text_size the size of the y tick labels (measurements).
#' @param axis_x_text_size the size of the x tick labels (loadings).
#' @param title_text_size the size of the title text.
#' @param color.pallette the color pallette that display colors for scores from low to high
#' @param num_color_bins the number of color bins that display in the heatmap.
#' @param margins the margins of the heatmap plot 
#' @param title the label display for title
#' @param border_color the color for cell border in heatmap

plot_loadings_heatmap <- function(input_coefficients, save_path, axis_y_text_size = NULL,axis_x_text_size=NULL,title_text_size=NULL,
                                  color.pallette=c("darkred","white","darkgreen"),num_color_bins = NULL,margins=NULL,
                                  title="Loadings Heatmap",border_color="grey") {
  if (is.null(axis_y_text_size)) {
    n_vars <- length(rownames(input_coefficients))
    if (n_vars <= 20) {
      axis_y_text_size <- 0.8
    }
    else {
      axis_y_text_size <- 0.7
    }
  }
  if (is.null(axis_x_text_size)) {
    n_loadings <- length(colnames(input_coefficients))
    if (n_loadings <= 5) {
      axis_x_text_size <- 0.8
    }
    else {
      axis_x_text_size <- 0.7
    }
  }
  if (is.null(margins)){
    margins=c(6,19)
  }
  if (is.null(num_color_bins)){
    num_color_bins = 25
  }
  if (is.null(title_text_size)){
    title_text_size = 0.8
  }
  my_palette <- colorRampPalette(color.pallette)(n = num_color_bins)
  
  par(cex.main=title_text_size) 
  p<-heatmap.2(input_coefficients, dendrogram = "none", Rowv = FALSE, Colv = FALSE,
               density.info="histogram",
               trace="none",
               col = my_palette,         
               main=title,
               cexRow=axis_y_text_size,cexCol=axis_x_text_size,margins=margins,
               # sepwidth=c(0.01,0.05),
               sepwidth=c(0,0),
               sepcolor=border_color,
               colsep=1:ncol(input_coefficients),
               rowsep=1:nrow(input_coefficients))
  
  if (!is.null(save_path)) {
    pdf(file.path(save_path,"heatmap.pdf" ))
    par(cex.main=title_text_size)
    p
    dev.off()
  }
  # methods::show(p)
}





#' @title
#' Summary of mfd object.
#' 
#' @description
#' Print a summary of an mfd object which prints its PPCA object's summary.
#' 
#' @param object mfd; the result object of the \code{mfd} function.
#' @param ... additional arguments affecting the summary produced.
summary.mfd <- function(object, ...) {
  # calling the method of pcaRes
  # Documentation: https://www.bioconductor.org/packages/devel/bioc/manuals/pcaMethods/man/pcaMethods.pdf
  print("Note: This summary shows diagnostic information of both the selected and not-selected PCs.")
  summary(object$ppca_object)
}


# === normalizing helper ===

#' @title
#' Normalize the data (helper)
#' 
#' @description
#' Normalize the variables in the data according to a certain strategy (helper function).
#'
#' @param data tibble; The data to transform.
#' @param columns vector of strings; the column names of variables to normalize.
#' @param type \code{"standardize"}, \code{"center"}, \code{"scale"}, \code{"largest_abs"}, \code{"scale_per_baseline"}, OR \code{"none"}; the type of normalization to apply(per grouping).
#'             \code{"standardize"}: normalize by subtracting the mean and dividing by SD;
#'             \code{"center"}: normalize by subtracting the column mean(s);
#'             \code{"scale"}: normalize with dividing by SD;
#'             \code{"largest_abs"}: normalize with dividing by largest absolute value column wise;
#'             \code{"scale_per_baseline"}: normalize with subtracting the baseline column mean(s) and then dividing by SD;
#'             \code{"scale_per_baseline_per_individual"}: normalize with subtracting the baseline score for each individual and then dividing by SD for each measurement;
#'             \code{"none"}: without any normalization.
#' @param group list of strings; the variable (column) by which to group for normalization (computes normalization statistics per variable and per group of rows).
#'                               If NULL: group by rows (default).
#' @param Time string; variable (column) in data that represents the longitudinal dimension. 
#' @param Time_base string or numeric; Value in data that represents the baseline in the longitudinal dimension. 
#' @return data tibble; The transformed data.
#'
#' @details
#' Possible extension:
#' \itemize{
#'   \item implement case "both" with median absolute distance
#' }
#'
#' Known issues:
#' \itemize{
#'   \item if sd(group measurements) = 0 -> NAs are produced during normalization
#' }
normalize_grouped <- function(data, columns,subject, type = "standardize", group = NULL,Time="AVISITN",Time_base=0) {
  # checks ---
  # check: data is tibble
  if (! class(data)[1] == "tbl_df") {
    stop("data is not a tibble.")
  }
  # check: group variable must exist
  if (!(is.null(group)) && !all((group %in% colnames(data)))) {  # all columns must be in data (logical and on vector)
    stop(paste0("Group variable(s) ", group, " not in colnames of data."))
  }
  # TODO check type of columns
  # check: columns must be in data
  if (!all(columns %in% colnames(data))) {  # all columns must be in data (logical and on vector)
    stop("Variable columns not in columns of data.")
  }
  # (checks end) ---
  
  # variables to group by for scaling
  if (!is.null(group)) {
    data <- data %>% dplyr::group_by(!!(dplyr::sym(group)))
  }
  
  if (type == "standardize") {
    data <- dplyr::mutate_at(data, columns, ~ as.numeric(scale(., center = TRUE, scale = TRUE)))   
  }
  else if (type == "center") {
    data <- dplyr::mutate_at(data, columns, ~ as.numeric(scale(., center = TRUE, scale = FALSE)))   
  }
  else if (type == "largest_abs") {
    data[,columns] <- t(t(data[,columns]) / apply(data[,columns], 2, function(v) max(abs(v), na.rm = T)))  
  }
  else if (type == "scale") {
    data <- dplyr::mutate_at(data, columns, ~ as.numeric(scale(., center = FALSE, scale = TRUE)))  
  }
  else if (type == "scale_per_baseline") {
    data$Time<-data[,Time]
    
    if (nrow(data %>% filter(Time==Time_base))==0) {
      stop("Need to specify the base.")
    }
    base_mean<-t(data %>% filter(Time==Time_base) %>% select(columns) %>% 
                   colMeans(.,na.rm = TRUE) %>% as.matrix)
    
    data_SD<-t(as.matrix(apply(data[,columns], 2, function(v) sd(v, na.rm = T))))
    
    for (i in columns){
      data_tmp<-(data[,i]-base_mean[,i])/data_SD[,i]
      data[,i]<-data_tmp}
  } 
  else if (type == "scale_per_baseline_per_individual") {
    data$Time<-data[,Time]
    
    if (nrow(data %>% filter(Time==Time_base))==0) {
      stop("Need to specify the base.")
    }
    
    data_base<-data %>% filter(Time==Time_base) %>% select(columns,subject) %>% unique()
    
    data_SD<-t(as.matrix(apply(data[,columns], 2, function(v) sd(v, na.rm = T))))
    
    for (i in columns){
      tmp_base<-data_base[,c(subject,i)] %>% dplyr::rename(base=i)
      
      data_tmp_v1<-data[,c(subject,Time,i)]%>% dplyr::rename(post=i) %>% full_join(tmp_base) %>% mutate(chg=post-base) %>% select(subject,Time,chg)
      colnames(data_tmp_v1) <- c(subject,Time,i)
      
      data_tmp_v2<-(data_tmp_v1[,i])/data_SD[,i]
      data[,i]<-data_tmp_v2
      
    }
  }
  
  else if (type == "none") {
    # do nothing
  } else {
    stop(paste0("Invalid variable type (== ", type, ") specified."))
  }
  
  # ungroup data
  data <- data %>% dplyr::ungroup()
  
  return(data)
}


# === normalizing ===


#' @title 
#' Normalize the data.
#' 
#' @description
#' Normalize the variables in the data according to a certain strategy (main function for specific strategies).
#'
#' @details
#' Possible Extensions:
#' \itemize{
#'   \item  third "weird" way of normalizing (see 'both' version in G.'s code)
#' }
#'
#' @param data tibble; The data to transform.
#' @param columns vector of strings; the column names of variables to normalize
#' @param type \code{"per_group_center_per_meas_standard"}, \code{"per_meas_standard"}, \code{"per_scale"}, \code{"per_largest_abs"}, \code{"per_center_baseline_scale"}, \code{"per_center_meas"}, OR \code{"none"}; the type of normalization to apply(per grouping).
#'             \code{"per_group_center_per_meas_standard"}: removes subject specific variability and centers the entire measurement, in G.'s code known as 'longit'
#'             \code{"per_meas_standard"}: normalize by subtracting the mean and dividing by SD;
#'             \code{"per_center_meas"}: normalize by subtracting the column mean(s);
#'             \code{"per_scale"}: normalize with dividing by SD;
#'             \code{"per_largest_abs"}: normalize with dividing by largest absolute value column wise;
#'             \code{"per_center_baseline_scale"}: normalize with subtracting the baseline score for each individual and then dividing by SD for each measurement;
#'             \code{"per_center_baseline_mean_scale"}: normalize with subtracting the baseline column mean(s) and then dividing by SD;
#'             \code{"none"}: without any normalization.
#' @param Time string; variable (column) in data that represents the longitudinal dimension. 
#' @param Time_base string or numeric; Value in data that represents the baseline in the longitudinal dimension. 
#' @param series_group string; if \code{type == "per_series_center_and_per_meas_center"}: The variable indicating the series to group by during normalization.
#' @return data tibble; The transformed data.
normalize <- function(data, columns,subject, type, series_group = NULL,Time="AVISITN",Time_base=0) {
  # checks ---
  if (! class(data)[1] == "tbl_df") {
    stop("data is not a tibble.")
  }
  # TODO check type of columns
  # check: columns must be in data
  if (!all(columns %in% colnames(data))) {  # all columns must be in data (logical and on vector)
    stop("Variable columns not in columns of data.")  # TODO !!! conduct this test in all functions where columns argument
  }
  # check: type correctly specified
  if (!type %in% c("per_group_center_per_meas_standard", 
                   "per_meas_standard",
                   "per_scale",
                   "per_largest_abs",
                   "per_center_baseline_scale",
                   "per_center_meas",
                   "none",
                   "standardize")) {
    stop("Variable type incorrectly specified.")
  }
  if ((type == "per_group_center_per_meas_standard") && (!series_group %in% colnames(data))) {
    stop("series_group column not present in data.")
  }
  # (checks end) ---
  if (type == "per_group_center_per_meas_standard") {
    data<-normalize_grouped(data = data, columns = columns, type = "center", group = series_group,Time=Time,Time_base=Time_base)
    data<-normalize_grouped(data = data, columns = columns, type = "standardize", group = NULL,Time=Time,Time_base=Time_base)
  }
  else if (type == "per_meas_standard") {
    data<-normalize_grouped(data = data, columns = columns, type = "standardize", group = NULL,Time=Time,Time_base=Time_base)
  }
  else if (type == "per_scale") {
    data<-normalize_grouped(data = data, columns = columns, type = "scale", group = NULL,Time=Time,Time_base=Time_base)
  }
  else if (type == "per_largest_abs") {
    data<-normalize_grouped(data = data, columns = columns, type = "largest_abs", group = NULL,Time=Time,Time_base=Time_base)
  }
  else if (type == "per_center_baseline_scale") {
    data<-normalize_grouped(data = data, columns = columns,subject=subject, type = "scale_per_baseline_per_individual", group = NULL,Time=Time,Time_base=Time_base)
  }
  else if (type == "per_center_baseline_mean_scale") {
    data<-normalize_grouped(data = data, columns = columns, type = "scale_per_baseline", group = NULL,Time=Time,Time_base=Time_base)
  }
  else if (type == "per_center_meas") {
    data<-normalize_grouped(data = data, columns = columns, type = "center", group = NULL,Time=Time,Time_base=Time_base)
  }
  else if (type == "none") {
    data<-normalize_grouped(data = data, columns = columns, type = "none", group = NULL,Time=Time,Time_base=Time_base)
  } 
  
  return(data)
}


# === select columns ===

#' @title
#' Select columns.
#' 
#' @description
#' Select specific columns of interest from the data.
#'
#' @param data tibble; The data to transform.
#' @param columns list of strings; The columns in \code{data} to select.
#' @return data tibble; The transformed data.
select_columns <- function(data, columns) {
  # check ---
  if (! class(data)[1] == "tbl_df") {
    stop("data is not a tibble.")
  }
  if (! all(columns %in% colnames(data))) {
    stop("One or more variables in columns not in data.")
  }
  # (checks end) ---
  
  data <- data[, columns, drop = F]  # TODO what does drop F do?
  
  return(data)
}

# === computing ranks of the data ===

#' @title
#' Transform variables to ranks.
#'
#' @description
#' Transform variables to ranks, i.e. each variable in columns is interpreted as ordinal.
#' Tie ranks will be averaged.
#'
#' @details
#' Potential extensions: 
#' \itemize{
#'   \item potentially add na.last, ties.method from internal function as input
#' }
#' 
#' @param data tibble; The data to be transformed.
#' @param columns vector of strings; the column names to transform to ranks.
#' @return data tibble; The transformed data.
as_ranks <- function(data, columns) {
  # checks ---
  if (! class(data)[1] == "tbl_df") {
    stop("data is not a tibble.")
  }
  if (! all(columns %in% colnames(data))) {
    stop("columns variables are not in data.")
  }
  # (checks end) ---
  
  # Function to compute the ranks of a measurement as decimal in interval [0.0; 1.0]
  #
  # TODO potential issue: what if min, max is missing
  # na.last: how to handle NAs -> "keep" assigns rank NA to NA values; ties.method: what is returned if two values are equal -> "average" returns average rank
  rank_dec_func <- function(x) {
    ranks <- rank(x, na.last = "keep", ties.method = "average")
    min_rank <- min(ranks, na.rm = T)
    max_rank <- max(x, na.rm = T)
    n_non_na_values <- length(ranks[!is.na(ranks)])
    dec_ranks <- (ranks - min_rank) / n_non_na_values  # na.last: how to handle NAs -> "keep" assigns rank NA to NA values; ties.method: what is returned if two values are equal -> "average" returns average rank
  }
  
  data <- dplyr::mutate_at(data, columns, rank_dec_func)
  
  return(data)
}


# === perform rotation ===

#' @title
#' Sparsity inducing rotation
#' 
#' @description
#' Performs a sparsity inducing rotation on matrix.
#' Helper function for \code{estimate_loadings(...)}.
#'
#' @param mat matrix; The matrix to rotate.
#' @param type \code{"promax"} or \code{"varimax"}: The type of rotation to perform.
#' @return mat matrix; The rotated matrix.
sparse_rotation <- function(mat, type) {
  # checks ---
  if (!is.matrix(mat)) {
    stop("mat is not a matrix.")
  }
  if (!(type %in% c("promax", "varimax"))) {
    stop("type is not a valid rotation type.")
  }
  # for rotation to be reasonable, mat should have at least two columns
  if (ncol(mat) < 2) {
    stop("mat has less than two columns so that rotation cannot be performed.")
  }
  # (checks end) ---
  rotmat <- switch(type,
                   promax = stats::promax(mat)$rotmat,
                   varimax = stats::varimax(mat)$rotmat)
  mat <- mat %*% rotmat  # perform the sparcity inducing rotation
  
  return(mat)
}


# === normalize A ===

#' @title
#' Normalize the loadings matrix \code{A}.
#'
#' @details
#' \code{type = "largest_abs"}: divide by largest absolute value of column
#'                        so that all values of mat are within interval [-1, 1];
#' \code{type = "orthonormal"}: obtain the orthonormal/singular vector decomposition matrix 
#'                        so that colSums(mat^2)=I;
#' \code{type = "L2"}: divide by the square root of column wise summation after square so that colSums(mat^2)=I;                        
#' \code{type = "none"}: keep original matrix without any standardization.
#' @param mat matrix; The matrix to transform.
#' @param type \code{"largest_abs"}; The performed standardization (could be extended with other ways of normalizing \code{A})
#' @return matrix; The transformed matrix.
normalize_loadings <- function(mat, type = "largest_abs") {
  # checks ---
  if (!is.matrix(mat)) {
    stop("mat is not a matrix.")
  }
  # if (!(type %in% c("largest_abs"))) {
  #   stop("type specified is not valid.")
  # }
  # (checks end) ---
  row_nms<-as.vector(rownames(mat))
  col_nms<-as.vector(colnames(mat))
  
  if (type == "largest_abs") {
    # increasing this to weaken effect of prior on z
    largest_abs_val <- 5
    # normalize A column-wise by the column (=factor) absolute maximums -> A will contain values between -1 and 1
    #   apply(...): computes absolute maximums of each column of A
    #   t(A) / ... : normalizes each column of A (division row-wise applied, that's why in the end final transformation)
    mat <- t(largest_abs_val * t(mat) / apply(mat, 2, function(v) max(abs(v), na.rm = T)))  # function t() returns transpose; na.rm means missing values are ignored
    mat <- sweep(mat, 2, apply(mat, 2, function(v) sign(v[which.max(abs(v))])), '*')
  } else if (type == "orthonormal"){
    ## find the orthonormal/singular Value decomposition matrix to get the colums(mat^2)=1
    mat<-as.data.frame(svd(mat)["u"])
  } else if (type == "L2"){
    ## find the matrix to get the colums(mat^2)=1
    denom<-apply(mat, 2, function(v) sqrt(sum((v^2), na.rm = T))) 
    mat_tmp<-NULL
    for (i in 1:ncol(mat))   {
      if(denom[i]==0){tmp<-rep(0,nrow(mat))} else {tmp<-mat[,i]/denom[i]}
      mat_tmp<-cbind(mat_tmp,tmp)
    } 
    mat<-mat_tmp
  } else if (type=="none"){
    mat<-mat  ## keep original matrix without any normalization
  } 
  
  rownames(mat)<-row_nms
  colnames(mat)<-col_nms
  mat<-as.matrix(mat)
  return(mat)
}


# === make A sparse with neuristic ===

#' @title
#' Sparsify a matrix.
#' 
#' @description
#' Set values below absolute threshold to 0.
#' 
#' @details
#' Heuristic to set all values in mat to 0 which are smaller in absolute value than thres
#' If combined with normalize_loadings(...) before:
#' Largest absolute value is 1 (values between -1 and 1), and
#' e.g. thres = 0.1 means everything smaller than 10% of highest absolute is set to 0
#'
#' @param mat matrix; The matrix to transform.
#' @param thres non-negative double; The absolute threshold below which values are set to 0.
#' @return mat matrix; The transformed matrix.
sparse_heur_thres <- function(mat, thres = 0.0) {
  if (!is.matrix(mat)) {
    stop("mat is not a matrix.")
  }
  if (!thres >= 0.0) {
    stop("thres must be non-negative double.")
  }
  mat[abs(mat) < thres] <- 0
  
  return(mat)
}



# === estimate A  ===

#' @title
#' Estimate \code{A}. 
#'
#' @description
#' Estimates the factor loadings matrix A with \emph{Probabilistic Principal Component Analysis (PPCA)}.
#' \code{A} is a (PxK) with the following properties:
#' \itemize{
#'   \item \code{A} is the loadings matrix, first resulting from a call to PPCA
#'     (and then further transformed as discussed below), thus,
#'     it is an orthogonal matrix which describes a basis of the original
#'     data in a lower-dimensional vector space.
#'   \item A is sparse through two operations:
#'      \enumerate{
#'         \item We apply a sparsity inducing rotation to it using the varimax/promax function.
#'         \item We set values below a certain absolute threshold to 0.
#'      }
#' }
#'
#' @details
#' When \code{type == "n_PCs"}, contrary to G.'s original script, we use the version
#' of estimating the n_PCs from the branch with var_explained, i.e. without the sign change,
#' mainly for consistency reasons (the user would not expect different operations in the two
#' versions except for how the PCs are selected.)
#'
#' Documentation for pcaMethods: https://www.bioconductor.org/packages/devel/bioc/manuals/pcaMethods/man/pcaMethods.pdf
#'
#' @param Y_in matrix; The data matrix to regress on. Dimensions: (NxP)
#' @param type \code{"var_explained"} or \code{"n_PCs"}: How to select the number of principal components.
#' @param var_explained double within interval (0.0, 1.0]: must be specified if \code{type == "var_explained"}, the amount of variance (cummulative R^2) that the chosen principal components corresponding to the largest eigenvalue explain at least.
#' @param n_PCs positive integer >= 2 or "max": must be specified if \code{type == "var_explained"}. 
#'               \itemize{
#'                   \item if integer: equals number of principal components used; note: maximum number of pinciple components chosen must be smaller than number of columns in Y_in
#'                   \item if "max": chooses the maximum possible number of principal components, usually the number of columns in Y_in minus 1
#'               }
#' @param rotation \code{"varimax"} or \code{"promax"}; the type of sparcity inducing rotation applied to the loadings matrix after PPCA.
#' @param conv_thres positive double; convergence threshold of PPCA.
#' @param ppca_seed positive integer; seed of PPCA
#'                                 }
#' @return list(A, ppca_object) list(matrix, ppca object); \code{A} is the transformed loadings matrix; ppca_object is the object from the ppca call.
#'                                                       \code{A} is a (PxK) matrix, where K is the number of chosen principal components.
estimate_loadings <- function(Y_in, 
                              type, 
                              var_explained = NULL, 
                              n_PCs = NULL, 
                              rotation = "promax", 
                              conv_thres = 1e-4, 
                              ppca_seed = 1) {
  # checks ---
  if (!is.matrix(Y_in)) {
    stop("Y_in is not a matrix.")
  }
  if (!type %in% c("var_explained", "n_PCs")) {
    stop("type must be either 'var_explained' or 'n_PCs'.")
  }
  # if type 'var_explained' selected -> var_exlained must be specified, n_PCs must be NULL
  # the following two checks must be done before checking for the values of n_PCs and var_explained
  if ((type == "var_explained") && (is.null(var_explained) || !is.null(n_PCs))) {
    stop("type is 'var_explained', i.e. only var_explained should be specified, n_PCs should be NULL.")
  }
  if ((type == "n_PCs") && (!is.null(var_explained) || is.null(n_PCs))) {
    stop("type is 'n_PCs', i.e. only n_PCs should be specified, var_explained should be NULL.")
  }
  # done after correct value checking
  if (type == "var_explained" && !( (var_explained > 0.0) && (var_explained <= 1.0) )) {
    stop("var_explained must be double in > 0.0 and <= 1.0.")
  }
  if (type == "n_PCs" && !(n_PCs == "max") && (!(n_PCs == as.integer(n_PCs)) || !(n_PCs >= 2) || !(n_PCs < dim(Y_in)[2]) )) {
    stop("n_PCs must be positive integer >= 2 and < P or 'max'.")
  }
  if (!rotation %in% c("varimax", "promax")) {
    stop("rotation must be either 'varimax' or 'promax'.")
  }
  if (!(conv_thres > 0.0)) {
    stop("conv_thres must be positive double.")
  }
  # if (!(sparse_abs_thres >= 0.0 && sparse_abs_thres <= 1.0)) {
  #   stop("sparse_abs_thres must be double between 0.0 and 1.0, including bounds.")
  # }
  if (!(is.na(ppca_seed)) && (!(ppca_seed == as.integer(ppca_seed)) || !(ppca_seed >= 1))) {
    stop("ppca_seed must be positive integer >= 1 or NA (random).")
  }
  # at least two variables must be present in group to reasonably estimate principal component of them
  if (ncol(Y_in) < 2) {
    stop("Less than two columns present in Y_in -> cannot reasonably estimate PCs.")
  }
  # if (!normalize_loadings_type %in% c("largest_abs", "none")) {
  #   stop("normalize_loadings_type is not valid.")
  # }
  # (checks end) ---
  
  # TODO check that at least two PCs selected
  # TODO: warning if sparsity_threshold is larger than largest absolute value -> no values removed
  # TODO: show info what fraction of values was set to 0 during sparcity inducing thing
  
  # empty result data types
  A <- NULL
  n_vars <- ncol(Y_in)
  
  # perform PPCA ---
  
  # Probabilistic PCA call -> https://www.rdocumentation.org/packages/pcaMethods/versions/1.64.0/topics/ppca
  # Returned object is of class pcaRes -> https://www.bioconductor.org/packages/devel/bioc/manuals/pcaMethods/man/pcaMethods.pdf (page 42)
  # model: Y =       U    %*%    t(V)
  #     (data)   (scores)    (loadings)
  # data passed to ppca call: (NxP)
  # ppca@loadings (V): (PxK)
  # ppca@scores (U): (NxK)
  # Note: dimensions of returned matrices are weird
  # TODO when could this try bracket throw an error?
  if (type == "var_explained") {
    tryer <- try({
      # perform ppca so to extract n_vars - 1 (maximum) pricipal components (at least one less than variables)
      # surpress warnings, since this call throws precision errors for some of the higher principal components
      suppressWarnings(ppca_object <- pcaMethods::ppca(Y_in, 
                                                       nPcs = n_vars - 1, 
                                                       seed = ppca_seed, 
                                                       threshold = conv_thres))
    })
    if(inherits(tryer, "try-error")){
      if(n_vars > 2) {
        # TODO why does -2 helps?
        # surpress warnings, since this call throws precision errors for some of the higher principal components
        suppressWarnings(ppca_object <- pcaMethods::ppca(Y_in, 
                                                         nPcs = n_vars - 2, 
                                                         seed = ppca_seed, 
                                                         threshold = conv_thres))
      }
      else {
        # next  # TODO what to do in this case?
        stop("PPCA error - Related to number of variables and maximum number of principal components.")
      }
    }
    # the minimum number of principal components (with largest eigenvalues) that explain n_PCs (fraction) of the variance
    # if (type == "var_explained") {
    n_PCs_keep <- min(c(ncol(ppca_object@loadings), 
                        sum(ppca_object@R2cum < as.numeric(var_explained)) + 1))
    # } else if (type == "n_PCs" && n_PCs == "max") {
    #   n_PCs_keep = dim(ppca_object@loadings)[2]  # keep all principal components estimated
    # }
    # 'n_PCs' and integer
    # else {
    #   n_PCs_keep = n_PCs
    # }
    # select the principal components (from the loadings matrix) corresponding to these principal components kept
    V <- ppca_object@loadings[, 1:n_PCs_keep, drop = F]
  } else if (type == "n_PCs") {
    if (n_PCs == "max") {
      tryer <- try({
        # perform ppca so to extract n_vars - 1 (maximum) pricipal components (at least one less than variables)
        # surpress warnings, since this call throws precision errors for some of the higher principal components
        suppressWarnings(ppca_object <- pcaMethods::ppca(Y_in, 
                                                         nPcs = n_vars - 1, 
                                                         seed = ppca_seed, 
                                                         threshold = conv_thres))
      })
      if(inherits(tryer, "try-error")){
        if(n_vars > 2) {
          # TODO why does -2 helps?
          # surpress warnings, since this call throws precision errors for some of the higher principal components
          suppressWarnings(ppca_object <- pcaMethods::ppca(Y_in, 
                                                           nPcs = n_vars - 2, 
                                                           seed = ppca_seed, 
                                                           threshold = conv_thres))
        }
        else {
          # next  # TODO what to do in this case?
          stop("n_PCs == 'max' didn't work. Try a (small) integer instead.")
        }
      }
    }
    else {
      ppca_object <- pcaMethods::ppca(Y_in, 
                                      nPcs = n_PCs, 
                                      seed = ppca_seed, 
                                      threshold = conv_thres)
    }
    V <- ppca_object@loadings
  }
  # perform rotation ---
  # perform sparcity inducing rotation  (not done: , if reasonable (more than one principal component kept))
  if(!(rotation == "none")) {  # n_PCs_keep > 1 &&
    Vrotraw <- sparse_rotation(V, rotation)  # perform the sparcity inducing rotation
  } else {
    # perform no rotation
    Vrotraw <- V
  }
  
  A <- Vrotraw  # A is the loadings matrix we want to estimate
  
  # # normalize A
  # # A is PxK matrix
  # if (normalize_loadings_type != "none") {
  #   A <- normalize_loadings(mat = A, type = normalize_loadings_type)
  # }
  # 
  # # make A more sparse with heuristic threshold
  # A <- sparse_heur_thres(mat = A, thres = sparse_abs_thres)
  
  return(list(A, ppca_object))
}


# === find unique missing patterns ===

#' Search for unique NA patterns in rows of input matrix.
#'
#' @param Y_in matrix; The matrix in which to find the unique missingness patterns.
#' @return list(unique_missing_patterns : list, indices_Y_in_list : list);
#'              \itemize{
#'                \item unique_missing_patterns: The unique missingness patterns found (as strings).
#'                \item indices_Y_in_list: list of elements where each element is list of indices in Y_in which have this missingness pattern.
#'              }
find_unique_missing_patterns <- function(Y_in) {
  # checks ---
  if (!is.matrix(Y_in)) {
    stop("Y_in is not a matrix.")
  }
  # (checks end) ---
  all_missing_patterns_num <- is.na(Y_in)  # same dimensions as Y_in, where every element indicates whether the corresponding element in Y is NA or not -> values: TRUE or FALSE
  mode(all_missing_patterns_num) = "numeric"  # compare ?mode: sets the storage mode of the object misall to numeric -> values: 1 or 0 (which indicate NA or not)
  all_missing_patterns_char <- NULL  # empty variable that is filled with the missing pattern of every row in all_missing_patterns_num
  for(i in 1:ncol(all_missing_patterns_num)) {
    # string concatenation of all columns of all_missing_patterns_num -> concatenates in vectorized fashion
    # i.e. all_missing_patterns_char has length of dim(all_missing_patterns_num)[1] and each element in all_missing_patterns_char correponds to the concatenated string of 0s and 1s in each row of all_missing_patterns_num
    all_missing_patterns_char <- paste0(all_missing_patterns_char, all_missing_patterns_num[, i])
  }
  # redefine rownames for internal processing only -> give each row unique ID
  # names(all_missing_patterns_char) <- rownames(all_missing_patterns_num)  # assign rownames of all_missing_patterns_num (e.g. specific ID or just numbers) to all_missing_patterns_char
  names(all_missing_patterns_char) <- c(1:length(all_missing_patterns_char))  # names are strings
  all_missing_patterns_char <- sort(all_missing_patterns_char)  # sorts all_missing_patterns_char in ascending order in terms of the missing pattern (i.e. string entry with most 0s first)
  unique_missing_patterns_char <- unique(all_missing_patterns_char)  # all unique missing patterns
  # the "start" index in all_missing_patterns_char of entries which have the same missing pattern
  # the positions (integers) in the sorted missing pattern vector all_missing_patterns_char where a unique missing pattern appears first (-> thus, an ascending integer vector)
  unique_missing_patterns_char_st <- match(unique_missing_patterns_char, all_missing_patterns_char)
  # the "end" index in all_missing_patterns_char of entries which have the same missing pattern
  # positions like unique_missing_patterns_char.st, but 1) leaving out first entry 2) adding last entry to be last position in all_missing_patterns_char 3) shiftig all other positions by 1 down
  if (length(unique_missing_patterns_char) > 1) {
    # regular case: at least two unique missing patterns
    unique_missing_patterns_char_en <- c(unique_missing_patterns_char_st[2:length(unique_missing_patterns_char_st)] - 1, length(all_missing_patterns_char))
  } else {
    # special case: only one unique missing pattern
    unique_missing_patterns_char_en <- c(length(all_missing_patterns_char))
  }
  # subj.t.keep <- rownames(Y)  # TODO why named keep?
  
  # list of list, where each lement list is rownames (row indices) that have same unique_missing_pattern
  indices_Y_in_list <- NULL  # empty result
  for (i in 1:length(unique_missing_patterns_char_st)) {
    # find start and end index of unique missing pattern rows in unique_missing_patterns_char (sorted)
    index_sorted_start <- unique_missing_patterns_char_st[i]
    index_sorted_end <- unique_missing_patterns_char_en[i]
    
    # find rownames -> indices in original unsorted Y_in that correspond to this unique missing patterns
    # since names are strings -> convert to numeric
    indices_Y_in_list <- c(indices_Y_in_list, list(names(all_missing_patterns_char)[index_sorted_start:index_sorted_end]))
  }
  
  return(list(unique_missing_patterns_char, indices_Y_in_list))
}


# === normalize_stage_1_2_results ===

#' @title
#' Normalize \code{A}, \code{Z} and \code{S}. 
#' 
#' @description
#' Normalize the result matrices \code{A} (stage 1), and \code{Z} and \code{S} (stage 2) per factor ("K dimension").
#'
#' @param A_in matrix; The loadings matrix.
#' @param Z_in matrix; The scores mean matrix.
#' @param S_in matrix: The scores standard error matrix.
#' @param type = "normalize_custom" : string; The type of normalization to perform.
#'          \itemize{
#'             \item \code{"normalize_custom"}: Normalize by statistics computed on \code{Z} and \code{S} (primary approach).
#'             \item \code{"scale"}: Scale each matrix (\code{A}, \code{Z}, \code{S}) (mean = 0, std = 1; column-wise)
#'          }
#' @return \code{list(A_sc, Z_sc, S_sc)}: List of \code{A_in}, \code{Z_in} and \code{S_in} after normalizing.
normalize_stage_1_2_results <- function(A_in, Z_in, S_in, type = "normalize_custom") {
  # checks ---
  if (!is.matrix(A_in)) {
    stop("A_in is not a matrix.")
  }
  if (!is.matrix(Z_in)) {
    stop("Z_in is not a matrix.")
  }
  if (!is.matrix(S_in)) {
    stop("S_in is not a matrix.")
  }
  if (!(type %in% c("normalize_custom", "scale"))) {
    stop("type specifies no valid normalization.")
  }
  # (checks end) ---
  if (type == "normalize_custom") {
    wmn <- colSums(Z_in / S_in^2, na.rm = T) / colSums(1 / S_in^2, na.rm = T)
    wsd <- sqrt(colSums(t(t(Z_in) - wmn)^2 / S_in^2, na.rm = T) / colSums(1 / S_in^2, na.rm = T))
    
    Z_sc <- t(t(Z_in) / wsd)
    S_sc <- t(t(S_in) / wsd)
    A_sc <- t(t(A_in) / wsd)
  }
  else if (type == "scale") {
    Z_sc <- scale(Z_in)
    S_sc <- scale(S_in)
    A_sc <- scale(A_in)
  }
  
  return(list(A_sc, Z_sc, S_sc))
}

# === estimate scores ===

#' @title
#' Estimate Z and S. 
#' 
#' @description
#' This function esitmates the scores matrix \code{Z} and its standard error matrix \code{S} corresponding to a loadings matrix A_in with \emph{least-squares linear regression}.
#' 
#' @details
#' A_in is the output of a call to estimate_loadings(...).
#'
#' Since we performed further sparcity inducing operations on the loadings matrix A of the PPCA call,
#' we cannot use the scores matrix of the resulting ppca_object directly. Therefore, this method estimates
#' the scores and their standard errors using least-squares linear regression. This is done in a
#' computationally efficient way by abusing the missingness pattern in \code{A_in} (rows with the same missing pattern
#' are estimated jointly).
#'
#' Further extensions:
#' \itemize{
#'   \item parallelize the fitting of rows with the same missingness pattern
#' }
#'
#' @param Y_in matrix; The original data marix to regress on. Dimensions: (NxP).
#' @param A_in matrix; The loadings matrix. Dimensions: (PxF).
#' @param type string: \code{"lin_reg"}: Perform a linear regression to estimate Z and S (currently only this type implemented).
#' @return \code{list(Z : matrix, S : matrix)}; The scores matrix (\code{Z} is mean, \code{S} is standard error).
#'         \itemize{
#'           \item \code{Z}: Dimensions: (NxK)
#'           \item \code{S}: Dimensions: (NxK)
#'         }
estimate_scores <- function(Y_in, A_in, type = "lin_reg") {
  # checks ---
  if (!is.matrix(Y_in)) {
    stop("Y_in is not a matrix.")
  }
  if (!is.matrix(A_in)) {
    stop("A_in is not a matrix.")
  }
  # if (!type %in% c("lin_reg")) {
  #   stop("Variable type must be 'lin_reg'.")
  # }
  if (!(dim(A_in)[1] == dim(Y_in)[2])) {
    stop("P dimension, which is first dimension of A_in and second dimension of Y_in, must match.")
  }
  # (checks end) ---
  #Calculate Z for each pattern of missingness in turn
  list[unique_missing_patterns_char, indices_Y_in_list] <- find_unique_missing_patterns(Y_in)
  # the result matrices (Z and S) of estimate_Z.R, initialized with NAs
  # Z: factor scores (weight estimates for each row)
  # S: standard deviations of weights (same for each row with same missingness pattern)
  # dimensions of Z: T by K; where T is number of (subject, visit) pairs, K is the number of principal components
  # TODO check dimensions
  # Z, S are reversed to transpose result Zc below! -> that's why transposed is assigned -> TODO make more consistent
  # Z, S has dimensions: (NxK)
  Z <- S <- matrix(NA, dim(Y_in)[1], dim(A_in)[2], dimnames = list(1:dim(Y_in)[1], 1:dim(A_in)[2]))
  
  # estimate the factor score for each (subject, visit) pair (= each row of Y) separately
  # IMPORTANT NOTE: if the missing pattern in Y is the same for any rows, we can parallelize the estimation of the factor scores
  #                 since the non-missing variables in such rows are the same, only they have different values.
  #                 That's why we loop over the unique missing value patterns of Y instead of the rows of Y directly.
  names(Y_in) <- c(1:nrow(Y_in))  # will be used for addressing
  dimnames(Y_in) <- list(1:dim(Y_in)[1], 1:dim(Y_in)[2])  # will be used for indexing (dimnames required)
  
  if (type == "a_posteriori") {
    # new version from George
    Y_cov_mat <- cov(Y_in, use = "pairwise.complete.obs")
    Y_cov_mat[is.na(Y_cov_mat)] <- 0
    eigenvalues_lambda <- eigen(Y_cov_mat)$values
    # see (12.46) in Bishop
    K <- dim(A_in)[2]
    sigma_square_ml <- mean(eigenvalues_lambda[(K + 1):length(eigenvalues_lambda)]) 
  }

  for(indices_Y_in in indices_Y_in_list){
    # pattern.inds <- names(missall.bin)[missun.st[missc]:missun.en[missc]]  # all subj.t identifiers (rownames of missall / Y) with the same missing pattern
    # Y_smp has rows with the same missingness pattern in Y_in
    Y_smp <- Y_in[indices_Y_in, , drop = F]  # drop = F ensures that the matrix is not flattened, even if pattern.inds is only referring to one row
    non_NA_indices_smp <- which(!is.na(Y_smp[1, ]))  # boolean vector referring to non-NA values in Y with the pattern missun[missc]
    # equ:   y_it =    A    z_it + e_it
    # dims:  (PxF)   (PxK)  (KxF)  (PxF)
    # A: output of estimate_loadings() -> is "already transposed" as first equation in stage 2.
    # P: number of original variables (!) which are not NA
    # K: number of principal components
    # F: number of rows in Y_in with this missingness pattern
    # thus, if P is reduced -> subselect COLUMNS of Y_smp (corresponds to row entries of y_it) and subselect ROWS of A
    Y_smp_non_NA <- Y_smp[, non_NA_indices_smp, drop = F]  # only select COLUMNS in Y_smp with non-NA values
    A_smp <- A_in[non_NA_indices_smp, , drop = F]  # only select ROWS in A which correspond to the non-NA values in Y_smp
    
    # --- POTENTIAL HELPER START
    if (type == "lin_reg") {
      At_A_inv <- MASS::ginv(t(A_smp) %*% A_smp)  # At_A_inv is symmetric, i.e. transpose does not change calculation
      dimnames(At_A_inv) <- list(colnames(A_smp), colnames(A_smp))  # name dimensions of At_A_inv according to column names of A
      Zc <- At_A_inv %*% t(A_smp) %*% t(Y_smp_non_NA)  # equation (2); TODO why not swapping the order of the terms according to transpose of the term?
      eps.tol <- 1e-10
      diag(At_A_inv)[abs(diag(At_A_inv)) < eps.tol] <- 0  ## Deal with small negative numbers caused by rounding errors -> sets them to 0
      Zc.se <- sqrt(diag(At_A_inv))  # standard  deviation of k-th factor is sqrt of k-th diagonal element of At_A_inv (see formula s_it(k) in stage 3)
      if(any(!is.finite(sqrt(diag(At_A_inv))))) {
        stop("Non-finite value on diagonal of ginv(t(A) %*% A).")  # if non-finite standard deviations estimated: stop  # TODO when might this happen? how to handle?
      }
    }
    else if (type == 'a_posteriori') {
      # Note: this implementation has another part above! check if (type == "a_posteriori") ... there.
      M <- t(A_smp) %*% A_smp
      M <- M + sigma_square_ml * diag(dim(M)[1])
      dimnames(M) <- list(colnames(A_smp), colnames(A_smp))  # name dimensions of M according to column names of A
      Zc <- MASS::ginv(M) %*% t(A_smp) %*% t(Y_smp_non_NA)  # equation (2); TODO why not swapping the order of the terms according to transpose of the term?
      eps.tol <- 1e-10
      At_A_inv <- MASS::ginv(M)
      diag(At_A_inv)[abs(diag(At_A_inv)) < eps.tol] <- 0  ## Deal with small negative numbers caused by rounding errors -> sets them to 0
      Zc.se <- sqrt(sigma_square_ml*diag(At_A_inv))  # standard  deviation of k-th factor is sqrt of k-th diagonal element of At_A_inv (see formula s_it(k) in stage 3)
      if (any(!is.finite(sigma_square_ml*sqrt(diag(At_A_inv))))) {
        stop("Non-finite value on diagonal of ginv(t(A) %*% A).")  # if non-finite standard deviations estimated: stop  # TODO when might this happen? how to handle?
      }
    }
    # --- POTENTIAL HELPER END
    
    facs.with.no.data <- colnames(A_smp)[colSums(A_smp != 0) == 0]  # columns (= factors) of A_smp which have at least one non-zero (before: NA) value
    Zc[facs.with.no.data, ] <- Zc.se[facs.with.no.data] <- NA  # set factors which do not have at least one non-zero value to NA -> computed mean and standard deviation is meaningless
    # assign computed
    # Zc has dimensions (K, N_ind), where N_ind is the number of observations with the same missing pattern
    Z[indices_Y_in, ] <- t(Zc)  # the estimated factor scores for every row separately
    S[indices_Y_in, ] <- rep(Zc.se, each = length(indices_Y_in))  # standard deviation is the same for each row with certain missing pattern
  }
  
  return(list(Z = Z, S = S, sigma_square_ml = sigma_square_ml))
}


#' @title 
#' Remove all NA rows
#' 
#' @description
#' Remove rows that contain only NAs in \code{Y_loadings}/\code{Y_scores}. 
#' The colums considered are those in \code{Y_loadings}/\code{Y_scores} (\code{vars} columns), 
#' and their corresponding rows are also removed in \code{data_loadings}/\code{data_scores}, which contain further columns, 
#' but the same number of rows.
#' 
#' @param Y_loadings matrix; The input matrix for estimating loadings.
#' @param Y_scores matrix; The input matrix for estimating scores.
#' @param data_loadings tibble; data frame containing data for loadings estimation.
#' @param data_scores tibble; data frame containing data for scores estimation.
#' @return \code{list(Y_loadings, Y_scores, data_loadings, data_scores)} with all NA rows (of \code{vars} columns) removed.
remove_all_NA_rows <- function(Y_loadings, Y_scores, data_loadings, data_scores) {
  # checks ---
  if (!is.matrix(Y_loadings)) {
    stop("Y_loadings is not a matrix.")
  }
  if (!is.matrix(Y_scores)) {
    stop("Y_scores is not a matrix.")
  }
  if (!tibble::is_tibble(data_loadings)) {
    stop("data_loadings is not a tibble.")
  }
  if (!tibble::is_tibble(data_scores)) {
    stop("data_scores is not a tibble.")
  }
  # (checks end) ---
  # if there are rows will only missing values in the data, these must be removed,
  # otherwise, in the "normalize_custom" version of normalize_A_Z_S.type,
  # S values will be 0 and this will cause NA/inf values, resulting in errors.
  if ((length(which(rowSums(is.na(Y_loadings) == 0) == 0)) > 0) || (length(which(rowSums(is.na(Y_scores) == 0) == 0)) > 0)) {
    # count before
    rows_before_loadings <- dim(Y_loadings)[1]
    rows_before_scores <- dim(Y_scores)[1]
    # remove rows which contain any missing values in 4 objects: Y_loadings, Y_scores, data_loadings, data_scores
    at_least_one_non_NA_rows_loadings <- !(rowSums(is.na(Y_loadings) == 0) == 0)
    at_least_one_non_NA_rows_scores <- !(rowSums(is.na(Y_loadings) == 0) == 0)
    Y_loadings <- Y_loadings[at_least_one_non_NA_rows_loadings, ]
    Y_scores <- Y_scores[at_least_one_non_NA_rows_scores, ]
    data_loadings <- data_loadings[at_least_one_non_NA_rows_loadings, ]
    data_scores <- data_scores[at_least_one_non_NA_rows_scores, ]
    # count after
    rows_after_loadings <- dim(Y_loadings)[1]
    rows_after_scores <- dim(Y_scores)[1]
    if ((rows_before_loadings-rows_after_loadings)!=(rows_before_scores-rows_after_scores)) {
      stop("Different number of rows removed from data_loadings and data_scores.")
    }
    warning(paste0(as.character(rows_before_loadings-rows_after_loadings), " rows with NAs in all vars columns in data. These rows are removed in the internal processing of mfd."))
  }
  return(list(Y_loadings, Y_scores, data_loadings, data_scores))
}
