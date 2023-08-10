#' @title
#' For returning lists of objects.
#' 
#' @description 
#' Together with \code{library(gsubfn)}, this function allows to assign to a list of values multiple return values from a function,
#' e.g. \code{list[a, b] <- functionReturningTwoValues()}, where \code{functionReturningTwoValues()} returns \code{list(a, b)}
#' 
#' @details
#' Compare: https://stackoverflow.com/questions/1826519/how-to-assign-from-a-function-which-returns-more-than-one-value
#' 
#' Source code copied from: https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
#' @param x unknown.
#' @param value unknown.
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}
list <- structure(NA,class="result")   # put inside function. correct?

#' Combine a list of lists into a single list.
#'
#' @param list_of_lists a list of lists.
#' @return combined_list a list.
combine_lists_as_one <- function(list_of_lists) {
  combined_list <- c()
  for (sublist in list_of_lists) {
    combined_list <- c(combined_list, sublist)
  }
  
  return(combined_list)
}


# remove rows and columns ==========================

#' @title
#' Remove rows and columns with high NA ratio.
#' 
#' @description
#' Remove both rows and columns which don't have a minimum fraction of non-NA values, only considering the \code{columns} specified.
#' 
#' @details
#' Ensures that the data is both in rows and columns not too sparse.
#' This helper function can be used during the preprocessing of the data.
#'
#' Could have been implemented with dplyer: https://sebastiansauer.github.io/sum-isna/, but tested and much slower.
#'
#' @param data tibble; The data to transform.
#' @param columns list of strings; The columns to which this function applies.
#' @param min_nonNA_fraction_per_col double in [0, 1]; Minimum fraction of non-\code{NA} values per column (applied to all \code{columns}).
#' @param min_nonNA_fraction_per_row double in [0, 1]; Minimum fraction of non-\code{NA} values per row (applied to all \code{columns}).
#' @return data tibble; The data to transform.
remove_obs_and_meas_NA_threshold <- function(data, columns, min_nonNA_fraction_per_row = 0.0, min_nonNA_fraction_per_col = 0.0) {
  # checks ---
  # check: data is tibble
  if (! class(data)[1] == "tbl_df") {
    stop("data is not a tibble.")
  }
  # check: columns must be existent
  if (is.null(columns)) {
    stop("Variable columns must not be NULL.")
  }
  if (!all(columns %in% colnames(data))) {  # all columns must be in data (logical and on vector)
    stop("Variable columns specifies columns in data that are not existent.")
  }
  if (!(min_nonNA_fraction_per_row >= 0.0 && min_nonNA_fraction_per_row <= 1.0 && min_nonNA_fraction_per_col >= 0.0 && min_nonNA_fraction_per_col <= 1.0)) {
    stop("Variables min_nonNA_fraction_per_row or min_nonNA_fraction_per_col not in continuous interval [0.0; 1.0].")
  }
  # (checks end) ---
  
  # only specified columns are used -> just used for internal processing
  data_cols <- dplyr::select(data, columns)
  
  # counting non-NA values
  row_sums <- rowSums(!is.na(data_cols))
  col_sums <- colSums(!is.na(data_cols))
  removed_rows <- which(!(row_sums/(dim(data_cols)[2]) >= min_nonNA_fraction_per_row))
  removed_cols <- which(!(col_sums/dim(data_cols)[1] >= min_nonNA_fraction_per_col))
  
  # rows are removed and warning
  if (length(removed_rows > 0)) {
    # warning(paste(c(length(removed_rows), " rows are removed, since they do not contain the minimum fraction (min_nonNA_fraction_per_row = ", min_nonNA_fraction_per_row, ") non-NA values.")))
    data <- dplyr::slice(data, -removed_rows)
  }
  
  # columns are removed and warning
  if (length(removed_cols) > 0) {
    # warning(paste(c("The following columns are removed, since they do not contain the minimum fraction (min_nonNA_fraction_per_col = ", min_nonNA_fraction_per_col, ") non-NA values: \n", paste(names(removed_cols), collapse=", " ))))
    data <- dplyr::select(data, -removed_cols)
  }
  
  return(data)
}

# === stratified sampling helper ===

#' @title
#' Stratified sampling
#' 
#' @description
#' Perform stratified subsampling in subgroups of the population according to a subsampling strategy.
#' 
#' @details
#' Similar behavior like function \code{stratified}, but ensuring a given total sample size.
#'
#' Possible extensions:
#' \itemize{
#'   \item extend to multiple variables in \code{group} and \code{sample_unit}
#' }
#'
#' @param data tibble; The data to transform.
#' @param group string; column in data indicating the subgroup a row belongs to.
#' @param total_sample_size integer; The size of the total sample to be produced.
#' @param sample_unit string; Variables in data indicating what to sample (sampling "rows" by default with value \code{NULL}).
#' @param subgroup_sample_strat \code{"equal"} or \code{"proportionate"}; The sampling strategy to use in the subgroups. 
#'                              \itemize{
#'                                \item \code{"equal"} chooses an equal number of sample units per subgroup, 
#'                                \item \code{"proportionate"} chooses rows according to the fraction of sample units from this subgroup in \code{data}.
#'                              }
#' @param with_replace boolean; whether to sample with or without replacement
#' @return data tibble; The transformed data.
stratified_sampling <- function(data, group, total_sample_size, sample_unit = NULL, subgroup_sample_strat = "proportionate", with_replace = FALSE) {
  # checks ---
  if (! class(data)[1] == "tbl_df") {
    stop("data is not a tibble.")
  }
  if (!all((group %in% colnames(data)))) {  # all columns must be in data (logical and on vector)
    stop("variable group specifies columns not in data")
  }
  if (!( (as.integer(total_sample_size) == total_sample_size) && (total_sample_size >= 1) && (total_sample_size >= length(group)) )) {
    stop("total_sample_size must be positive int and with at least one sample per group.")
  }
  if ( (!is.null(sample_unit)) && (!all(sample_unit %in% colnames(data))) ) {
    stop("Specified sample_unit not in columns of data (set to NULL for rows).")
  }
  if (!(subgroup_sample_strat %in% c("proportionate", "equal"))) {
    stop("subgroup_sample_strat chosen incorrectly, must be 'proportionate' or 'equal'.")
  }
  if (!(typeof(with_replace) == "logical")) {
    stop("with_replace must be boolean.")
  }
  # number of units that can be sampled from
  if (is.null(sample_unit)) {
    # all rows are sample units
    n_total_sample_units <- nrow(data)
  } else {
    n_total_sample_units <- nrow(dplyr::distinct(data, !!dplyr::sym(sample_unit)))  # TODO implement as for n_total_sample_units
  }
  if (n_total_sample_units < total_sample_size) {
    stop(paste0("total_sample_size (", total_sample_size, ") is smaller than number of sampling units in population (", n_total_sample_units, ")"))
  }
  # TODO for each group, check if enough sample units are there -> depending on sampling strategy
  # (checks end) ---
  
  # number of groups to sample from
  n_groups <- dim(dplyr::distinct(data[, group]))[1]  # alternative way: data %>% select(group) %>% dplyr::n_distinct()
  # vector to be filled with number of samples per group, initialized with zeros
  n_samples_per_group <- as.integer(n_groups)
  
  distinct_subgroups <- dplyr::distinct(data, !!dplyr::sym(group))
  distinct_subgroups <- dplyr::pull(distinct_subgroups, !!(dplyr::sym(group)))  # as vector
  
  if (subgroup_sample_strat == "proportionate") {
    i <- 1
    # n_total_sample_units <- select(data, sample_unit) %>% dplyr::n_distinct()  # already done above
    for (gr in distinct_subgroups) {
      # special case to ensure the correct size of the total sample
      if (i == n_groups) {
        n_samples_per_group[n_groups] <- total_sample_size - sum(n_samples_per_group[1:n_groups-1])
      } else {
        group_ratio <- dplyr::filter(data, !!dplyr::sym(group) == gr)
        if (!is.null(sample_unit)) {
          group_ratio <- group_ratio %>% dplyr::select(sample_unit) %>% dplyr::n_distinct() / n_total_sample_units
        }
        else {
          group_ratio <- group_ratio %>% dplyr::n_distinct() / n_total_sample_units
        }
        
        n_samples_per_group[i] <- floor(group_ratio * total_sample_size)
      }
      i <- i + 1
    }
  }
  else if (subgroup_sample_strat == "equal") {
    n_samples_per_group[1:n_groups-1] <- floor(total_sample_size / n_groups)
    n_samples_per_group[n_groups] <- total_sample_size - sum(n_samples_per_group[1:n_groups-1])  # last group ensures that total sample size is correct
  }
  
  # do the sampling ---
  data_sampled <- NULL  # result dataframe with subsampled data
  i <- 1
  for (gr in distinct_subgroups) {
    data_sampling_subgroup <- dplyr::filter(data, !!dplyr::sym(group) == gr)
    if (!is.null(sample_unit)) {
      sampled_units <- sample(dplyr::pull(dplyr::distinct(data_sampling_subgroup[, sample_unit]), !!(dplyr::sym(sample_unit))), n_samples_per_group[i], replace = with_replace)
      data_add <- dplyr::filter(data_sampling_subgroup, !!(dplyr::sym(sample_unit)) %in% sampled_units)
      data_sampled <- dplyr::bind_rows(data_sampled, data_add)
      i <- i+1
    }
    else {
      sampled_units <- sample(dplyr::pull(dplyr::distinct(data_sampling_subgroup)), n_samples_per_group[i], replace = with_replace)
    }
    
  }
  
  # ungroup the data
  data_sampled <- data_sampled %>% dplyr::ungroup()
  
  return(data_sampled)
}


# === sampling ===

#' @title
#' Sample the rows of the data.
#' 
#' @description
#' Wrapper function for sampling of subjects from data.
#' Uses helper function \code{stratified_sampling(...)}.
#'
#' @param data tibble; The data to sample from.
#' @param series_group string; The variable that indicates series to sample (e.g. \code{"subj_ID"}).
#' @param subgroup string; The variable that indicates the subgroups from which to sample.
#' @param total_sample_size The total number of series to draw.
#' @param type \code{"equal_n_series_from_subgroups"} or \code{"proportionate_n_series_from_subgroups"}: The type of stratified sampling to perform.
#'             \itemize{
#'               \item \code{"equal_n_series_from_subgroups"}: Samples an equal amount of series from each subgroup.
#'               \item \code{"proportionate_n_series_from_subgroups"}: Samples amount of series from each indicatino proportionate to their fraction among all subgroups.
#'             }
#' @return data tibble; The sampled data (with less series).
sampling <- function(data, series_group, subgroup, total_sample_size, type = "equal_n_series_from_subgroups") {
  # checks ---
  if (! class(data)[1] == "tbl_df") {
    stop("data is not a tibble.")
  }
  if (!( (as.integer(total_sample_size) == total_sample_size) && (total_sample_size >= 1) )) {
    stop("total_sample_size must be positive int.")
  }
  if (!type %in% c("equal_n_series_from_subgroups", "proportionate_n_series_from_subgroups")) {
    stop("type specification is not valid.")
  }
  if (!(all(c(series_group, subgroup) %in% colnames(data)))) {
    stop("Columns series_group or subgroup not present in data.")
  }
  # (checks end) ---
  if (type == "equal_n_seres_from_subgroups") {
    data <- stratified_sampling(data = data, group = subgroup, total_sample_size = total_sample_size, sample_unit = series_group, subgroup_sample_strat = "equal")
  }
  else if (type == "proportionate_n_subj_from_subgroups") {
    data <- stratified_sampling(data = data, group = subgroup, total_sample_size = total_sample_size, sample_unit = series_group, subgroup_sample_strat = "proportionate")
  }
  
  return(data)
}


# === convert_factor_to_numeric ===

#' @title 
#' Convert the (factor) input data to numeric data.
#' 
#' @description
#' Convert any factor columns (unordered or ordered categorical) to numeric. 
#' 
#' @details
#' The following conversions are conducted: 
#' - ordered factors -> converted to [0, 1, 2, 3, ...] (assumption of equal difference between category values)
#' - unordered factor (binary, unordered categorical with or without baseline) -> converted to treatment contrasts matrix
#' - numeric, character and other columns are not converted (simply binded to the result)
#' 
#' @param data tibble; the data that contains factor columns to convert to a purely numeric tibble
#' @return data_conv tibble; converted data tibble, with potentially more columns than data, but each (set of) columns representing one original column in data
convert_factor_to_numeric <- function(data) {
  # checks ---
  if (!tibble::is_tibble(data)) {
    stop("data is not a tibble.")
  }
  # (checks end) ---
  cols <- colnames(data)
  col_types <- sapply(data, class)
  data_conv <- NULL
  for (i in 1:length(cols)) {
    col <- cols[i]
    col_type <- col_types[[i]]  # e.g. "character", c("ordered", "factor")
    data_col <- data %>% select(col)
    # case: factor
    if ("factor" %in% col_type) {
      if (length(unique(na.omit(data_col))) > 20) {
        warning(paste0(col, " has ", as.character(length(unique(data_col))), " levels. Are you sure it is a factor?"))
      }
      # sub-case: ordered factor
      if ("ordered" %in% col_type) {
        print(paste0("Preprocessing ordered factor column ", col, " ..."))
        data_col_ord_cat <- convert_ordered_categorical_to_numeric(col, data_col %>% pull())
        data_conv <- dplyr::bind_cols(data_conv, data_col_ord_cat)
      }
      # sub-case: unordered factor (binary, unordered categorical with or without baseline)
      else {
        print(paste0("Preprocessing unordered factor (binary, unordered categorical with or without baseline) column ", col, " ..."))
        data_col_unord_cat <- dplyr::as_tibble(stats::model.matrix(~., data=(stats::model.frame(~., data=data_col, na.action=NULL)) ))
        # remove (Intercept) column
        data_col_unord_cat <- data_col_unord_cat %>% select(-"(Intercept)")
        data_conv <- dplyr::bind_cols(data_conv, data_col_unord_cat)
      }
    }
    # case: numeric or character or something else
    else {
      data_conv <- dplyr::bind_cols(data_conv, data_col)
    }
  }
  return(data_conv)
}


# === convert_ordered_categorical_to_numeric ===

#' @title 
#' Ordered categorical to numeric.
#' 
#' @description 
#' Convert an ordered categorical variable column to a numeric representation.
#' Maps an ordered categorical varialbe to [0, 1, 2, 3, ...].
#' 
#' @details
#' This assumes that the difference between category value are equal.
#' 
#' @param col character; the name of the column to convert
#' @param data_col factor vector; the data column to convert 
#' @return data_col_conv tibble; the converted data column
convert_ordered_categorical_to_numeric <- function(col, data_col) {
  # checks ---
  if (!is.character(col)) {
    stop("col is not a character.")
  }
  if (!is.factor(data_col) || !("ordered" %in% class(data_col))) {
    stop("data_col is not an ordered factor.")
  }
  # (checks end) ---
  data_col_conv <- data_col
  levels(data_col_conv) <- as.character(0:(length(levels(data_col_conv))-1))
  data_col_conv <- as.numeric(data_col_conv)-1  # added one, which we subtract again
  data_col_conv <- tibble::tibble(data_col_conv)
  colnames(data_col_conv) <- col
  return(data_col_conv)
}


# === convert_unordered_categorical_to_numeric

#' @title 
#' Make one-hot encoding.
#' 
#' @description 
#' Converts an unordered categorical variable to a numeric representation.
#' 
#' @details
#' Converts an unordered categorical variable to a one-hot encoding.
#' 
#' @param col character; the name of the column to convert
#' @param data_col factor vector; the data column to convert 
#' @return data_col_conv tibble; the converted data column
convert_unordered_categorical_to_numeric <- function(col, data_col) {
  # (checks) ---
  if (!is.character(col)) {
    stop("col is not a character.")
  }
  if (!is.factor(data_col) || ("ordered" %in% class(data_col))) {
    stop("data_col is not an unordered factor.")
  }
  # (checks end) ---
  # create a one-hot encoding
  data_col_conv <- stats::model.matrix(~0+data_col)   # stats::as.formula(paste0("~ 0 + ", col))
  column_names <- colnames(data_col_conv)  # data_colblue, data_colgreen etc.
  column_names <- str_remove(column_names, "data_col")
  column_names <- paste0(col, "_", column_names)
  colnames(data_col_conv) <- column_names
  data_col_conv <- as.data.frame(data_col_conv)
  data_col_conv[column_names] <- sapply(data_col_conv[column_names], as.numeric)
  return(data_col_conv)
}


# === convert_unordered_categorical_with_baseline_to_numeric

# convert_unordered_categorical_with_baseline_to_numeric <- function(col, data_col, col_info) {
#   # set what level is baseline
#   data_col <- factor(data_col)
#   data_col <- relevel(data_col, col_info)
#   # default is treatment contrast
#   data_col_conv <- stats::model.matrix(~data_col)   # stats::as.formula(paste0("~ 0 + ", col))
#   column_names <- colnames(data_col_conv)  # data_colblue, data_colgreen etc.
#   column_names <- str_remove(column_names, "data_col")
#   column_names <- paste0(col, "_", column_names)
#   colnames(data_col_conv) <- column_names
#   data_col_conv <- as.data.frame(data_col_conv)
#   data_col_conv[column_names] <- sapply(data_col_conv[column_names], as.numeric)
#   return(data_col_conv)
# }


# === convert_binary_to_numeric ===


#' #' @title 
#' #' 
#' #' @description 
#' #' 
#' #' @details
#' #' 
#' #' @param col character; 
#' #' @param data_col character vector; 
#' #' @param col character; 
#' #' @return data_col_conv tibble; 
#' convert_binary_to_numeric <- function(col, data_col) {
#'   unique_vals <- unique(data_col)
#'   data_col_conv <- data_col
#'   data_col_conv[data_col_conv == unique_vals[1]] = 0
#'   data_col_conv[data_col_conv == unique_vals[2]] = 1
#'   data_col_conv <- tibble::tibble(as.numeric(data_col_conv))
#'   colnames(data_col_conv) <- col
#'   return(data_col_conv)
#' }
