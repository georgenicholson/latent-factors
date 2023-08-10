
# Transform numeric columns with normal score transformation:
ns.transform <- function(y) {
  ns <- qqnorm(y, plot.it = FALSE)
  ns <- tbl_df(data.frame(x = ns$x, y = ns$y)) %>% 
    group_by(factor(y)) %>% 
    summarise(x = mean(x), y = mean(y))
  ns.transform.fun <- approxfun(ns$y, ns$x, rule = 2)
  yt <- ns.transform.fun(y)
  return(yt)
}


normal_quantile_map <- function(y) {
  y_quant <- rank(y, na.last = "keep") / (length(y))
  y_out <- qnorm(y_quant)
  return(y_out)
}

