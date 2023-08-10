data_to_conv <- dplyr::as_tibble(data.frame(bin = factor(c("N", "Y", "N", "N", "Y"), levels=c("N", "Y")), 
                                            ord_cat = factor(c("no response", "high response", "light response", "light response", "no response"), levels = c("no response", "light response", "high response"), ordered=TRUE), 
                                            num=1:5, 
                                            char_col = c("a", "b", "c", "d", "e"), 
                                            unord_cat = factor(c("green", NA, "blue", "green", "red"), ordered=FALSE), 
                                            unord_cat_w_base = factor(c("no_event", "heartattack", "stroke", "incident", "stroke"), levels=c("no_event", "stroke", "incident", "heartattack"), ordered=FALSE )))


test_that("input: col is character", {
  col <- 'unord_cat'
  data_col <- data_to_conv %>% select(col) %>% pull()
  expect_error(convert_unordered_categorical_to_numeric(col=col, data_col=data_col), NA)
  
  expect_error(convert_unordered_categorical_to_numeric(col=1, data_col=data_col), "col is not a character.")
})

test_that("input: data_col is not an unordered factor.", {
  col <- 'unord_cat'
  data_col <- data_to_conv %>% select(col) %>% pull()
  expect_error(convert_unordered_categorical_to_numeric(col=col, data_col=data_col), NA)
  
  expect_error(convert_unordered_categorical_to_numeric(col=1, data_col=data_col), "col is not a character.")
})


# test_that("input: data_col is not an unordered factor.", {
#   col <- 'unord_cat'
#   data_col <- data_to_conv %>% select(col) %>% pull()
#   expect_error(convert_unordered_categorical_to_numeric(col=col, data_col=data_col), NA)
#   
#   col <- 'bin'
#   data_col <- data_to_conv %>% select(col) %>% pull()
#   expect_error(convert_unordered_categorical_to_numeric(col=col, data_col=data_col), NA)
#   
#   col <- 'unord_cat_w_base'
#   data_col <- data_to_conv %>% select(col) %>% pull()
#   expect_error(convert_unordered_categorical_to_numeric(col=col, data_col=data_col), NA)
#   
#   col <- 'ord_fac'
#   data_col <- data_to_conv %>% select(col) %>% pull()
#   expect_error(convert_unordered_categorical_to_numeric(col=col, data_col=data_col), "data_col is not an unordered factor.")
#   
#   col <- 'char_col'
#   data_col <- data_to_conv %>% select(col) %>% pull()
#   expect_error(convert_unordered_categorical_to_numeric(col=col, data_col=data_col), "data_col is not an unordered factor.")
# })

