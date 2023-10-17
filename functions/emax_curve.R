emax_curve <- function(t_eval,
                          time_of_final_week_effect = 16,
                          size_of_final_week_effect = 1,
                          time_of_half_max_effect = .5){
  y <- size_of_final_week_effect * (1 + time_of_half_max_effect / time_of_final_week_effect) / 
    (1 + time_of_half_max_effect / t_eval)
  return(y)
  
}

