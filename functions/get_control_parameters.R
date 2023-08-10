#########################################################################
# Function to get plotting parameters
get_control_parameters <- function() {
  arm_all <- c("Placebo", "75mg", "150mg", "300mg", "Abatacept")
  n_fac <- 3
  
  control <- list(
    n_fac = n_fac,
    estimate_types = c("treatment_means", "contrasts_vs_placebo"),
    treatment_means = list(title_text = "Treatment mean trajectories",
                           arm_include_label = "arm_all",
                           ylim_curr = c(-1, 1) * 1.25),
    contrasts_vs_placebo = list(title_text = "Active treatment arms relative to placebo",
                                arm_include_label = "arm_no_plac",
                                ylim_curr = c(-1, .25)),
    col_heat = c(rgb(0, 0, 1, alpha = seq(1, 0, len = 500)), rgb(1, 0, 0, alpha = seq(0, 1, len = 500))),
    col_arm = viridis::viridis(n = length(arm_all), alpha = .5), 
    col_arm_opaque = viridis::viridis(n = length(arm_all), alpha = .5),
    arm_all = c("Placebo", "75mg", "150mg", "300mg", "Abatacept"),
    cexax = 1.3,
    fac_labs = paste0("Factor ", 1:n_fac),
    single_studies = c("PsA_F2312", "PsA_F2342", "RA_F2302", "RA_F2309"),
    study_indic = c("PsA_F2312" = "PsA", "PsA_F2342" = "PsA", "RA_F2302" = "RA", "RA_F2309" = "RA"),
    study_labs = c('F2312' = 'FUTURE-2',
                   'F2309' = 'NURTURE-1',
                   'F2302' = 'REASSURE',
                   'F2342' = 'FUTURE-5'),
    do_list = list('PsA_F2312' = 
                     list(stage1 = 'F2312', 
                          stage2 = 'F2312'),
                   'PsA_F2342' = 
                     list(stage1 = 'F2342', 
                          stage2 = 'F2342'),
                   'RA_F2302' = 
                     list(stage1 = 'F2302', 
                          stage2 = 'F2302'),
                   'RA_F2309' = 
                     list(stage1 = 'F2309', 
                          stage2 = 'F2309'),
                   'PsA_meta' = 
                     list(stage1 = c('F2312', 'F2342'), 
                          stage2 = c('F2312', 'F2342')),
                   'RA_meta' = 
                     list(stage1 = c('F2302', 'F2309'), 
                          stage2 = c('F2302', 'F2309')),
                   'PsA_F2312_cross_indication_loadings' = 
                     list(stage1 = c('F2312', 'F2342', 'F2302', 'F2309'), 
                          stage2 = 'F2312'),
                   'PsA_F2342_cross_indication_loadings' = 
                     list(stage1 = c('F2312', 'F2342', 'F2302', 'F2309'), 
                          stage2 = 'F2342'),
                   'RA_F2302_cross_indication_loadings' = 
                     list(stage1 = c('F2312', 'F2342', 'F2302', 'F2309'), 
                          stage2 = 'F2302'),
                   'RA_F2309_cross_indication_loadings' = 
                     list(stage1 = c('F2312', 'F2342', 'F2302', 'F2309'), 
                          stage2 = 'F2309'),
                   'PsA_meta_cross_indication_loadings' = 
                     list(stage1 = c('F2312', 'F2342', 'F2302', 'F2309'), 
                          stage2 = c('F2312', 'F2342')),
                   'RA_meta_cross_indication_loadings' = 
                     list(stage1 = c('F2312', 'F2342', 'F2302', 'F2309'), 
                          stage2 = c('F2302', 'F2309'))
    )
  )
  names(control$col_arm)<-names(control$col_arm_opaque) <- arm_all
  return(control)
}

