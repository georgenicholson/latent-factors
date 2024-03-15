
#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}

control <- get_control_parameters_simulations(parameter_set = 1)


time_of_final_week_effect <- control$time_pt_eval
time_of_half_max_effect <- control$time_of_half_max_effect
size_of_final_week_effect_placebo <- -control$size_of_final_week_effect_placebo
size_of_final_week_effect_active_effect <- -(control$size_of_final_week_effect_active - control$size_of_final_week_effect_placebo)
xseq <- seq(0, time_of_final_week_effect, len = 100)
y1 <- emax_curve(t_eval = xseq, 
                    time_of_final_week_effect = time_of_final_week_effect,
                    time_of_half_max_effect = time_of_half_max_effect,
                    size_of_final_week_effect = -size_of_final_week_effect_placebo)
y2 <- emax_curve(t_eval = xseq, 
                    time_of_final_week_effect = time_of_final_week_effect,
                    time_of_half_max_effect = time_of_half_max_effect,
                    size_of_final_week_effect = -(size_of_final_week_effect_placebo + size_of_final_week_effect_active_effect))
dir.create("plots", showWarnings = FALSE)

cexax <- .65
pdf(file = control$emax_example_plot_file, 6, 4)
par(oma = c(1, 1, 1, 7))
matplot(xseq, cbind(y1, y2), ty = "l",
        xlab = "Week",
        ylab = "",
        lty = 1,
        las = 1,
        xaxs = "i",
        yaxs = "r",
        yaxt = "s")
axis(side = 4,
     at = -c(size_of_final_week_effect_placebo, size_of_final_week_effect_placebo + size_of_final_week_effect_active_effect),
     labels = c("Placebo, week 16", "Active, week 16"),
     las = 2,
     cex.axis = cexax)
axis(side = 1,
     at = time_of_half_max_effect,
     labels = expression(EC[50]),
     las = 1)
abline(v = time_of_half_max_effect,
       lty = 2)
par(xpd = NA)
arrows(x0 = time_of_final_week_effect * 1.01,
       x1 = time_of_final_week_effect * 1.01,
       y0 = -size_of_final_week_effect_placebo,
       y1 = -(size_of_final_week_effect_placebo + size_of_final_week_effect_active_effect),
       code = 3,
       length = .05, 
       angle = 20,
       col = 1
      )
par(xpd = FALSE)
mtext(side = 4,
      text = "ATE",
      at = -size_of_final_week_effect_placebo - .5 * size_of_final_week_effect_active_effect,
      line = .25,
      las = 1,
      cex = cexax)



lines(x = rep(time_of_final_week_effect * 1.01, 2), 
      y = -c(size_of_final_week_effect_placebo,
            size_of_final_week_effect_placebo + size_of_final_week_effect_active_effect),
      col = 2)
dev.off()


