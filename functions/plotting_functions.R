

#########################################################################
# Function to add letter to plot
add_letter <- function(letter = "a", x_fac = .05, y_fac = .05, 
                       cex = 1.5, ...) {
  par(xpd = NA)
  text(x = par("usr")[1] - diff(par("usr")[1:2]) * x_fac, 
       y = par("usr")[4] + diff(par("usr")[3:4]) * y_fac,
       label = paste0("(", letter, ")"), 
       cex = cex, 
       ...)
  par(xpd = F)
}

#########################################################################
# Function to plot scores
score_plot_types <- c("treatment_means", "contrasts_vs_placebo")
plot_scores <- function(control, 
                        analysis_results_list, 
                        analysis_plot,
                        fac_plot,
                        estimate_type = c("treatment_means", "contrasts_vs_placebo")[2],
                        add_legend = FALSE,
                        arms_in = NULL,
                        legend_y_mult = 1,
                        add_zero_horiz_line = TRUE,
                        ...) {
  if (is.null(arms_in)) {
    arms_in <- analysis_results_list[[analysis_plot[1]]][[control[[estimate_type]]$arm_include_label]]
  }
  m_pl <- analysis_results_list[[analysis_plot[1]]][[fac_plot]][[estimate_type]]$m[, arms_in, drop = F]
  l_pl <- analysis_results_list[[analysis_plot[1]]][[fac_plot]][[estimate_type]]$l[, arms_in, drop = F]
  u_pl <- analysis_results_list[[analysis_plot[1]]][[fac_plot]][[estimate_type]]$u[, arms_in, drop = F]
  t_unique <- as.numeric(rownames(m_pl))
  matplot(x = t_unique, 
          y = m_pl, 
          ylim = control[[estimate_type]]$ylim_curr,
          ty = "n", 
          xlab = "", 
          ylab = "", 
          las = 1,
          yaxs = "i",
          xaxs = "i", 
          cex.axis = control$cexax,
          ...)
  for (analysis_plot_curr in analysis_plot) {
    m_pl <- analysis_results_list[[analysis_plot_curr]][[fac_plot]][[estimate_type]]$m[, arms_in, drop = F]
    l_pl <- analysis_results_list[[analysis_plot_curr]][[fac_plot]][[estimate_type]]$l[, arms_in, drop = F]
    u_pl <- analysis_results_list[[analysis_plot_curr]][[fac_plot]][[estimate_type]]$u[, arms_in, drop = F]
    for (arm_curr in arms_in) {
      lines(x = t_unique, 
            y = m_pl[, arm_curr], 
            col = control$col_arm[arm_curr], 
            lwd = 2,
            lty = match(analysis_plot_curr, analysis_plot))
      polygon(x = c(t_unique, rev(t_unique)),
              y = c(l_pl[, arm_curr], rev(u_pl[, arm_curr])),
              col = control$col_arm[arm_curr],
              border = 1)
      lines(x = t_unique, y = m_pl[, arm_curr], col = 1, lwd = 2,
            lty = match(analysis_plot_curr, analysis_plot))
      lines(x = t_unique, y = m_pl[, arm_curr], col = control$col_arm_opaque[arm_curr], lwd = 2,
            lty = match(analysis_plot_curr, analysis_plot))
    }
  }
  if (add_zero_horiz_line) {
    abline(h = 0)
  }
  if (add_legend) {
    par(xpd = NA)
    legend_title <- switch(estimate_type,
                           treatment_means = "Treatment arm means",
                           contrasts_vs_placebo = "Average treatment effects")
    legend(x = max(t_unique) * 1.15, 
           # y = mean(par("usr")[3:4]), 
           y = par("usr")[3] + (par("usr")[4] - par("usr")[3]) * legend_y_mult,
           legend = paste0(arms_in, ifelse(estimate_type == "contrasts_vs_placebo", " - Placebo", "")), 
           col = 1,#control$col_arm[arms_in],
           pch = 22, 
           pt.bg = control$col_arm[arms_in], 
           cex = 1.5,
           pt.cex = 1.5,
           yjust = .5,
           title = legend_title,
           title.adj = 0.5
          )
    par(xpd = F)
  }
}

#########################################################################
# Function to add scale bar
add_scale_bar <- function (control,
                           color_vec, 
                           legend_range, 
                           y_at_fac = 0.05, 
                           y_width_fac = 0.025,
                           x_shrink_fac = .75) {
  y_lim <- par("usr")[3:4]
  x_lim <- par("usr")[1:2]
  x_lim_legend <- x_lim + .5 * c(1, -1) * diff(x_lim) * (1 - x_shrink_fac)
  
  y_wid <- diff(y_lim) * y_width_fac
  y_at <- y_lim[1] + diff(y_lim) * (1 + y_at_fac)
  n_col <- length(color_vec)
  mini_xlims <- c()
  mini_ylims <- y_at + c(0, y_wid)
  par(xpd = NA)
  for (j in 1:n_col) {
    mini_xlims[1] <- x_lim_legend[1] + ((j - 1) / n_col) * diff(x_lim_legend)
    mini_xlims[2] <- x_lim_legend[1] + ((j) / n_col) * diff(x_lim_legend)
    polygon(x = mini_xlims[c(1, 2, 2, 1, 1)], 
            y = mini_ylims[c(1, 1, 2, 2, 1)],
            col = color_vec[j], border = NA)
  }
  labs <- c(-1, 0, 1)#pretty(legend_range)
  axis(side = 3, 
       labels = labs, 
       at = x_lim_legend[1] + (labs - legend_range[1]) / diff(legend_range) * diff(x_lim_legend),
       pos = mini_ylims[2], 
       cex.axis = control$cexax)
  par(xpd = F)
}


#########################################################################
# Function to plot loadings



plot_loadings <- function(control, 
                          analysis_results_list, 
                          analysis_plot,
                          add_legend = TRUE, 
                          add_endpt_labels = TRUE,
                          y_at_fac = .025,
                          rotate = FALSE,
                          add_numbers = FALSE,
                          ...) {
  A_curr <- analysis_results_list[[analysis_plot]]$mfd_object$A
  if (rotate) {
    xpl <- 1:n_endpt
    ypl <- 1:control$n_fac
    zpl <- A_curr
  } else {
    xpl <- 1:control$n_fac
    ypl <- 1:n_endpt
    zpl <- t(A_curr)
  }
  image(xpl, ypl, zpl,
        col = control$col_heat, 
        zlim = c(-1, 1), 
        xaxt = "n", 
        yaxt = "n",
        xlab = "",
        ylab = "",
        axes = T,
        ...)
  if (add_numbers) {
    for (i in 1:nrow(zpl)) {
      for (j in 1:ncol(zpl)) {
        text(x = xpl[i], y = ypl[j], labels = round(zpl[i, j], 2))
      }
    }
  }
  if (add_endpt_labels) {
    axis(side = ifelse(rotate, 1, 2), 
         at = 1:n_endpt, 
         labels = endpts, 
         las = 2, 
         cex.axis = control$cexax)
  }
  axis(side = ifelse(rotate, 2, 1), 
       at = 1:control$n_fac, 
       labels = control$fac_labs, 
       las = 2, 
       cex.axis = control$cexax)
  if (add_legend) {
    add_scale_bar(control = control,
                  color_vec = control$col_heat, 
                  legend_range = c(-1, 1), 
                  y_at_fac)
  }
}

#########################################################################
# Add indication and study label
add_study_lab <- function(control, 
                          study_to_plot, 
                          eps_lab = .05) {
  indication <- ifelse(grepl("RA", study_to_plot), "RA", "PsA")
  studies <- paste0(paste(control$study_labs[control$do_list[[study_to_plot]]$stage2], collapse = " & "))
  study_lab <- paste0(indication, " ", studies)
  # mtext(side = 4, text = indication, cex = 1, las = 2, line = 2, 
        # at = par("usr")[3] + (par("usr")[4] - par("usr")[3]) * (0.5 + eps_lab))
  # mtext(side = 4, text = studies, cex = control$cexax, las = 2, line = 2, 
        # at = par("usr")[3] + (par("usr")[4] - par("usr")[3]) * (0.5 - eps_lab))
  mtext(side = 4, text = studies, cex = control$cexax, las = 2, line = 2, 
        at = par("usr")[3] + (par("usr")[4] - par("usr")[3]) * 0.5)
}

#########################################################################
# Add indication and study label
add_study_lab_top <- function(control = control,
                              study_to_plot) {
  indication <- ifelse(grepl("RA", study_to_plot), "RA", "PsA")
  studies <- paste0(paste(control$study_labs[control$do_list[[study_to_plot]]$stage2], collapse = " & "))
  study_lab <- paste0(indication, " ", studies)
  mtext(side = 3, text = indication, cex = 1, line = 1.6) 
  mtext(side = 3, text = studies, cex = .85, line = .4) 
}




