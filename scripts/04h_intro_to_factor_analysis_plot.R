source("scripts/04_plotting_intro.R")


study_plot <- "F2342"
subj_unique <- unique(d_time$USUBJID[d_time$study == study_plot])
n_miss <- c()
for (subj_plot in subj_unique) {
  d_time_sub <- d_time[which(d_time$USUBJID == subj_plot & d_time$AVISITN %in% c(0, 2, 4, 8, 12, 16)), ]
  ymat <- d_time_sub[, nms_params_v2$PARAMCD]
  n_miss[subj_plot] <- sum(is.na(ymat))
}



table(n_miss)
subj_plot <- subj_unique[which(n_miss == 8)[1]]#"CAIN457F2206_0153_00002"
week_plot <- c(0, 2, 4, 8, 12, 16)
d_time_sub <- d_time[which(d_time$USUBJID == subj_plot & d_time$AVISITN %in% week_plot), ]
study <- d_time_sub$study[1]
indication <- control$study_indic[grep(study, names(control$study_indic))]
mfd_obj <- analysis_results_list[[paste0(indication, "_", study_plot)]]$mfd_object
Z <- mfd_obj$Z[which(mfd_obj$data_out$USUBJID == subj_plot), ]
S <- mfd_obj$S[which(mfd_obj$data_out$USUBJID == subj_plot), ]

ymat <- d_time_sub[, nms_params_v2$PARAMCD]
rownames(ymat) <- d_time_sub$AVISITN
n_meas <- nrow(nms_params_v2)




graphics.off()
if (export_plots){
  pdf(file = paste0("figures/intro_to_factor_models.pdf"), 12, 8)
}

upper_band_wid <- c(.1)
right_band_wid <- .25
left_cen <- .12
upper_levels <- c(.75, .875)
right_cen <- .75#1 - right_band_wid - left_cen
arrow_col <- grey(.7)
arrow_length <- .25
arrow_width <- arrow_length * .75

##########################
# Add raw plots
upper_centres_x <- seq(left_cen, right_cen, length.out = n_meas)
upper_centres_y <- rep(upper_levels, length.out = n_meas)
names(upper_centres_y) <- names(upper_centres_x) <- nms_params_v2$PARAMCD
load_plot_hei <- .3
load_plot_wid <- .6
load_centre_x <- mean(upper_centres_x)
load_centre_y <- upper_levels[1] / 2
upper_plot_wid <- diff(base::range(upper_centres_x)) / n_meas 
upper_plot_hei <- upper_plot_wid * 1.15
par(mar = c(0, 0, 0, 0) * 2)
coords <- list()
plot(0, axes = F, ty = "n")
par(xpd = NA)
meas_codes <- nms_params_v2$PARAMCD
for (meas in meas_codes) {
  coords[[meas]] <- c(upper_centres_x[meas] + c(-1, 1) * upper_plot_wid / 2, 
                      upper_centres_y[meas] + c(-1, 1) * upper_plot_hei / 2)
                      
  par(fig = coords[[meas]], new = T)
  ypl <- ymat[as.character(week_plot), meas]
  par(mar = c(0, .5, 0, 0))
  cexax <- .75
  plot(week_plot, ypl, 
       ty = "b", 
       axes = T,
       xaxt = "n",
       yaxt = "s",
       ylab = "",
       xlab = "",
       bty = "l",
       cex = .6,
       pch = 19,
       las = 1,
       cex.axis = cexax)
  axis(side = 1, 
       at = c(0, 15), 
       labels = NA,
       cex.axis = cexax,
       tcl = -.25
       )
  mtext(side = 1, 
        at = c(0, 15), 
        text = c(0, 15),
        cex = cexax)
  mtext(side = 1, 
        line = .05, 
        text = "Week", 
        cex = cexax)
}
text(x = grconvertX(mean(upper_centres_x), from = "ndc"),
     y = grconvertY(upper_levels[2] + (1 - upper_levels[2]) * .65, from = "ndc"),
     labels = "Inputs: raw data from a single patient, Y")


# text_transform <- "Transform and take weighted average to generate factor scores (loading weights below)"
# llegend(x = grconvertX(left_cen / 2, from = "ndc"),
#      y = grconvertY(upper_levels[1] - (1 - upper_levels[1]) * .6, from = "ndc"),
#      labels = "Weighted average")

##########################
# Add scores
par(mar = c(0, 0, 0, 0))
fac_coords <- list()
score_y_offset <- .21
score_plot_hei <- .14
score_plot_wid <- .1
for (fac_num in 1:3) {
  fac_coords[[fac_num]] <- c(1 - right_band_wid * .4 + c(-1, 1) * score_plot_wid / 2,
                             load_centre_y + c(-1, 0, 1)[fac_num] * score_y_offset + c(-1, 1) * score_plot_hei / 2)
  par(fig = fac_coords[[fac_num]], new = T)
  
  plot(week_plot, Z[, fac_num], ty = "n", 
       ylim = base::range(c(c(Z + 2 * S), c(Z - 2 * S)), na.rm = T),
       xaxt = "n",
       yaxt = "n",
       las = 1,
       ylab = "",
       xlab = "")
  lines(week_plot, Z[, fac_num], lwd = 2)
  lines(week_plot, Z[, fac_num] + 2 * S[, fac_num], lwd = 1)
  lines(week_plot, Z[, fac_num] - 2 * S[, fac_num], lwd = 1)
  axis(side = 1, 
       at = c(0, 15), 
       labels = NA,
       cex.axis = cexax,
       tcl = -.25)
  mtext(side = 1, 
       at = c(0, 15), 
       text = c(0, 15),
       cex = cexax)
  axis(side = 4, 
       cex.axis = cexax,
       las = 2)
  mtext(side = 3, 
        text = paste0("Factor ", fac_num), 
        line = .15,
        cex = .9)
  mtext(side = 1, 
        line = .05, 
        text = "Week", 
        cex = cexax)
}
mtext(side = 3, text = "Outputs: scores, Z", line = 2)

##########################
# Add loadings heatmap
coords$loadings <- c(coords[[meas_codes[1]]][1], coords[[meas_codes[length(meas_codes)]]][2],
                     load_centre_y + c(-1, 1) * load_plot_hei / 2)
coords$loadings
par(fig = coords$loadings, new = T)

par(mar = c(0, 0, 0, 0))
plot_loadings(control, 
              analysis_results_list, 
              analysis_plot = "PsA_F2312_cross_indication_loadings",
              add_legend = F, 
              rotate = TRUE, 
              add_numbers = T)

for (meas in meas_codes) {
  xends <- grconvertX(rep(mean(coords[[meas]][1:2]), 2), from = "ndc", to = "user")
  yends <- grconvertY(c(coords[[meas]][3] - .04, coords$loadings[4]), from = "ndc", to = "user")
  shape::Arrows(x0 = xends[1], x1 = xends[2],
                y0 = yends[1], y1 = yends[2],
                lwd = 1, 
                arr.type = "triangle",
                col = arrow_col,
                arr.col = arrow_col,
                arr.adj = 1,
                arr.length = arrow_length,
                arr.width = arrow_width)
}

for (fac_num in 1:3) {
  yends <- c(fac_num, grconvertY(mean(fac_coords[[fac_num]][3:4]), from = "ndc", to = "user"))
  xends <- c(n_meas + .5, grconvertX(fac_coords[[fac_num]][1], from = "ndc", to = "user"))
  shape::Arrows(x0 = xends[1], x1 = xends[2],
                y0 = yends[1], y1 = yends[2],
                lwd = 1, 
                arr.type = "triangle",
                col = arrow_col,
                arr.col = arrow_col,
                arr.adj = 1,
                arr.length = arrow_length,
                arr.width = arrow_width)
  
}

text_transform <- "Take weighted average of transformed data using loadings, A, below"

legend(x = grconvertX(mean(upper_centres_x), from = "ndc"),
       y = grconvertY(upper_levels[1] - (1 - upper_levels[1]) * .6, from = "ndc"),
       legend = text_transform, 
       bg = "white",
       xjust = .5,
       yjust = .5,
       box.col = "white")

if (export_plots){
  dev.off()
}




