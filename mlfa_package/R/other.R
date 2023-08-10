


#' @title
#' Calculate the latent trajectory weights.
#' 
#' @description
#' Estimates the latent trajectory weights transformed from matrix \code{A} that estimate scores matrix \code{Z} and its standard error matrix \code{S}.
#'
#' @param mfd_object mfd; the result object of \code{mfd}.

calculate_latent_weights <-function(mfd_object){
  
  A_tmp<-mfd_object$A
  colnames(A_tmp)<-paste0("Latent weights",1:ncol(A_tmp))
  # obtain the direct coefficient that calculate Z score from Y
  At_A_inv <- MASS::ginv(t(A_tmp) %*% A_tmp)  
  dimnames(At_A_inv) <- list(colnames(A_tmp), colnames(A_tmp))  # name dimensions of At_A_inv according to column names of A
  # because Zc <- At_A_inv %*% t(A_tmp) %*% t(Y_smp_non_NA), the direct coefficient from Y to Z is B=At_A_inv %*% t(A_tmp)
  B<-At_A_inv %*% t(A_tmp)
  
  return(B)
}





#' @title
#' Plot coefficients of \code{A}.
#' 
#' @description
#' Plot coefficients of the loadings matrix \code{A} (variables by (sparse) principal components).
#'
#' @param input_coefficients matrix; the result \code{A} matrix from \code{mfd} or the matrix from \code{calculate_latent_weights}, where the rows are the variables and the columns are different loadings.
#' @param save_path string; the path where to save the file to. Must end with \code{".pdf"}.
#' @param axis_y_text_size the size of the y tick labels (measurements).
#' @param add_dashline number; the shrinkage level to coefficient \code{A} matrix.
#' @param color_dashline;string; the color for the vertical shrinkage line.
#' @param by.facet the option to select whether by facet to display or not.
plot_coefficients <- function(input_coefficients, save_path=NULL, axis_y_text_size = NULL,by.facet="N",
                              add_dashline=NULL,color_dashline="black") {
  if (is.null(axis_y_text_size)) {
    n_vars <- length(rownames(input_coefficients))
    if (n_vars <= 20) {
      axis_y_text_size <- 6
    }
    else {
      axis_y_text_size <- 4
    }
  }
  
  data<-melt(data = input_coefficients, id.vars = c(rownames(input_coefficients)))
  colnames(data)<-c("Y","shape_num","Z")
  data$Factor<-as.factor(data$shape_num)
  
  if(by.facet=="N")  {
    p <- ggplot2::ggplot(data) +
      geom_point(aes(x=Z,  y = as.character(Y),fill=Factor,color=Factor),shape = data$shape_num) +
      theme_bw()+ geom_vline(xintercept = 0, linetype="dashed", 
                             color = "black", size=0.5)
    p <- p + ggplot2::ggtitle("Coefficient Plot") + 
      ggplot2::xlab("Coefficients") + 
      ggplot2::ylab("Measurements") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),strip.text.y = element_text(size = axis_y_text_size))  
    
    
    if(!is.null(add_dashline)){
      p<-p+geom_vline(xintercept = c(abs(add_dashline)*(-1),abs(add_dashline)), linetype="dotted", 
                      color = color_dashline, size=0.5)
    }
  } else {
    if(!is.null(add_dashline)){
      data<-data %>% mutate(x1=abs(add_dashline),x2=abs(add_dashline)*(-1),center=0)}
    
    p <- ggplot2::ggplot(data) +
      geom_point(aes(x=Z,  y = as.character(Y),fill=Factor,color=Factor),shape = data$shape_num) +
      theme_bw()+
      facet_grid(cols = vars(data$shape_num)) 
    
    p <- p + ggplot2::ggtitle("Coefficient Plot") + 
      ggplot2::xlab("Coefficients") + 
      ggplot2::ylab("Measurements") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),strip.text.y = element_text(size = axis_y_text_size))  
    
    
    if(!is.null(add_dashline)){
      p<-p +
        geom_vline(data =data,aes(xintercept=center), linetype="dotdash") +
        geom_vline(data =data,aes(xintercept=x1),color = color_dashline, linetype="dotted") +
        geom_vline(data =data,aes(xintercept=x2),color = color_dashline, linetype="dotted")
    }
  }
  
  if (!is.null(save_path)) {
    pdf(file.path(save_path,"coefficient_plot.pdf" ))
    p
    dev.off()
  }
  # methods::show(p)
  return(p)
}





#' @title
#' Longitudinal trajectories (wrapper function)
#' 
#' @description
#' Plot the longitudinal z trajectory of a \code{factor}, corresponding to an \code{mlfa} object of a certain \code{type}.
#'
#' @param mfd_object mfd; the result object of \code{mfd}.
#' @param save_path string; the path where to save the file to. Must end with \code{".pdf"}.
#' @param factor integer; the factor whose z profile shall be plotted.
#' @param time string; the column in the input dataframe to \code{mfd} that indicates the longitudinal dimension.
#' @param subject string; the column in the input dataframe to \code{mfd} that indicates the series.
#' @param type integer; type 1 plots each group of z trajectories side by side, type 2 plots z trajectories grouped by colour in one figure
#' @param colour_by string; the column in the input dataframe to \code{mfd} to colour the z trajectories by.
#'                                 By default \code{NULL}, meaning that no colouring is applied.
#' @param subject_data tibble = NULL; A tibble with a column \code{mfd$subject} and various other columns that contain subject specific information.
#'                             This information can be used e.g. for the colour_by column.
#'                             Default: \code{NULL} (not merging any data).
#' @param selected_subjects character vector; a vector containing the series of which the corresponding z profile will be plotted.
#'                        By default, all series are plotted (no selection).
#'                        Each entry in this vector appears in a column \code{subject} of the input data frame to the \code{mfd} call that indicates the series.
#'                        Multiple rows in the input data frame correspond to one series.
#' @param alpha decimal number; the transparency level of the lines
#' @param line.size number; the size of the line in the plot   
#' @param axis_text_size number; size of the text;                     
#' @return ggplot object.

plot_z_trajectory <- function(mfd_object, save_path, factor, time="AVISITN",subject="USUBJID", type = 1, colour_by = "TRT01P",subject_data = NULL, 
                              selected_subjects = NULL, alpha=0.5,line.size=1,axis_text_size=15) {
  data <- mfd_object$data_out
  # merge in subject_data
  data <- dplyr::inner_join(data, subject_data, by=c(mfd_object$subject))
  if (!is.null(selected_subjects)) {
    data <- data %>% dplyr::filter(!!dplyr::sym(subject) %in% selected_subjects)
  }
  
  if (type == 1) {
    z_chosen <- paste0("z", as.character(factor))
    grs <- (data %>% select(colour_by) %>% na.omit() %>% distinct %>% pull())
    colours <- rainbow(length(grs))
    data_temp <- data %>% select(time, z_chosen,subject,colour_by)
    p_sub <- ggplot2::ggplot(mapping=ggplot2::aes_string(x=time, y=z_chosen, group=subject,colour=colour_by))
    
    if (!is.null(colour_by)) {
      p_sub <- p_sub + ggplot2::geom_line(data=data_temp,  size=line.size, alpha=alpha) + 
        scale_color_manual(values=colours)+ facet_wrap(colour_by) #facet_grid(colour_by)
    }
    else {  ##ggplot2::aes_string(colour=colour_by),
      p_sub <- p_sub + ggplot2::geom_line(data=data_temp,  size=line.size, alpha=alpha)
    }
    p_sub <- p_sub + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p_sub <- p_sub + ggplot2::ggtitle(paste0(z_chosen, " trajectories")) + ggplot2::xlab(time) + ggplot2::ylab(z_chosen) +
      theme_bw()+ 
      theme(axis.text.x = element_text(size=axis_text_size), axis.text.y = element_text(size=axis_text_size)) +
      theme(axis.title.x = element_text(size=axis_text_size), axis.title.y = element_text(size=axis_text_size)) +
      theme(legend.position = "bottom", legend.text=element_text(size=axis_text_size)) +
      theme(legend.key.width = unit(2, "cm")) + 
      theme(strip.text.y = element_text(size = axis_text_size), strip.text.x = element_text(size = axis_text_size))
    
    p<-p_sub
  }
  else if (type == 2) {
    z_chosen <- paste0("z", as.character(factor))
    data_temp<-data %>% filter(!is.na(colour_by))
    p <- ggplot2::ggplot(mapping=ggplot2::aes_string(x=time, y=z_chosen, group=subject))
    if (!is.null(colour_by)) {
      p <- p + ggplot2::geom_line(data=data_temp, alpha=alpha, ggplot2::aes_string(colour=colour_by), size=2)
    }
    else {
      p <- p + ggplot2::geom_line(data=data_temp, alpha=alpha, size=line.size)
    }
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p <- p + ggplot2::ggtitle(paste0(z_chosen, " trajectories by ", subject)) + ggplot2::xlab(time) + ggplot2::ylab(z_chosen)  +
      theme_bw()+ 
      theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) +
      theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) +
      theme(legend.position = "bottom", legend.text=element_text(size=18)) +
      theme(legend.key.width = unit(2, "cm")) + 
      theme(strip.text.y = element_text(size = 18), strip.text.x = element_text(size = 18))
  }
  if (!is.null(save_path)) {
    ggplot2::ggsave(file=file.path(save_path,"z_by_patients_plot.pdf"), dpi=2000)
  }
  methods::show(p)
  return(p)
}







#' @title
#' Longitudinal trajectories (wrapper function)
#' 
#' @description
#' Plot the longitudinal z trajectory of a \code{factor}, corresponding to an \code{mlfa} object of a certain \code{type}.
#'
#' @param input data frame; the result object of \code{mfd}.
#' @param save_path string; the path where to save the file to. Must end with \code{".pdf"}.
#' @param time string; the column in the input dataframe to \code{mfd} that indicates the longitudinal dimension.
#' @param subject string; the column in the input dataframe to \code{mfd} that indicates the series.
#' @param groups string; the column in the input dataframe to \code{mfd} to colour the z trajectories by.
#'                                 By default \code{NULL}, meaning that no colouring is applied.
#' @param subject_data tibble = NULL; A tibble with a column \code{mfd$subject} and various other columns that contain subject specific information.
#'                             This information can be used e.g. for the groups column.
#'                             Default: \code{NULL} (not merging any data).
#' @param facet.var string; the column name in the subject_data data that indicates the facet variable.
#' @param selected_subjects character vector; a vector containing the series of which the corresponding z profile will be plotted.
#'                        By default, all series are plotted (no selection).
#'                        Each entry in this vector appears in a column \code{subject} of the input data frame to the \code{mfd} call that indicates the series.
#'                        Multiple rows in the input data frame correspond to one series.
#' @param labels vector; the label name for each lading.
#' @type type of plot; if type is "mean", the plot shows the mean scores; if the type is "CI", the plot shows the mean score with confidence interval.                        
#' @return ggplot object.
plot_z_by_group <- function(mfd_object, save_path,groups=NULL,time="AVISITN",subject="USUBJID",subject_data = NULL,type="CI",
                            selected_subjects = NULL,xlab="Time",ylab="Latent Scores",font.size=18,x.breaks=c(0,2,4,8,12,16),
                            facet.var=NULL,
                            ncol=3,ylim=c(0,3),color.theme="Paired",labels=c("Domain1","Domain2","Domain3")) {
  data <- mfd_object$data_out
  data <- dplyr::inner_join(data, subject_data, by=c(mfd_object$subject))
  if (!is.null(selected_subjects)) {
    data <- data %>% dplyr::filter(!!dplyr::sym(subject) %in% selected_subjects)
  }
  
  if (is.null(facet.var)) {
    data_tmp<-melt(data = data %>% select(time,subject,colnames(mfd_object$Z),groups)%>% as.data.frame(),id.vars = c(time,subject,groups), measure.vars =colnames(mfd_object$Z))
    colnames(data_tmp)<-c("time","subject","groups","variable","value")
    ##
    data_tmp$groups<-as.factor(data_tmp$groups)
    if (!is.null(labels)){
      data_tmp$variable<-factor(data_tmp$variable, levels=paste0("z", seq(1:length(unique(data_tmp$variable)))),labels = labels)
    }
    if(type=="mean") { 
      p<-ggplot2::ggplot(data_tmp, aes_string(x=data_tmp$time, color=data_tmp$groups, shape=data_tmp$groups, 
                                              group=data_tmp$groups, y=data_tmp$value)) +
        labs(x=xlab,y=ylab)+   
        theme_bw() + 
        theme(axis.text.x = element_text(size=font.size), axis.text.y = element_text(size=font.size)) +
        theme(axis.title.x = element_text(size=font.size), axis.title.y = element_text(size=font.size)) +
        theme(legend.position = "bottom", legend.text=element_text(size=font.size),legend.title = element_blank()) +
        theme(legend.key.width = unit(2, "cm")) + 
        stat_summary(fun.y=mean, geom="line", size=2) + 
        facet_wrap(.~data_tmp$variable,ncol=ncol) +
        
        theme(strip.text.y = element_text(size = font.size), strip.text.x = element_text(size = font.size))
      if (!is.null(color.theme)){
        p<-p+
          scale_color_brewer(palette=color.theme)
      }
    }else{
      
      tmp_mean<-data_tmp %>%
        group_by(time,groups,variable) %>%
        summarise_at(vars(value), funs(mean(., na.rm=TRUE)))  
      
      colnames(tmp_mean)<-c("time","groups","variable","mean")
      
      
      stderr <- function(x, na.rm=FALSE) {
        if (na.rm) 
          x <- na.omit(x)
        sqrt(var(x)/length(x))
      }
      
      tmp_se<-data_tmp %>%
        group_by(time,groups,variable) %>%
        summarise_at(vars(value), funs(stderr)) 
      
      colnames(tmp_se)<-c("time","groups","variable","se")
      
      tmp_se<-tmp_se %>% full_join(tmp_mean) %>%
        mutate(upper=mean+1.96*se)%>%
        mutate(lower=mean-1.96*se)
      
      
      p<-ggplot2::ggplot(data=tmp_se, aes(y=mean, x=time,  group=groups)) +
        geom_line(aes(colour=groups),size = 2) +
        geom_ribbon(aes(ymin=lower,ymax=upper,fill=groups), alpha=0.3, linetype=0) +
        xlab(xlab)+ylab(ylab)+
        scale_x_continuous(breaks=x.breaks, labels=x.breaks)+
        theme_bw() + 
        theme(axis.text.x = element_text(size=font.size), axis.text.y = element_text(size=font.size)) +
        theme(axis.title.x = element_text(size=font.size), axis.title.y = element_text(size=font.size)) +
        theme(legend.position = "bottom", legend.text=element_text(size=font.size),legend.title = element_blank()) +
        theme(legend.key.width = unit(2, "cm")) + 
        facet_wrap(.~variable,ncol=ncol) +
        theme(strip.text.y = element_text(size = font.size), strip.text.x = element_text(size = font.size))
      
      if (!is.null(color.theme)){
        p<-p+
          scale_fill_brewer(palette=color.theme) +
          scale_color_brewer(palette=color.theme)
      }
      
    }
    
  } else {
    
    data_tmp<-melt(data = data %>% select(time,subject,facet.var,colnames(mfd_object$Z),groups)%>% as.data.frame(),id.vars = c(time,subject,facet.var,groups), measure.vars =colnames(mfd_object$Z))
    colnames(data_tmp)<-c("time","subject","facet","groups","variable","value")
    ##
    data_tmp <- data_tmp %>% mutate(groups=as.factor(groups),facet=as.factor(facet))
    
    if (!is.null(labels)){
      data_tmp$variable<-factor(data_tmp$variable, levels=paste0("z", seq(1:length(unique(data_tmp$variable)))),labels = labels)
    }
    if(type=="mean") { 
      p<-
        ggplot2::ggplot(data_tmp, aes_string(x=data_tmp$time, color=data_tmp$groups, shape=data_tmp$groups, 
                                             group=data_tmp$groups, y=data_tmp$value)) +
        labs(x=xlab,y=ylab)+  
        facet_grid(facet~variable) +
        theme_bw() + 
        theme(axis.text.x = element_text(size=font.size), axis.text.y = element_text(size=font.size)) +
        theme(axis.title.x = element_text(size=font.size), axis.title.y = element_text(size=font.size)) +
        theme(legend.position = "bottom", legend.text=element_text(size=font.size),legend.title = element_blank()) +
        theme(legend.key.width = unit(2, "cm")) + 
        stat_summary(fun.y=mean, geom="line", size=2) + 
        theme(strip.text.y = element_text(size = font.size), strip.text.x = element_text(size = font.size))
      
      if (!is.null(color.theme)){
        p<-p+
          scale_color_brewer(palette=color.theme)
      }
    }else{
      
      tmp_mean<-data_tmp %>%
        group_by(time,groups,variable,facet) %>%
        summarise_at(vars(value), funs(mean(., na.rm=TRUE)))  
      
      colnames(tmp_mean)<-c("time","groups","variable","facet","mean")
      
      
      stderr <- function(x, na.rm=FALSE) {
        if (na.rm) 
          x <- na.omit(x)
        sqrt(var(x)/length(x))
      }
      
      tmp_se<-data_tmp %>%
        group_by(time,groups,variable,facet) %>%
        summarise_at(vars(value), funs(stderr)) 
      
      colnames(tmp_se)<-c("time","groups","variable","facet","se")
      
      tmp_se<-tmp_se %>% full_join(tmp_mean) %>%
        mutate(upper=mean+1.96*se)%>%
        mutate(lower=mean-1.96*se)
      
      
      p<-ggplot2::ggplot(data=tmp_se, aes(y=mean, x=time,  group=groups)) +
        geom_line(aes(colour=groups),size = 2) +
        geom_ribbon(aes(ymin=lower,ymax=upper,fill=groups), alpha=0.3, linetype=0) +
        xlab(xlab)+ylab(ylab)+
        scale_x_continuous(breaks=x.breaks, labels=x.breaks)+
        facet_grid(facet~variable)+
        theme_bw() + 
        theme(axis.text.x = element_text(size=font.size), axis.text.y = element_text(size=font.size)) +
        theme(axis.title.x = element_text(size=font.size), axis.title.y = element_text(size=font.size)) +
        theme(legend.position = "bottom", legend.text=element_text(size=font.size),legend.title = element_blank()) +
        theme(legend.key.width = unit(2, "cm")) + 
        theme(strip.text.y = element_text(size = font.size), strip.text.x = element_text(size = font.size))
      
      if (!is.null(color.theme)){
        p<-p+
          scale_fill_brewer(palette=color.theme) +
          scale_color_brewer(palette=color.theme)
      }
      
    }
    
    
  }
  
  if (!is.null(ylim)){
    p<-p+ylim(ylim)
  }
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(file=file.path(save_path,"z_by_groups_plot.pdf"), dpi=2000)
  }
  
  return(p)
}



#' @title
#' Plot heatmap of data availability of each measurement along the time.
#' 
#' @description
#' Plot heatmap of data availability of each measurement along the time by category.
#'
#' @param input_data data frame; the data frame which includes the columns of measurements and time for each subject information.
#' @param Time character; the column name of Time in the \code{input_data} dataset.
#' @param category character; the column name of factor variable that needs to facet in the \code{input_data} dataset.
#' @param id.var.names vector; the column names which are not included to visualize the data availability in the \code{input_data} dataset.
#' @param font.size.y number; the font size of y axis labels.
#' @param font.size.x number; the font size of x axis labels.
#' @param ncol number; the number of categories displayed in the figure column.
#' @param color.low character; the color displayed for the lowest frequency count in the figure.
#' @param color.high character; the color displayed for the highest frequency count in the figure.
plot_data_availability_heatmap <- function(input_data=d_PsA, Time="AVISITN",category="STUDYID",id.var.names=c("STUDYID","USUBJID","AVISITN","TRT01P"),
                                           font.size.y=9,font.size.x=9,ncol=3,color.low = "navy", color.high = "yellow") {
  
  input_data$Time<-input_data[,Time]
  input_data$category<-input_data[,category]
  id.var.names_v2<-unique(c(id.var.names,"Time","category"))
  
  d_melt <- melt(input_data,id.vars = id.var.names_v2) %>% group_by(category,variable,Time) %>% summarise(non_na_count = sum(!is.na(value))) %>%  dplyr::rename(PARAMCD=variable) %>% mutate(Time=as.factor(Time))
  
  dt<-spread(d_melt,"Time","non_na_count")
  dt[is.na(dt)] <- 0  ## put missing counts as 0
  dtm<- melt(dt %>% as.data.frame(),id.vars = c("category","PARAMCD"))%>% dplyr::rename(Time=variable)
  
  
  p<-ggplot(dtm, aes(x=Time, y=PARAMCD)) +
    facet_wrap( ~ category,ncol = ncol) +
    geom_tile(aes(fill = value)) +
    scale_fill_continuous(low = color.low, high = color.high)+
    xlab("Time")+ylab("Frequencies")+
    theme_bw() + 
    theme(axis.text.x = element_text(size=font.size.x), axis.text.y = element_text(size=font.size.y, colour= "black"),
          axis.title.x = element_text(size=font.size.x), axis.title.y = element_text(size=font.size.y,colour= "black"),
          legend.position = "bottom", legend.text=element_text(size=12),legend.key.width = unit(2, "cm"),
          strip.text.y = element_text(size = font.size.y),strip.text.x = element_text(size = font.size.x))  
  
  return(p)
}



##


#' @title
#' Data rank check for GAMM
#' 
#' @description
#' check the rank of the data and return data without problem
#'
#' @param data data.frame; data used to fit the gamm function \code{fit_gamm}.
#' @param strata character vector; main effects to be used in the gamm function \code{fit_gamm}.
#' @param p integer vector; which latent factors to be checked.
#' @return data.frame.

check_data_rank = function(data,strata, p = 1){
  
  for(y in paste0("z", 1:p)){
    #fit a lm model to check deficient rank
    formula_lm <- paste0(y, "~ ")
    
    for (i in 1:length(strata)) {
      formula_lm <- paste0(formula_lm, strata[i], " + ")
    }
    formula_lm = gsub("\\+ $","",formula_lm)
    formula_lm <- stats::as.formula(  formula_lm)
    gm <- lm(formula_lm, data=data)
    # Identify which columns constitute a valid design matrix:
    valid.cols <- names(coef(gm)[!is.na(coef(gm))])[-1]
    
    data = as.data.frame(data)
    for( cc in strata[order(nchar(strata), decreasing = T)]){
      valid.cols <- gsub(cc, "", valid.cols)
      ref_lev = levels(data[,cc])[1]
      data[,cc][!(as.character(data[,cc]) %in% valid.cols)] <- ref_lev
      data[,cc] = droplevels(data[,cc])
    }
  }
  
  
  return(data)
}








#' @title
#' Multivariate Longitudinal Factor Analysis using GAMM (internal invoked).
#' 
#' @description 
#' Modified version of \code{mlfa} function. Fit \code{gamm} models marginally for each factor (stage 3) that allow to predict longitudinal latent trajectories.
#'
#' @param data data.frame; merged from A fitted mfd object (result object of stage 1 and 2) and subject_data.
#' @param df int; The number of degrees of freedom of the fitted splines. The higher \code{df}, the wigglier the fitted spline.
#' @param factor int; which latent factor to fit
#' @param strata character vector; The strata variable(s) of the \code{gamm} model.
#' @param group string; The group variable of the \code{gamm} model. Usually \code{mfd$subject}.
#' @param AR1_group string; The AR1 auto-cor variable of the \code{gamm} model. Usually \code{mfd$subject}.
#' @return \code{gamm_object} object; an object of S3 class \code{mlfa}.
#'               # Primary (necessary inputs)
gamm_fit_internal <- function(data, df, factor = 1, time, strata, group, AR1_group= 1, numerical_fix = NULL) {
  
  formula <- paste0(paste0("z", as.character(factor), " ~ ", "s(", time, ", k = ", as.character(df), ")"))
  for (i in 1:length(strata)) {
    formula <- paste0(formula, " + ", strata[i], " + ", "s(", time, ", k = ", as.character(df), ", by = factor(", strata[i], ", ordered = TRUE))")
  }
  if(!is.null(numerical_fix)){
    for (i in 1:length(numerical_fix)) {
      formula <- paste0(formula, " + ", numerical_fix[i])
    }
  }
  formula <- stats::as.formula(formula)
  ##random intercept
  random_list <- list()
  for(gg in group){
    random_list[[gg]] <- stats::as.formula(~1)
  }
  
  data = data %>% filter_at(vars(paste0("s", factor)), ~ . != 0)
  sv = data %>% pull(paste0("s", as.character(factor)))
  data$var = sv^2
  
  gamm_object <- mgcv::gamm(
    formula = formula,
    random = random_list,
    weights =varFixed(~var),
    correlation=corCAR1(form=as.formula(paste0("~", time, "|", AR1_group))),
    data = data,
    niterPQL=50
  )
  
  # stored for plotting
  gamm_object$time <- time
  gamm_object$strata <- strata
  gamm_object$group <- group
  gamm_object$factor <- factor
  
  return(gamm_object)
}





#' @title
#' Multivariate Longitudinal Factor Analysis using GAMM (Wrapper).
#' 
#' @description 
#' Modified version of \code{mlfa} function. Fit \code{gamm} models marginally for each factor (stage 3) that allow to predict longitudinal latent trajectories.
#'
#' @param mfd \code{mfd} object; A fitted mfd object (result object of stage 1 and 2).
#' @param subject_data tibble; A tibble with a column \code{mfd$subject} and various other columns that contain subject specific information.
#' @param strata character vector; The strata variable(s) of the \code{gamm} model.
#' @param trt_levels character vector; level of treatment factor
#' @param fit_factors int vector; which latent factor to be fitted
#' @return list; contains a list of fitted gamm (\code{gamm_object} object; an object of S3 class \code{mlfa}) and the rank-checked data.

fit_gamm = function(mfd, subject_data,
                    strata, trt_levels, fit_factors = NULL){
  ## # of z factors:
  p = mfd$fac_dim
  ## name of time variables
  time = mfd$time
  ## create auto correlation within the same subjects, name of subject id:
  subject=mfd$subject
  AR1_group = subject
  
  ## random effect of subject, name of subject id:
  group = subject
  
  
  
  if(is.null(fit_factors)){
    fit_factors = 1:p
  }
  
  #merge data
  subject_data = subject_data %>% mutate( TRT = factor(.data[[TRT]], levels = trt_levels))
  data <- dplyr::inner_join(mfd$data_out, subject_data, by=c(mfd$subject)) %>%
    mutate_at(vars(unique(c(strata, AR1_group, group)) ), as.factor )%>%
    filter(!is.na(z1)) %>%
    filter_at(vars(paste0("s", fit_factors)), ~ . != 0) %>% # remove invalid values (s = 0)
    mutate_if(is.factor, droplevels)
  
  ## # of time points -1
  df = length(unique(data$AVISITN))-1
  
  data_r = check_data_rank(data, strata = strata, p = p)
  
  model_list = list()
  
  
  for(i in fit_factors){
    
    model_list[[i]] = gamm_fit_internal(data = data_r, df = df, factor = i,
                                        time = time, strata = strata, group = group,
                                        AR1_group = AR1_group)
  }
  
  
  return(list(model_list = model_list, used_data = data_r))
  
}






#' @title
#' Treatment comparison of trajectory (wrapper function)
#' 
#' @description
#' Plot the trt comparison of the longitudinal z trajectory of a \code{factor} and the 95% CI, corresponding to an \code{mlfa} object of a certain \code{type}.
#'
#' @param model_list list of fitted \code{gamm}; the result object of \code{fit_gamm}.
#' @param data tibble; data used to get the gamm model; the result object of \code{fit_gamm}.
#' @param save_path string; the path where to save the file to. Must end with \code{".pdf"}.
#' @param time string; the column in the input dataframe to \code{mfd} that indicates the longitudinal dimension.
#' @param subject string; the column in the input dataframe to \code{mfd} that indicates the series.
#' @param STUDYID string; the column in the input dataframe to \code{mfd} that indicates the study id.
#' @param TRT string; the column in the input dataframe to \code{mfd} that indicates the treatment arm.
#' @param contrast_name string; name of the comparison contrast.
#' @param contrast int vector; contrast vector with the corresponding trt arm name (need to be consistent with the trt levels in the input data).
#' @param pp int; how many replicates are simulated to get the Bayesian credible interval
#' @param labels vector; the label name for each lading.
#' @return list; contain a ggplot object and a summary table.

plot_trt_comparison <- function(model_list, data,
                                time="AVISITN",subject="USUBJID", STUDYID = "STUDYID", TRT = "newTRT",
                                contrast_name = NULL, contrast, pp = 1000,
                                font.size=18,x.breaks=c(0,2,4,8,12,16),
                                xlab = "Time",
                                ncol=3, color.theme="Paired",labels= NULL,
                                save_path = NULL) {
  
  strata = model_list[[1]]$strata
  time = model_list[[1]]$time
  # check if there are studies contain two trt
  contrast_levels = names(contrast)[contrast != 0]
  level_sutdy = data %>% select(TRT, STUDYID) %>% table()
  used_studies = which(colSums(level_sutdy * (rownames(level_sutdy) %in% contrast_levels) != 0) == length(contrast_levels))
  
  if(length(used_studies) == 0){
    cat("There is NO study used the two levels!")
    return(NULL)
  }
  used_studies = names(used_studies)
  
  pd = data %>% select(time, strata) %>% distinct() %>% as.data.frame() %>% mutate_if(is.factor, droplevels)
  #pd = pd[order(pd[,TRT], pd[,time]),]  
  
  p = length(model_list)
  
  trajectories = list()
  for(i in 1:p){
    ##design matrix for each z score
    Xp <- predict(model_list[[i]]$gam, pd, type="lpmatrix",  newdata.guaranteed=TRUE)
    br <- MASS::mvrnorm(pp, coef(model_list[[i]]$gam), vcov(model_list[[i]]$gam, unconditional = T))
    trajectories[[i]] <-Xp%*%t(br) %>% as.data.frame()
    names(trajectories[[i]]) <- paste0(paste0("z",i,"_"),1:pp)
  }
  
  ##summarize and plot
  sample_size = data %>% dplyr::rename(time = time, TRT = TRT, STUDYID = STUDYID) %>%
    group_by(TRT, time, STUDYID) %>% tally()
  
  # combine data
  pd = pd %>% dplyr::rename(time = time, TRT = TRT, STUDYID = STUDYID)
  
  pd_full <- bind_cols(pd, trajectories)
  gg.dat <- pd_full %>% gather(key="simulation", value="value", -c(1:ncol(pd)))
  gg.dat2 = gg.dat %>% filter(STUDYID %in% used_studies) %>%
    mutate(constrast_sep_value = (value*contrast[as.character(TRT)])) %>%
    group_by(simulation, time, STUDYID)  %>%
    left_join(sample_size) %>%
    summarise(delta = sum(constrast_sep_value), n_size = sum( 1*(contrast[as.character(TRT)]!= 0)/n ))%>%
    dplyr::mutate(variable=substr(simulation, 1, regexpr("_", simulation)-1))
  
  
  if (is.null(labels)){
    labels = paste0("Domain ", 1:length(unique(gg.dat2$variable)))
  }
  gg.dat2$variable<-factor(gg.dat2$variable, levels=paste0("z", seq(1:length(unique(gg.dat2$variable)))),labels = labels)
  
  
  mean.gg.dat = gg.dat2 %>%
    ungroup() %>%
    group_by(variable, time, STUDYID) %>%
    summarise(mean=mean(delta), var = var(delta), sd = sd(delta), se = sd(delta)*sqrt(n_size)[1],
              lwr=quantile(delta, .025), upr=quantile(delta, 0.975)) %>%
    mutate(contrast = contrast_name)
  
  if(length(used_studies) > 1){
    mean.gg.dat = mean.gg.dat %>%
      ungroup() %>%
      group_by(mean.gg.dat, time) %>%
      mutate(wi = (1/var)/sum(1/var), pooled_est = sum(mean*wi), pooled_var = (1/sum(1/var)),
             pooled_lwr = pooled_est+qnorm(0.025)*sqrt(pooled_var), pooled_upr = pooled_est+qnorm(0.975)*sqrt(pooled_var))
    
    
    p_CI <-ggplot2::ggplot(data=mean.gg.dat, aes(y= pooled_est, x=time)) +
      geom_line(size = 2) +
      geom_ribbon(aes(ymin= pooled_lwr,ymax=  pooled_upr), alpha=0.3, linetype=0)
    
    
  }else{
    p_CI <-ggplot2::ggplot(data=mean.gg.dat, aes(y= mean, x=time)) +
      geom_line(size = 2) +
      geom_ribbon(aes(ymin = lwr,ymax = upr), alpha=0.3, linetype=0)
  }
  
  p_CI <-p_CI +
    geom_hline(yintercept = 0, col = "red", size = 1, linetype = "dashed") +
    xlab(xlab)+ylab(contrast_name)+
    scale_x_continuous(breaks=x.breaks, labels=x.breaks)+
    theme_bw() + 
    theme(axis.text.x = element_text(size=font.size), axis.text.y = element_text(size=font.size)) +
    theme(axis.title.x = element_text(size=font.size), axis.title.y = element_text(size=font.size)) +
    theme(legend.position = "bottom", legend.text=element_text(size=font.size),legend.title = element_blank()) +
    theme(legend.key.width = unit(2, "cm")) + 
    facet_wrap(.~variable,ncol=ncol) +
    theme(strip.text.y = element_text(size = font.size), strip.text.x = element_text(size = font.size))
  
  
  if(!is.null(save_path)){
    ggsave(p_CI, file=file.path(save_path), dpi=2000)
  }
  
  return(list(p_CI = p_CI, summary_table = mean.gg.dat))
  
  
}
