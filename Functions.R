# Functions cluster / IPD paper
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(magrittr)
library(table1)

library(splines)
library(locfit)
library(gam)

library(metamisc)
library(ggplot2)
library(latex2exp)


#' get_calmod: function to get the calibration information at cluster-specfic data using the loess, locfit and natural splines methods
#'
#' @param data # evaluated data (y and pi taken from a prediciton model)
#' @param new_df # new data with the pi range of expected risks
#' @param df  degress of freedom of the natural splines function
#' @param alpha confidence interval
#'
#' @return datset with expected observed probability pr along with CI for each method.
#'
get_calmod <- function(data, new_df, df = 3, alpha = 0.05){
  
  out <- NULL
  for (i in unique(data$model)){
    data_i <- data|>filter(model==i)
    new_df_i <- new_df|>filter(model==i)
    
    cal_df <- data.frame(y  = as.numeric(data_i$y==1), lpi = data_i$lpi)
    # loess
    cal_mod <- loess(y ~ lpi, degree = 2, data = cal_df,control=loess.control(surface="interpolate"))
    cal_pred <- predict(cal_mod,newdata=new_df_i,se=TRUE)
    pr <- cal_pred$fit
    pr[pr <= 0] <- 0
    pr[pr >= 1] <- 1
    lpr <- qlogis(pr)
    lpr[is.infinite(lpr)|is.na(lpr)] <- NA
    lpr.se <- cal_pred$se.fit/(pr*(1-pr))
    lpr.se[is.infinite(lpr.se)] <- NA
    new_df_i$lpr_loess <- lpr
    new_df_i$lprse_loess <- lpr.se
    
    # gam
    formula <- paste0("y ~ ns(lpi,",df,")")
    cal_mod <-  suppressWarnings(try( gam::gam(as.formula(formula), data = cal_df, family = binomial("logit"))))
    cal_pred2 <- predict(cal_mod,newdata = new_df_i, type = "link", se = TRUE)
    cal_pred2$fit[new_df_i$lpi < min(cal_df$lpi)|new_df_i$lpi > max(cal_df$lpi)] <- NA
    cal_pred2$se.fit[new_df_i$lpi < min(cal_df$lpi)|new_df_i$lpi > max(cal_df$lpi)] <- NA
    new_df_i$lpr_splines <- cal_pred2$fit
    new_df_i$lprse_splines <- cal_pred2$se.fit
    
    
    # locfit
    cal_mod <- locfit(y ~ lpi, family="binomial",data = cal_df, link="logit")
    cal_pred3 <- predict(object=  cal_mod, newdata = new_df_i, tr = function(x){x}, se=TRUE, band="local")
    cal_pred3$fit[new_df_i$lpi < min(cal_df$lpi)|new_df_i$lpi > max(cal_df$lpi)] <- NA
    cal_pred3$se.fit[new_df_i$lpi < min(cal_df$lpi)|new_df_i$lpi > max(cal_df$lpi)] <- NA
    new_df_i$lpr_locfit <- cal_pred3$fit
    new_df_i$lprse_locfit <- cal_pred3$se.fit
    
    new_df_i <- new_df_i |> 
      gather(v, value, lpr_loess:lprse_locfit) |>
      separate(v, c("var", "method")) |> 
      pivot_wider(names_from = var, values_from = value)|>
      mutate(pr = plogis(lpr),
             lower = plogis(lpr+qnorm(alpha/2)*lprse),
             upper = plogis(lpr+qnorm(1-alpha/2)*lprse),
             model = i)
    
    out <- rbind(out, new_df_i)
  }
  
  return(out)
}


#' get_dens: get the density of pi
#'
#' @param data # evaluated data (y and pi taken from a prediciton model)
#' @param new_df # new data with the pi range of expected risks
#'
#' @return datatable with densities of pi used in calibration plot
#'

get_dens <- function(data, new_df){
  out <- NULL
  
  for (i in unique(data$model)){
    data_i <- data|>filter(model==i)
    new_df_i <- new_df|>filter(model==i)
    data_i$pi_bin <- cut(data_i$pi, breaks=c(0,new_df_i$pi))
    get_midpoint <- function(cut_label) {
      mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
    }
    
    data_i$mid_bin <- sapply(data_i$pi_bin, get_midpoint)
    
    databins_i <- data.frame(table(data_i$y,data_i$mid_bin))|>
      mutate(Freq = 0.1*Freq/max(Freq))|>
      rename_with(~ c("y","pi","height"))|>
      mutate(pi = as.numeric(as.character(pi)),
             height = ifelse(y==0,-height-0.15,height-0.15))|>
      pivot_wider(names_from = y, values_from = height)|>
      rename_with(~ c("pi","min","max"))|>
      mutate(model = i)
    out <- rbind(out, databins_i)
  }
  return(out)
}

cal_plot_ind <- function(plot_data,pi_dens,alpha_i = 0.1, alpha_m = 0.3){
  plot_data$model<-paste0("Model ",plot_data$model)
  pi_dens$model<-paste0("Model ",pi_dens$model)
  meth_u<-unique(plot_data$method)
  alpha_val <- rep(alpha_i,length(meth_u))
  names(alpha_val)<-meth_u
  names_alpha<-names(alpha_val)
  alpha_val[startsWith(names_alpha, 'marg')] <- alpha_m
  
  p0<- ggplot(plot_data,aes(x=pi,y=pr))+
    geom_abline(intercept = 0, slope = 1,linetype=2) +
    geom_line(linewidth = 0.5,aes(group=cluster,color=cluster))+
    geom_ribbon(aes(ymin=lower, ymax=upper,fill=cluster,group=cluster),alpha=0.2,colour = NA) + 
    labs( x = TeX(r"( Estimated probability, $\pi_{k}$ )"),
          y = TeX(r"( Observed probability, $\rho_{k}$ )"),
          fill="Validation\ndataset",
          color="Validation\ndataset") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 12),legend.position = "right")+
    scale_y_continuous(limits = c(-0.25,1))+
    scale_alpha_manual(values=alpha_val)+
    geom_segment(aes(x=0,y=-0.15,xend=1,yend=-0.15),color="black",linewidth = 0.1)+
    geom_segment(data=pi_dens,aes(x=pi,y=min,xend=pi,yend=max),color="black",linewidth = 1)+
    annotate("text", x=0, y=-0.1, label= "1",size=2) + 
    annotate("text", x=0, y=-0.2, label= "0",size=2) +
    facet_grid(method~ model)+theme_classic()
  
  return(p0)
  
}

#' cal_plot plot calibration plot
#'
#' @param plot_data data with plot information
#' @param pi_dens  data with density information
#' @param cluster  logical displayed by cluster?
#'
#' @return Calibration plot
#'

cal_plot <- function(plot_data,pi_dens,cluster = FALSE, alpha_i = 0.1, alpha_m = 0.3){
  plot_data$model<-paste0("Model ",plot_data$model)
  pi_dens$model<-paste0("Model ",pi_dens$model)
  meth_u<-unique(plot_data$method)
  alpha_val <- rep(alpha_i,length(meth_u))
  names(alpha_val)<-meth_u
  names_alpha<-names(alpha_val)
  alpha_val[startsWith(names_alpha, 'marg')] <- alpha_m
  
  cbp2 <- c("#F0E442","#56B4E9", "#009E73","#440154FF",
            "#0072B2", "#D55E00","#CC79A7")
  names(cbp2)<-c("marg_stack", 
                 "marg_gee",
                 "marg_1stfix",
                 "marg_1stmar",
                 "marg_2stpar",
                 "marg_2stpar",
                 "marg_2stpw")
  
  if(any(startsWith(names_alpha, 'ind'))){
    cbp1 <- rep("#999999",length(names_alpha))
    cbp1[startsWith(names_alpha, 'marg')] <-cbp2[ names_alpha[startsWith(names_alpha, 'marg')]]
  }else if(any(startsWith(names_alpha, 'marg'))){
    cbp1 <- cbp2
  }else{
    cbp1 <- c("#440154FF","#21908CFF","#FDE725FF")
  }
  
  p1<- ggplot(plot_data,aes(x=pi,y=pr))+
    geom_abline(intercept = 0, slope = 1,linetype=2) +
    geom_line(linewidth = 0.5,aes(group=method,color=method))+
    geom_ribbon(aes(ymin=lower, ymax=upper,fill=method,group=method,alpha=method),colour = NA) + 
    labs( x = TeX(r"( Estimated probability, $\pi_{k}$ )"),
          y = TeX(r"( Observed probability, $\rho_{k}$ )")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 12),legend.position = "bottom")+
    scale_y_continuous(limits = c(-0.25,1))+
    scale_alpha_manual(values=alpha_val)+
    geom_segment(aes(x=0,y=-0.15,xend=1,yend=-0.15),color="black",linewidth = 0.1)+
    geom_segment(data=pi_dens,aes(x=pi,y=min,xend=pi,yend=max),color="black",linewidth = 1)+
    annotate("text", x=0, y=-0.1, label= "1",size=2) + 
    annotate("text", x=0, y=-0.2, label= "0",size=2) +
    scale_fill_manual(values = cbp1)+
    scale_colour_manual(values = cbp1)
  if(cluster){
    p1<-p1+ facet_grid(.~ model)+theme_classic()
  }
  
  if(any(startsWith(names_alpha, 'ind'))){
    p1<- p1+ theme(legend.position = "none")
  }
  return(p1)
}


# Marginal Stacked calibration plots ----

#' Marginal Stacked calibration plots
#'
#' @param pred_data calibration information (yij, pij) dataset for all the i studies
#' @param new_df new dataset with the given p_k model estimated prediction values.
#' @param cal.mod  calibration model form
#' @param alpha  confidence level
#'
#' @return dataset with calibration plot point estimated values (rho_k,p_k)  and CI for the stacked method. 
#' 
marg_stack <- function(pred_data, new_df, cal.mod = "y ~ ns(lpi,df=3)", alpha = 0.05){
  out <- NULL
  
  for (i in unique(pred_data$model)){
    
    pred_data_i <- pred_data|>filter(model==i)
    new_df_i <- new_df|>filter(model==i)
    
    if (cal.mod == "y ~ ns(lpi,df=3)"){ #gams
      cal_mod_stack <-  gam::gam(as.formula(cal.mod), data = pred_data_i, family = binomial("logit"))
      fprediction <- predict(cal_mod_stack, newdata = new_df_i, type="link", se=TRUE)
      lpr <- fprediction$fit
      lprse <- fprediction$se.fit
      
    }else if (cal.mod =="loess"){ #loess
      cal_mod_stack <- loess(y ~ lpi, degree = 2, data = pred_data_i,control=loess.control(surface="interpolate"))
      cal_pred <- predict(cal_mod_stack,newdata=new_df_i,se=TRUE)
      pr <- cal_pred$fit
      pr[pr <= 0] <- 0
      pr[pr >= 1] <- 1
      lpr <- qlogis(pr)
      lpr[is.infinite(lpr)|is.na(lpr)] <- NA
      lprse <- cal_pred$se.fit/(pr*(1-pr))
      lprse[is.infinite(lprse)] <- NA
      
    }else{ #locfit
      
      cal_mod <- locfit(y ~ lpi, family="binomial",data = pred_data_i, link="logit")
      cal_pred <- predict(object=  cal_mod, newdata = new_df_i, tr = function(x){x}, se=TRUE, band="local")
      # for avoiding extrapolation outside of the input prediction range
      cal_pred$fit[new_df_i$lpi < min(pred_data_i$lpi)|new_df_i$lpi > max(pred_data_i$lpi)] <- NA
      cal_pred$se.fit[new_df_i$lpi < min(pred_data_i$lpi)|new_df_i$lpi > max(pred_data_i$lpi)] <- NA
      lpr <- cal_pred$fit
      lprse <- cal_pred$se.fit
    }
    
    marg_stack_i <- new_df_i |>
      mutate(pr = plogis(lpr),
             lower = plogis(lpr+qnorm(alpha/2)*lprse),
             upper = plogis(lpr+qnorm(1-alpha/2)*lprse),
             method = "marg_stack",
             model = i)
    
    out <- rbind(out, marg_stack_i)
  }
  
  return(out)
}



# Marginal GEE calibration plots ----

#' Marginal GEE calibration plots
#'
#' @param pred_data calibration information (yij, pij) dataset for all the i studies
#' @param new_df new dataset with the given p_k model estimated prediction values.
#' @param cal.mod  calibration model form
#' @param alpha  confidence level
#'
#' @return dataset with calibration plot point estimated values (rho_k,p_k)  and CI for the GEE method. 
#' 

marg_gee <- function(pred_data, new_df, cal.mod = "y ~ ns(lpi,df=3)", alpha = 0.05){
  library(glmtoolbox)
  library(doBy)
  pred_data$cluster.num <- as.numeric(pred_data$cluster)
  
  out <- NULL
  
  for (i in unique(pred_data$model)){
    
    pred_data_i <- pred_data|>filter(model==i)
    new_df_i <- new_df|>filter(model==i)
    
    pred_data_i <- orderBy(~cluster.num, pred_data_i)
    
    gee_mod <- function(x){
      tryCatch(
        # This is what I want to do...
        {
          mod <- glmgee(as.formula(cal.mod), 
                        id = cluster.num, 
                        data = x,
                        corstr="Exchangeable",
                        family = binomial(link = "logit"))
          return(mod)
        },
        # ... but if an error occurs, tell me what happened: 
        warning = function(w) {
          return(NA)
        }
      )
    }
    
    cal_mod_gee <-  gee_mod( x = pred_data_i)
    if(!all(is.na(cal_mod_gee))){
      pred_gee  <-  as.data.frame(predict(cal_mod_gee, newdata = new_df_i,se.fit = T))
      
      marg_gee_i <- new_df_i |>
        mutate(pr = plogis(pred_gee[,1]),
               lower = plogis(pred_gee[,1]+ qnorm(alpha/2)*pred_gee[,2]),
               upper = plogis(pred_gee[,1]+ qnorm(1-alpha/2)*pred_gee[,2]),
               method = "marg_gee",
               model = i)
      out <- rbind(out,marg_gee_i)
    }
    
  }
  
  return(out)
}  


# Marginal 1 step method calibration plots ----

#' Marginal 1 step method calibration plots
#'
#' @param pred_data calibration information (yij, pij) dataset for all the i studies
#' @param new_df new dataset with the given p_k model estimated prediction values.
#' @param cal.mod  calibration model form
#' @param alpha  confidence level
#'
#' @return dataset with calibration plot point estimated values (rho_k,p_k)  and CI for the 1 step method method. 
#' 

marg_1stp <- function(pred_data, new_df, cal.mod ="y ~  ns(lpi,df=3)", alpha = 0.05){
  library("GLMMadaptive")
  out <- NULL
  
  for (i in unique(pred_data$model)){
    
    pred_data_i <- pred_data|>filter(model==i)
    new_df_i <- new_df|>filter(model==i)
    
    cal_mod_mix <- mixed_model(fixed = as.formula(cal.mod),
                               random = ~ 1| cluster, 
                               data = pred_data_i, family = binomial("logit"))
    
    # Fixed effect prediction (b=0)
    mix_fix  <-  predict(cal_mod_mix, newdata = new_df_i, type_pred = "link", type = "mean_subject", se.fit = TRUE)
    
    marg_mixfix <- new_df_i |>
      mutate(pr = plogis(mix_fix$pred),
             lower = plogis(mix_fix$pred+ qnorm(alpha/2)*mix_fix$se.fit),
             upper = plogis(mix_fix$pred+ qnorm(1-alpha/2)*mix_fix$se.fit),
             method = "marg_1stfix",
             model = i)
    
    out <- rbind(out, marg_mixfix)
    
    # Population averaged predictions (marginalized coefficients)
    mix_mar <- predict(cal_mod_mix, newdata = new_df_i, type_pred = "link", type = "marginal", se.fit = TRUE)
    
    marg_mixmar <- new_df_i |>
      mutate(pr = plogis(mix_mar$pred),
             lower = plogis(mix_mar$pred + qnorm(alpha/2)*mix_mar$se.fit),
             upper = plogis(mix_mar$pred + qnorm(1-alpha/2)*mix_mar$se.fit),
             method = "marg_1stmar",
             model = i)
    
    out <- rbind(out, marg_mixmar)
  }
  
  return(out)
}  



marg_1stp_s <- function(pred_data, new_df, cal.mod="y ~  ns(lpi,df=3)", alpha = 0.05){
  library("GLMMadaptive")
  out <- NULL
  
  for (i in unique(pred_data$model)){
    
    pred_data_i <- pred_data|>filter(model==i)
    new_df_i <- new_df|>filter(model==i)
    
    cal_mod_mix <- mixed_model(fixed = as.formula(cal.mod),
                               random = ~ ns(lpi,df=3)| cluster, 
                               data = pred_data_i, family = binomial("logit"))
    
    # Fixed effect prediction (b=0)
    mix_fix  <-  predict(cal_mod_mix, newdata = new_df_i, type_pred = "link", type = "mean_subject", se.fit = TRUE)
    
    marg_mixfix <- new_df_i |>
      mutate(pr = plogis(mix_fix$pred),
             lower = plogis(mix_fix$pred+ qnorm(alpha/2)*mix_fix$se.fit),
             upper = plogis(mix_fix$pred+ qnorm(1-alpha/2)*mix_fix$se.fit),
             method = "marg_1stfix",
             model = i)
    
    out <- rbind(out, marg_mixfix)
    
    # Population averaged predictions (marginalized coefficients)
    mix_mar <- predict(cal_mod_mix, newdata = new_df_i, type_pred = "link", type = "marginal", se.fit = TRUE)
    
    marg_mixmar <- new_df_i |>
      mutate(pr = plogis(mix_mar$pred),
             lower = plogis(mix_mar$pred + qnorm(alpha/2)*mix_mar$se.fit),
             upper = plogis(mix_mar$pred + qnorm(1-alpha/2)*mix_mar$se.fit),
             method = "marg_1stmar",
             model = i)
    
    out <- rbind(out, marg_mixmar)
  }
  
  return(out)
}  


# Marginal 2 step method parameters calibration plots ----

#' Marginal 2 step parameters method calibration plots
#'
#' @param pred_data calibration information (yij, pij) dataset for all the i studies
#' @param new_df new dataset with the given p_k model estimated prediction values.
#' @param df new  degress of freedom of the natural splines
#' @param alpha  confidence level
#'
#' @return dataset with calibration plot point estimated values (rho_k,p_k)  and CI for the 2 step parameters method method. 
#' 
marg_2stpar <- function(pred_data, new_df, df = 3, alpha = 0.05){
  out <- NULL
  for (i in unique(pred_data$model)){
    
    pred_data_i <- pred_data|>filter(model==i)
    new_df_i <- new_df|>filter(model==i)
    
    #  Set a weight to observable data ----
    pred_data_i$weights <- 1
    
    #  Augmented dataset to assure same range of pi on all clusters----
    maxq <- quantile(pred_data_i$lpi, 0.95)
    minq <- quantile(pred_data_i$lpi, 0.05) 
    
    # add observations close to the boundaries to each cluster with a small weight
    bounds <- pred_data_i|>
      filter(lpi >= maxq|lpi <= minq)|>
      mutate( weights = 0.0000001 )
    
    nameclus <- unique(pred_data_i$cluster)
    bound_data <- do.call("rbind", replicate(n=length(nameclus), bounds,simplify = FALSE))
    bound_data$cluster <- rep( nameclus, each=nrow(bounds))
    aug_data <- rbind(pred_data_i,bound_data)
    
    #  Calculate the knots ----
    
    # Using fixed points
    knots_all <-with(aug_data, attr(ns(lpi, df = df),"knots"))
    
    # Data to be predicted ----
    n_weights <- rep(1,nrow(new_df)) 
    
    # 1 step: Run model and predict ----
    model_info <- aug_data|>
      nest(data = -cluster)|>
      mutate(weights = map(data, ~.x%>%pull(weights)),
             model = map2(data, weights, ~tryCatch(
               expr = {gam::gam(y ~ ns(lpi, knots = knots_all), weights = .y, 
                                family = binomial("logit"), data=.x)},
               error = function(e){NA})),
             n_vec = map(data,function(x){x%>%filter(weights==1)%>%nrow()}),         
             coef  = map(model,~.x$coefficients),
             cov_vec = map(model, ~ as.vector(vcov(.x)[lower.tri(vcov(.x), diag = T)])))
    
    # 2 step: Apply multivariate random-effects meta-analysis ----
    m_n <- do.call("rbind",model_info$n_vec) 
    m_coef <- do.call("rbind",model_info$coef) 
    m_cov <- do.call("rbind",model_info$cov_vec)   
    fit <- suppressWarnings(try(mixmeta::mixmeta( m_coef, m_cov, method = "reml",
                                                  control = list(hessian = TRUE)), silent = TRUE))
    
    
    # Get predictions        
    mm <- model.matrix(~ns(lpi, knots = knots_all),new_df_i) # Get model matrix for prediction values 
    tr_obs <- as.numeric(mm %*% as.vector(fit[["coefficients"]]))
    tr_se <- sqrt(diag(mm %*% as.matrix(vcov(fit)) %*% t(mm)))
    
    marg_2stpar_i <- new_df_i |>
      mutate(pr = plogis(tr_obs),
             lower = plogis(tr_obs + qnorm(alpha/2)*tr_se),
             upper = plogis(tr_obs + qnorm(1-alpha/2)*tr_se),
             method = "marg_2stpar") 
    out <- rbind(out,marg_2stpar_i)
  }
  return(out)   
}

# Marginal 2 step method predicted values calibration plots ----

#' Marginal 2 step  predicted values method calibration plots
#'
#' @param pred_data calibration information (yij, pij) dataset for all the i studies
#' @param new_df new dataset with the given p_k model estimated prediction values.
#' @param cal.mod  calibration model form
#' @param alpha  confidence level
#'
#' @return dataset with calibration plot point estimated values (rho_k,p_k)  and CI for the 2 step  predicted values method method. 
#' 

marg_2stpw <- function(pred_data, new_df, cal.mod = "y~ns(lpi,df=3)", alpha = 0.05){
  
  out <- NULL
  
  for (i in unique(pred_data$model)){
    
    pred_data_i <- pred_data|>filter(model==i)
    new_df_i <- new_df|>filter(model==i)
    
    # 1-step predict and estimate observed probabilities ---- 
    
    model_info <- pred_data_i|>
      nest(data = -cluster)|>
      mutate(model = map(data, ~tryCatch(
        expr = {if (cal.mod == "y~ns(lpi,df=3)"){gam::gam(as.formula(cal.mod),
                                                          family = binomial("logit"), data=.x)}
          else if (cal.mod =="loess"){ #loess
            loess(y ~ lpi, degree = 2, data =.x,control=loess.control(surface="interpolate"))}
          else{
            locfit(y ~ lpi, family="binomial",data =.x, link="logit")}} ,
        error = function(e){NA})))|>
      mutate(mod_predict = map2(model,data,function(x,y){
        # reduce data to predict at model range
        new_data_i<-new_df_i|>filter(lpi>=min(y$lpi)&lpi<=max(y$lpi))
        if(cal.mod == "y~ns(lpi,df=3)"){
          pred <- predict(x, newdata = new_data_i, type="link",se = TRUE)
          new_data_i$tr_obs <- pred$fit
          new_data_i$tr_se <-  pred$se.fit
        }else if( cal.mod == "loess"){
          pred <- predict(x,newdata=new_data_i,se=TRUE)
          pr <- pred$fit
          pr[pr <= 0] <- 0
          pr[pr >= 1] <- 1
          lpr <- qlogis(pr)
          lpr[is.infinite(lpr)|is.na(lpr)] <- NA
          lprse <- pred$se.fit/(pr*(1-pr))
          lprse[is.infinite(lprse)] <- NA
          new_data_i$tr_obs <- lpr
          new_data_i$tr_se <-  lprse
        }else{
          pred <- predict(object=x, newdata = new_data_i, tr = function(x){x}, se=TRUE, band="local")
          new_data_i$tr_obs <- pred$fit
          new_data_i$tr_se <-  pred$se.fit
        }
        
        return(new_data_i)
      }))|>
      select(cluster,mod_predict)|>
      unnest(mod_predict)|>
      nest(data = -lpi)
    
    
    # 2step: Apply univariate random-effects meta-analysis
    
    marg_2stpw_i <-  model_info%>%
      mutate(mod_data = map(data,function(x){
        metamod <- tryCatch(
          expr = { metafor::rma(yi=as.vector(x$tr_obs), sei=as.vector(x$tr_se),
                                method="REML",test= "knha")},
          error = function(e){ tryCatch( expr = {metafor::rma(yi=as.vector(x$tr_obs),    
                                                              sei=as.vector(x$tr_se), 
                                                              method="REML",
                                                              test= "knha",control=list(stepadj=0.5))},
                                         error = function(e){NA})})}))|>
      mutate(out_data=map(mod_data,function(x){
        resp <- data.frame(
          pr = plogis(as.numeric(x$beta)),
          lower = plogis(as.numeric(x$ci.lb)), # with t distribution instead of Z
          upper = plogis(as.numeric(x$ci.ub)))
      }))|>
      select(c(lpi,out_data))|>
      unnest(col=c(out_data))|>
      mutate(method = "marg_2stpw",
             model = i,
             pi =plogis(lpi)
      )|>
      select(c("model", "pi","lpi","pr","lower","upper","method"))
    out <- rbind(out,marg_2stpw_i)
  }
  return(out)
}


