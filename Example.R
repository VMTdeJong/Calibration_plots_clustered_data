source("Functions.R") # functions to plot overall calibration curves
source("Load data.R") # load the data sets


#Table to display cluster characteristics ----
data <- imp |>
  filter(Imputation == 1 & clus!="Oudega study")|>
  mutate(cluster = factor(ifelse( as.numeric(clus) > 8,
                                  paste0("Study ", as.numeric(clus)-1), 
                                  paste0("Study ", as.numeric(clus))) , levels = paste0("Study ",1:12)),
         dvt = as.factor(ifelse(dvt==2,"yes","no")),
         y = ifelse(dvt=="yes",1,0),
         sex = as.factor(ifelse(sex==2,"male","female")),
         malign = factor(as.factor(ifelse(malign==2,"active","no-active")),levels=c("no-active","active")),
         surg= as.factor(ifelse(surg==2,"yes","no")),
         vein= as.factor(ifelse(vein==2,"yes","no")),
         calfdif3=as.factor(ifelse(calfdif3==2,"yes","no")))



tab1<-table1(~ dvt + ddimdich + calfdif3 + oachst + sex + notraum + vein + malign + surg | cluster,
             data = data, 
             caption = "Baseline characteristics of the DVT dataset.")
t1kable(tab1, booktabs = TRUE)


# CPM predictions ----

betas4 <- c(-4.25, 2.46, 1.15, 0.72, 0.72,    0,  0,   0,    0)
betas8 <- c(-5.02, 2.43, 1.15, 0.76, 0.71, 0.67,0.53,0.50, 0.42)

data <- data |>
  group_split(cluster)|>
  as_tibble_col()|> 
  rename( pred_data = value)|> # make prediction on each study dataset
  mutate( pred_data = map(pred_data, function(x){
    xmat <- as.matrix(model.matrix(~ ddimdich + calfdif3 + oachst + sex + notraum + vein + malign + surg, x))
    x$lpi_4<- as.vector(xmat %*% betas4)
    x$lpi_8 <-as.vector(xmat %*% betas8)
    x<- pivot_longer(x,cols = starts_with("lpi"),names_to = "model",names_prefix = "lpi_",values_to = "lpi")
    x$pi <-plogis(x$lpi)
    x}))|>
  mutate( cluster = map(pred_data,function(x){as.character(x[[1,"cluster"]])}))



# Cluster-specific plots ----

## Set the range of prediction ----
## Stack pred_data  
pred_data <- data|>
  select(pred_data)|>
  unnest(cols = c(pred_data))

## Get the range of prediction
list_new_df <- lapply(unique(pred_data$model), function(i){
  pred_data_i <- pred_data|>filter(model==i)
  new_df<- data.frame(model=i,pi=seq(min(pred_data_i$pi),max(pred_data_i$pi),length.out = 100 ) )})

new_df <- do.call(rbind,list_new_df)
new_df$lpi <-qlogis(new_df$pi)


## Plot density predicted probability vs real value ----

#data %<>% mutate(pi_dens = map(pred_data,~get_dens(data=.x,new_df=new_df))) # get density of pi at cluster level
pi_dens <- get_dens(data=pred_data,new_df=new_df) # get density of pi at marginal level
data %<>% mutate(plot_cal = map(pred_data,~get_calmod(data=.x,new_df=new_df,df=3))) # get the observed pr with CI at cluster level

plot_data <- data|>
  select(c(cluster, plot_cal))|>
  unnest(cols = c(cluster, plot_cal)) |>
  mutate(cluster = factor(cluster, levels = paste("Study", 1:12)))


cal_plot_ind(plot_data,pi_dens)

ggsave("figures/fig_1_cluster_specific.pdf")


# Overall calibration curves ----

## Stack ----
 # gam
dvt_marg_stack <- marg_stack(pred_data, new_df)

plot_marg_stack <- plot_data|>
  filter(method == "splines")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_stack)

cal_plot(plot_data = plot_marg_stack, pi_dens, cluster = TRUE)

ggsave("figures/fig_2_marg_stack_gam.pdf")

 # loess
dvt_marg_stackloe <- marg_stack(pred_data, new_df, cal.mod = "loess")

plot_margloe <- plot_data|>
  filter(method == "loess")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_stackloe)

cal_plot(plot_data = plot_margloe, pi_dens, cluster = TRUE)

ggsave("figures/fig_a1_stack_loess.pdf")

 # locfit
dvt_marg_stacklof <-  marg_stack(pred_data, new_df, cal.mod = "locfit")
plot_marglof <- plot_data|>
  filter(method == "locfit")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_stacklof)

cal_plot(plot_data = plot_marglof, pi_dens, cluster = TRUE)

ggsave("figures/fig_a2_stack_locfit.pdf")

## GEE ----
dvt_marg_gee <- marg_gee(pred_data,new_df)

plot_marg_gee <- plot_data|>
  filter(method == "splines")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_gee)

cal_plot(plot_data = plot_marg_gee, pi_dens, cluster = TRUE)

ggsave("figures/fig_3_gee_gam.pdf")

## 1step ----

dvt_marg_1stp <- marg_1stp(pred_data,new_df)

plot_marg_1stp <- plot_data|>
  filter(method == "splines")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_1stp)

cal_plot(plot_data = plot_marg_1stp, pi_dens, cluster = TRUE)

ggsave("figures/fig_4_1step_gam.pdf")

## 2step----

### parameters----
dvt_marg_2stpar <- marg_2stpar(pred_data, new_df)    

plot_marg <- plot_data|>
  filter(method == "splines")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_2stpar)

cal_plot(plot_data = plot_marg, pi_dens, cluster = TRUE)

ggsave("figures/fig_5_2step_gam.pdf")


plot_marg_loe <- plot_data|>
  filter(method == "loess")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_2stpar)

cal_plot(plot_data = plot_marg_loe, pi_dens, cluster = TRUE)

ggsave("figures/fig_a3_2step_loess.pdf")


plot_marg_loc <- plot_data|>
  filter(method == "locfit")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_2stpar)

cal_plot(plot_data = plot_marg_loc, pi_dens, cluster = TRUE)

ggsave("figures/fig_a4_2step_locfit.pdf")

### predictions----
dvt_marg_2stpw <- marg_2stpw(pred_data, new_df)    

plot_marg <- plot_data|>
  filter(method == "splines")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_2stpw)

cal_plot(plot_data = plot_marg, pi_dens, cluster = TRUE)  

ggsave("figures/fig_6_2step_pred_gam.pdf")


# loess
marg_2stpwloe <- marg_2stpw(pred_data, new_df, cal.mod = "loess")     

plot_margloe <- plot_data|>
  filter(method == "loess")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(marg_2stpwloe)

cal_plot(plot_data = plot_margloe, pi_dens, cluster = TRUE)  

ggsave("figures/fig_a5_2step_pred_loess.pdf")

# locfit
dvt_marg_2stpwloc <- marg_2stpw(pred_data, new_df, cal.mod = "locfit")  
plot_margloc <- plot_data|>
  filter(method == "locfit")|>
  mutate(method = paste0("ind_", cluster))|>
  bind_rows(dvt_marg_2stpwloc)

cal_plot(plot_data = plot_margloc, pi_dens, cluster = TRUE)

ggsave("figures/fig_a6_2step_pred_locfit.pdf")

## all---

marg_all <- do.call("rbind", 
                    list(dvt_marg_stack, 
                         dvt_marg_gee, 
                         dvt_marg_1stp_s, 
                         dvt_marg_2stpar, 
                         dvt_marg_2stpw))|>
  mutate(method = factor(method, 
                         levels = c("marg_stack", 
                                    "marg_gee", 
                                    "marg_1stfix", 
                                    "marg_1stmar", 
                                    "marg_2stpar",
                                    "marg_2stpw")))

cal_plot(plot_data = marg_all, 
         pi_dens,
         cluster = TRUE, 
         alpha_i = 0, 
         alpha_m = 0.1)  

ggsave("figures/fig_7_all.pdf")
