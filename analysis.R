library(cmdstanr)
library(rstan)
library(gridExtra)
library(ggplot2)
library(deSolve)
library(latex2exp)

## Load New York City Cases and raw data visualization
#https://github.com/nychealth/coronavirus-data
nyc_data <- read.csv("boroughs-case-hosp-death.csv")
# outcome y1 will be cases starting March 8th, 2020 to align with Apple mobility data and state of emergency declared in NYC
n_obs=90
y1 = rowSums(cbind(nyc_data$BK_CASE_COUNT,
                   nyc_data$QN_CASE_COUNT,
                   nyc_data$MN_CASE_COUNT,
                   nyc_data$BX_CASE_COUNT,
                   nyc_data$SI_CASE_COUNT))[7:(7+n_obs-1)]
y1full = rowSums(cbind(nyc_data$BK_CASE_COUNT,
                       nyc_data$QN_CASE_COUNT,
                       nyc_data$MN_CASE_COUNT,
                       nyc_data$BX_CASE_COUNT,
                       nyc_data$SI_CASE_COUNT))[7:(7+n_obs-1)]
N = sum(2559903, 2253858, 1628706, 1418207, 476143)

## Mobility Time Series Data: 
# load Apple mobility data
applemobilitytrends <- read.csv("applemobilitytrends.csv")
# filter by NYC data
nyc_region = which(applemobilitytrends$region=="New York")
# remove row labels. dat now contains driving, walking, transit % change
dat = applemobilitytrends[nyc_region,-c(1:7)]
# outcome for analysis will be first 100 days of transit change after March 7th
start_date = which(colnames(dat)=="X2020.03.08") # start on 8th
y2 = as.numeric((dat[3,start_date:(start_date+n_obs-1)]))/100
df = data.frame("Day"=1:n_obs,
                "Cases" = y1,
                "Apple Transit Mobility" = y2)
plot1 <- ggplot(data=df, aes(x = Day, y = Cases)) + ggtitle("Covid-19 Case Counts")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal() +
  geom_point(size=2)+theme_minimal(base_size = 15)+ylim(c(0,8500))+
  theme(plot.title = element_text(size=30, hjust=.5))
plot1

plot2 <- ggplot(data=df, aes(x = Day, y = Apple.Transit.Mobility)) + ggtitle("Apple Transit Mobility")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal() +
  geom_point(size=2)+theme_minimal(base_size = 15)+ylim(c(0,1))+
  theme(plot.title = element_text(size=30, hjust=.5))
plot2

###### Inference with vanilla SIR model ###### 
model_data <- list(
  n_obs = n_obs,
  N = N, 
  t0=0,
  cases = df$Cases,
  movement_data=1-df$Apple.Transit.Mobility)
file <- file.path("sir.stan")
sir.mod <- cmdstan_model(file)
sir_mcmc <- sir.mod$sample(data = model_data,
                           max_treedepth = 10,adapt_delta = .8,
                           chains = 1,
                           iter_warmup = 500,
                           iter_sampling = 500,
                           refresh = 100, parallel_chains = getOption("mc.cores", 4))
stanfit <- rstan::read_stan_csv(sir_mcmc$output_files())
a=as.data.frame(stanfit)
findstr = paste0("pred_cases\\[\\d+\\]")
locs = grepl(findstr, names(a))
sir_pred <- cbind(as.data.frame(summary(
  stanfit, pars = names(a)[locs], probs = c(0.05, 0.5, 0.95))$summary), 1:model_data$n_obs)
colnames(sir_pred) <- make.names(colnames(sir_pred)) 


###### Inference with SLIR model ######  
file <- file.path("model.stan")
slir.mod <- cmdstan_model(file)
slir_mcmc <- slir.mod$sample(data = model_data,seed=123,
                             max_treedepth = 10,adapt_delta = .8,
                             chains = 4, 
                             iter_warmup = 500,
                             iter_sampling = 500,
                             refresh = 250, parallel_chains = getOption("mc.cores", 4))

## Diagnostics of SLIR model fit
stanfit.slir <- rstan::read_stan_csv(slir_mcmc$output_files())
summary(stanfit.slir, pars = c("R0","a","b","gamma","phi"))$summary
traceplot(stanfit.slir, pars = c("R0","a","b","gamma"))

###### Plot predictive distributions  ###### 
a=as.data.frame(stanfit.slir)
findstr = paste0("pred_cases\\[\\d+\\]")
locs = grepl(findstr, names(a))
slir_pred <- cbind(as.data.frame(summary(
  stanfit.slir, pars = names(a)[locs], probs = c(0.05, 0.5, 0.95))$summary), 1:model_data$n_obs)
colnames(slir_pred) <- make.names(colnames(slir_pred)) 
p = ggplot() + ggtitle("NYC Cases")
datFull2 = y1full
sir_pred[,'X50.'][50:n_obs]=NA
p + geom_line(data = sir_pred, mapping = aes(x=1:model_data$n_obs,y=X50., color="orange"), size=3)+
  geom_line(data = slir_pred, mapping = aes(x=1:model_data$n_obs,y=X50., color="cornflowerblue"), size=3)+
  geom_ribbon(data = slir_pred,aes(x=1:model_data$n_obs,ymin= X5., ymax=X95.),alpha=.3, fill="cornflowerblue")+
  xlab("Days") + ylab("Counts")+
  geom_point(data = df, aes_string(x="Day",y="Cases"), size=3,alpha=.5)+
  theme_minimal(base_size = 20)+theme(legend.position="bottom")+
  scale_color_identity(name = "",
                       breaks = c("orange", "cornflowerblue"),
                       labels = c("SIR Fit", "SLIR Prior Predictive"),
                       guide = "legend")+
  guides(colour = guide_legend(override.aes = list(size = 3)))+
  coord_cartesian(xlim=c(0,length(y1)),
                  ylim = c(0,max(slir_pred$X95.)),
                  expand = TRUE,
                  default = FALSE,
                  clip = "on"
  )+theme(plot.title = element_text(size=30, hjust=.5))

findstr = paste0("pred_movement\\[\\d+\\]")
locs = grepl(findstr, names(a))
slir_pred <- cbind(as.data.frame(summary(
  stanfit.slir, pars = names(a)[locs], probs = c(0.05, 0.5, 0.95))$summary), 1:model_data$n_obs)
colnames(slir_pred) <- make.names(colnames(slir_pred))
p = ggplot() + ggtitle("NYC Transit Reduction")
datFull2 = model_data$movement_data
p + geom_line(data = slir_pred, mapping = aes(x=1:model_data$n_obs,y=X50., color="cornflowerblue"), size=3)+
  geom_ribbon(data = slir_pred,aes(x=1:model_data$n_obs,ymin= X5., ymax=X95.),alpha=.3, fill="cornflowerblue")+
  xlab("Days") + ylab("% in L Compartment")+
  geom_point(data = as.data.frame(model_data), aes_string(1:n_obs,y="movement_data"), size=3,alpha=.5)+
  scale_color_identity(name = "",
                       breaks = c("cornflowerblue"),
                       labels = c("SLIR Fit"),
                       guide = "legend")+
  theme_minimal(base_size = 20)+theme(legend.position="bottom")+
  guides(colour = guide_legend(override.aes = list(size = 3)))+
  theme(plot.title = element_text(size=30, hjust=.5))


findstr = paste0("R_effective\\[\\d+\\]")
locs = grepl(findstr, names(a))
slir_pred <- cbind(as.data.frame(summary(
  stanfit.slir, pars = names(a)[locs], probs = c(0.05, 0.5, 0.95))$summary), 1:model_data$n_obs)
colnames(slir_pred) <- make.names(colnames(slir_pred)) 
p = ggplot() + ggtitle(TeX('$\\mathit{R}_0$ and $\\mathit{R}_t$'))#
datFull2 = model_data$movement_data
labs = list(TeX('$\\mathit{R}_0'),TeX('$\\mathit{R}_t$'))
nm = names(df)[-1]
p + geom_line(data = slir_pred, mapping = aes(x=1:model_data$n_obs,y=X50., color="re"), size=1.5)+
  geom_ribbon(data = slir_pred,aes(x=1:model_data$n_obs,ymin= X5., ymax=X95.),alpha=.3, fill="darkorange")+
  xlab("Days") + ylab(TeX('$\\mathit{R}$'))+ylim(c(0,7))+
  geom_point(aes(x=0,y=mean(slir_mcmc$draws("R0")),color='r0'),size=3)+
  geom_segment(aes(x=0,xend=0,y=quantile(slir_mcmc$draws("R0"),probs = .025), yend=quantile(slir_mcmc$draws("R0"),probs = .975),color='r0'),size=5,alpha=.3)+
  theme_minimal(base_size = 20)+
  scale_colour_manual(name = '', 
                      values =c('r0'='firebrick2','re'='darkorange'),labels=labs)+
  theme(legend.position="bottom")+geom_hline(aes(yintercept = 1),alpha=.6,linetype=2,size=1)+
  guides(colour = guide_legend(override.aes = list(size = 5,linetype=c(0,0))))+
  theme(plot.title = element_text(size=30, hjust=.5))

