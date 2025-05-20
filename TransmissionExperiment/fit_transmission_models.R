library(tidyverse)
library(rstan)
library(loo)
library(bayesplot)
library(gridExtra)

options(mc.cores = parallel::detectCores()) #run cores in parallel
rstan_options(auto_write = TRUE) #auto-write compiled code to hard drive

source("functions.R")

`%ni%` <- Negate(`%in%`)

input_data <- read_csv("transmission_fitting/data/STAN_input_data.csv")

input_data <- input_data %>% filter(strain %ni% c('H07','MEW'), alive>0) %>% select(-isolate)

snpv_data <- input_data  %>% filter(capsid == "SNPV") %>% mutate(isolate = as.numeric(as.factor(strain)))
mnpv_data <- input_data %>% filter(capsid == "MNPV") %>% mutate(isolate = as.numeric(as.factor(strain)))

snpv_data2 <- snpv_data %>% mutate(isolate = isolate + 5)

input_data2 <- rbind(mnpv_data,snpv_data2)
input_data3 <- input_data2 %>% mutate(isolate = ifelse(capsid == "SNPV",isolate -5, isolate))

n_obs <- nrow(input_data2)
n_trees <- length(unique(input_data2$tree))

df_M <- list(N = n_obs, virus = input_data2$virus, 
               alive = input_data2$alive, 
               pathogen = input_data2$P0)

df_M2 <- list(N = n_obs, virus = input_data2$virus, 
             alive = input_data2$alive, 
             pathogen = input_data2$P0,
             morph = as.numeric(as.factor(input_data2$capsid)))

df_M2xT <- list(N = n_obs, virus = input_data2$virus, T = n_trees,
              alive = input_data2$alive, 
              pathogen = input_data2$P0,
              morph = as.numeric(as.factor(input_data2$capsid)),
              trees = as.numeric(as.factor(input_data2$tree)))

df_SxM <- list(N = n_obs, I = 10, virus = input_data2$virus, 
               alive = input_data2$alive, 
               pathogen = input_data2$P0,
               isolate = input_data2$isolate)

df_SxMxT_all <- list(N = n_obs, T = n_trees, I = 5,
                 virus = input_data3$virus,
                 alive = input_data3$alive,
                 pathogen = input_data3$P0,
                 trees = as.numeric(as.factor(input_data3$tree)),
                 isolate = input_data3$isolate,
                 morph = as.numeric(as.factor(input_data3$capsid)))

################################
#
# Parameter values for each model
#
################################

mod0_params <- c('nu_bar','C')
mod1_params <- c('mu','muC',paste0('nu_bar[',seq(1:10),']'),paste0('C[',seq(1:10),']'))
SNPV_mod3_params <- c('mu_S','muC_S',paste0('nu_bar_S[',seq(1:5),']'),paste0('C_S[',seq(1:5),']'))
MNPV_mod3_params <- c('mu_M','muC_M',paste0('nu_bar_M[',seq(1:5),']'),paste0('C_M[',seq(1:5),']'))

SNPV_mod5_params <- c('mu_S[1]','mu_S[2]','muC_S[1]','muC_S[2]',paste0('nu_bar_S[',seq(1:5),',1]'),paste0('C_S[',seq(1:5),',1]'),
                      paste0('nu_bar_S[',seq(1:5),',2]'),paste0('C_S[',seq(1:5),',2]'))
MNPV_mod5_params <- c('mu_M[1]','mu_M[2]','muC_M[1]','muC_M[2]',paste0('nu_bar_M[',seq(1:5),',1]'),paste0('C_M[',seq(1:5),',1]'),
                      paste0('nu_bar_M[',seq(1:5),',2]'),paste0('C_M[',seq(1:5),',2]'))

SNPV_mod6_params <- c('mu_S','muC_S[1]','muC_S[2]',paste0('nu_bar_S[',seq(1:5),']'),paste0('C_S[',seq(1:5),',1]'),
                      paste0('C_S[',seq(1:5),',2]'))
MNPV_mod6_params <- c('mu_M','muC_M[1]','muC_M[2]',paste0('nu_bar_M[',seq(1:5),']'),paste0('C_M[',seq(1:5),',1]'),
                      paste0('C_M[',seq(1:5),',2]'))

SNPV_mod7_params <- c('mu_S[1]','mu_S[2]','muC_S',paste0('nu_bar_S[',seq(1:5),',1]'),paste0('C_S[',seq(1:5),']'),
                      paste0('nu_bar_S[',seq(1:5),',2]'))

MNPV_mod7_params <- c('mu_M[1]','mu_M[2]','muC_M',paste0('nu_bar_M[',seq(1:5),',1]'),paste0('C_M[',seq(1:5),']'),
                      paste0('nu_bar_M[',seq(1:5),',2]'))

################################
#
#
#
# STARTING OVER WITH CORRECT K
#
#
#
################################

niter = 20000
nchains = 5
nwarmup = niter*0.25

fit_mod1_K <- rstan::stan(file = 'models/Mod1_S_K.stan',
                          data = df_SxM,
                          iter = niter,
                          chains = nchains,
                          warmup = nwarmup,
                          refresh=0, control = list(adapt_delta = 0.95,
                                                    max_treedepth = 12)) 

mod1_df <- data.frame(summary(fit_mod1_K)$summary[mod1_params,
                            c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod1_df$params <- rownames(mod1_df)

# write_csv(mod1_df,"transmission_fitting/output/mod1_params.csv")
# 
# pdf("transmission_fitting/figures/mod1_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod1_K, nu_bars = mod1_params[c(1,3:12)], Cs = mod1_params[c(2,13:22)])
# dev.off()

fit_mod0_K <- rstan::stan(file = 'models/Mod0_S_K.stan',
                          data = df_M,
                          iter = niter,
                          chains = nchains,
                          warmup = nwarmup,
                          refresh=0, control = list(adapt_delta = 0.95,
                                                    max_treedepth = 12)) 

summary(fit_mod0_K)$summary[c(mod0_params),
                            c('mean','sd','se_mean','25%','75%')]

mod0_df <- data.frame(summary(fit_mod0_K)$summary[mod0_params,
                                                  c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])

mod0_df$params <- rownames(mod0_df)

# write_csv(mod0_df,"transmission_fitting/output/mod0_params.csv")
# 
# pdf("transmission_fitting/figures/mod0_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod0_K, nu_bars = mod0_params[1], Cs = mod0_params[2])
# dev.off()

fit_mod33_K <- rstan::stan(file = 'models/Mod33_SxMxT_K.stan',
                               data = df_SxMxT_all,
                               iter = niter,
                               chains = 3,
                               warmup = nwarmup,
                               refresh=0, control = list(adapt_delta = 0.95,
                                                         max_treedepth = 12)) 

summary(fit_mod33_K)$summary[c(SNPV_mod3_params,MNPV_mod3_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod33_df <- data.frame(summary(fit_mod33_K)$summary[c(SNPV_mod3_params,MNPV_mod3_params),
                                                  c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod33_df$params <- rownames(mod33_df)

# write_csv(mod33_df,"transmission_fitting/output/mod33_params.csv")
# 
# 
# pdf("transmission_fitting/figures/mod33_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod33_K, nu_bars = c(SNPV_mod3_params[c(1,3:7)],MNPV_mod3_params[c(1,3:7)]), Cs = c(SNPV_mod3_params[c(2,8:12)],MNPV_mod3_params[c(2,8:12)]))
# dev.off()

fit_mod55_K <- rstan::stan(file = 'models/Mod55_SxMxT_K.stan',
                               data = df_SxMxT_all,
                               iter = niter,
                               chains = nchains,
                               warmup = nwarmup,
                               refresh=0, control = list(adapt_delta = 0.95,
                                                         max_treedepth = 12)) 

mod55_df <- data.frame(summary(fit_mod55_K)$summary[c(SNPV_mod5_params,MNPV_mod5_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod55_df$params <- rownames(mod55_df)

# write_csv(mod55_df,"transmission_fitting/output/mod55_params.csv")
# 
# pdf("transmission_fitting/figures/mod55_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod55_K, nu_bars = c(SNPV_mod5_params[c(1,2,5:9,15:19)],MNPV_mod5_params[c(1,2,5:9,15:19)]),
#              Cs = c(SNPV_mod5_params[c(3:4,10:14,20:24)],MNPV_mod5_params[c(3:4,10:14,20:24)]))
# dev.off()

fit_mod66_K <- rstan::stan(file = 'models/Mod66_SxMxT_K.stan',
                               data = df_SxMxT_all,
                               iter = niter,
                               chains = nchains,
                               warmup = nwarmup,
                               refresh=0, control = list(adapt_delta = 0.95,
                                                         max_treedepth = 12)) 
mod66_df <- data.frame(summary(fit_mod66_K)$summary[c(SNPV_mod6_params,MNPV_mod6_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod66_df$params <- rownames(mod66_df)

# write_csv(mod66_df,"transmission_fitting/output/mod66_params.csv")
# 
# pdf("transmission_fitting/figures/mod66_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod66_K, nu_bars = c(SNPV_mod6_params[c(1,4:8)],MNPV_mod6_params[c(1,4:8)]),
#              Cs = c(SNPV_mod6_params[c(2,3,9:18)],MNPV_mod6_params[c(2,3,9:18)]))
# dev.off()

fit_mod77_K <- rstan::stan(file = 'models/Mod77_SxMxT_K.stan',
                               data = df_SxMxT_all,
                               iter = niter,
                               chains = nchains,
                               warmup = nwarmup,
                               refresh=0, control = list(adapt_delta = 0.95,
                                                         max_treedepth = 12)) 


mod77_df <- data.frame(summary(fit_mod77_K)$summary[c(SNPV_mod7_params,MNPV_mod7_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod77_df$params <- rownames(mod77_df)

# write_csv(mod77_df,"transmission_fitting/output/mod77_params.csv")
# 
# pdf("transmission_fitting/figures/mod77_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod77_K, nu_bars = c(SNPV_mod7_params[c(1,2,4:8,14:18)],MNPV_mod7_params[c(1,2,4:8,14:18)]),
#              Cs = c(SNPV_mod7_params[c(3,9:13)],MNPV_mod7_params[c(3,9:13)]))
# dev.off()

fit_mod35_K <- rstan::stan(file = 'models/Mod35_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = 3,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod35_K)$summary[c(SNPV_mod5_params,MNPV_mod3_params),c('mean','sd','se_mean','2.5%','97.5%','2.5%','97.5%')]

mod35_df <- data.frame(summary(fit_mod35_K)$summary[c(SNPV_mod5_params,MNPV_mod3_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod35_df$params <- rownames(mod35_df)
# 
# write_csv(mod35_df,"transmission_fitting/output/mod35_params.csv")
# 
# pdf("transmission_fitting/figures/mod35_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod35_K, nu_bars = c(SNPV_mod5_params[c(1,2,5:9,15:19)],MNPV_mod3_params[c(1,3:7)]),
#              Cs = c(SNPV_mod5_params[c(3:4,10:14,20:24)],MNPV_mod3_params[c(2,8:12)]))
# dev.off()

fit_mod36_K <- rstan::stan(file = 'models/Mod36_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = 3,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod36_K)$summary[c(SNPV_mod6_params,MNPV_mod3_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod36_df <- data.frame(summary(fit_mod36_K)$summary[c(MNPV_mod3_params,SNPV_mod6_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod36_df$params <- rownames(mod36_df)

# write_csv(mod36_df,"transmission_fitting/output/mod36_params.csv")
# 
# pdf("transmission_fitting/figures/mod36_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod36_K, nu_bars = c(SNPV_mod6_params[c(1,4:8)],MNPV_mod3_params[c(1,3:7)]),
#              Cs = c(SNPV_mod6_params[c(2,3,9:18)],MNPV_mod3_params[c(2,8:12)]))
# dev.off()

fit_mod37_K <- rstan::stan(file = 'models/Mod37_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = 3,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod37_K)$summary[c(SNPV_mod7_params,MNPV_mod3_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod37_df <- data.frame(summary(fit_mod37_K)$summary[c(MNPV_mod3_params,SNPV_mod7_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod37_df$params <- rownames(mod37_df)

# write_csv(mod37_df,"transmission_fitting/output/mod37_params.csv")
# 
# pdf("transmission_fitting/figures/mod37_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod37_K, nu_bars = c(SNPV_mod7_params[c(1,2,4:8,14:18)],MNPV_mod3_params[c(1,3:7)]),
#              Cs = c(SNPV_mod7_params[c(3,9:13)],MNPV_mod3_params[c(2,8:12)]))
# dev.off()


fit_mod53_K <- rstan::stan(file = 'models/Mod53_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = 3,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod53_K)$summary[c(SNPV_mod3_params,MNPV_mod5_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod53_df <- data.frame(summary(fit_mod53_K)$summary[c(MNPV_mod5_params,SNPV_mod3_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod53_df$params <- rownames(mod53_df)

# write_csv(mod53_df,"transmission_fitting/output/mod53_params.csv")
# 
# pdf("transmission_fitting/figures/mod53_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod53_K, nu_bars = c(SNPV_mod3_params[c(1,3:7)],MNPV_mod5_params[c(1,2,5:9,15:19)]),
#              Cs = c(SNPV_mod3_params[c(2,8:12)],MNPV_mod5_params[c(3:4,10:14,20:24)]))
# dev.off()

fit_mod63_K <- rstan::stan(file = 'models/Mod63_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = 3,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod63_K)$summary[c(SNPV_mod3_params,MNPV_mod6_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod63_df <- data.frame(summary(fit_mod63_K)$summary[c(MNPV_mod6_params,SNPV_mod3_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod63_df$params <- rownames(mod63_df)

# write_csv(mod63_df,"transmission_fitting/output/mod63_params.csv")
# 
# pdf("transmission_fitting/figures/mod63_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod63_K, nu_bars = c(SNPV_mod3_params[c(1,3:7)],MNPV_mod6_params[c(1,4:8)]),
#              Cs = c(SNPV_mod3_params[c(2,8:12)],MNPV_mod6_params[c(2,3,9:18)]))
# dev.off()

fit_mod73_K <- rstan::stan(file = 'models/Mod73_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = 3,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod73_K)$summary[c(SNPV_mod3_params,MNPV_mod7_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod73_df <- data.frame(summary(fit_mod73_K)$summary[c(MNPV_mod7_params,SNPV_mod3_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod73_df$params <- rownames(mod73_df)

# write_csv(mod73_df,"transmission_fitting/output/mod73_params.csv")
# 
# pdf("transmission_fitting/figures/mod73_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod73_K, nu_bars = c(SNPV_mod3_params[c(1,3:7)],MNPV_mod7_params[c(1,2,4:8,14:18)]),
#              Cs = c(SNPV_mod3_params[c(2,8:12)],MNPV_mod7_params[c(3,9:13)]))
# dev.off()


fit_mod56_K <- rstan::stan(file = 'models/Mod56_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = nchains,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod56_K)$summary[c(SNPV_mod6_params,MNPV_mod5_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod56_df <- data.frame(summary(fit_mod56_K)$summary[c(MNPV_mod5_params,SNPV_mod6_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod56_df$params <- rownames(mod56_df)

# write_csv(mod56_df,"transmission_fitting/output/mod56_params.csv")
# 
# pdf("transmission_fitting/figures/mod56_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod56_K, nu_bars = c(SNPV_mod6_params[c(1,4:8)],MNPV_mod5_params[c(1,2,5:9,15:19)]),
#              Cs = c(SNPV_mod6_params[c(2,3,9:18)],MNPV_mod5_params[c(3:4,10:14,20:24)]))
# dev.off()

fit_mod57_K <- rstan::stan(file = 'models/Mod57_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = nchains,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod57_K)$summary[c(SNPV_mod7_params,MNPV_mod5_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod57_df <- data.frame(summary(fit_mod57_K)$summary[c(MNPV_mod5_params,SNPV_mod7_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod57_df$params <- rownames(mod57_df)

# write_csv(mod57_df,"transmission_fitting/output/mod57_params.csv")
# 
# pdf("transmission_fitting/figures/mod57_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod57_K, nu_bars = c(SNPV_mod7_params[c(1,2,4:8,14:18)],MNPV_mod5_params[c(1,2,5:9,15:19)]),
#              Cs = c(SNPV_mod7_params[c(3,9:13)],MNPV_mod5_params[c(3:4,10:14,20:24)]))
# dev.off()

fit_mod65_K <- rstan::stan(file = 'models/Mod65_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = nchains,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod65_K)$summary[c(SNPV_mod5_params,MNPV_mod6_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod65_df <- data.frame(summary(fit_mod65_K)$summary[c(MNPV_mod6_params,SNPV_mod5_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod65_df$params <- rownames(mod65_df)

# write_csv(mod65_df,"transmission_fitting/output/mod65_params.csv")
# 
# pdf("transmission_fitting/figures/mod65_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod65_K, nu_bars = c(SNPV_mod5_params[c(1,2,5:9,15:19)],MNPV_mod6_params[c(1,4:8)]),
#              Cs = c(SNPV_mod5_params[c(3:4,10:14,20:24)],MNPV_mod6_params[c(2,3,9:18)]))
# dev.off()

fit_mod75_K <- rstan::stan(file = 'models/Mod75_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = nchains,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

summary(fit_mod75_K)$summary[c(SNPV_mod5_params,MNPV_mod7_params),c('mean','sd','se_mean','2.5%','97.5%')]

mod75_df <- data.frame(summary(fit_mod75_K)$summary[c(MNPV_mod7_params,SNPV_mod5_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod75_df$params <- rownames(mod75_df)

# write_csv(mod75_df,"transmission_fitting/output/mod75_params.csv")
# 
# pdf("transmission_fitting/figures/mod75_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod75_K, nu_bars = c(SNPV_mod5_params[c(1,2,5:9,15:19)],MNPV_mod7_params[c(1,2,4:8,14:18)]),
#              Cs = c(SNPV_mod5_params[c(3:4,10:14,20:24)],MNPV_mod7_params[c(3,9:13)]))
# dev.off()

fit_mod67_K <- rstan::stan(file = 'models/Mod67_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = nchains,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

mod67_df <- data.frame(summary(fit_mod67_K)$summary[c(MNPV_mod6_params,SNPV_mod7_params),
                                                    c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod67_df$params <- rownames(mo67_df)

write_csv(mod67_df,"transmission_fitting/output/mod67_params.csv")

pdf("transmission_fitting/figures/mod67_K.pdf",height = 8, width = 16)
plot_summary(mod = fit_mod67_K, nu_bars = c(SNPV_mod7_params[c(1,2,4:8,14:18)],MNPV_mod6_params[c(1,4:8)]),
             Cs = c(SNPV_mod7_params[c(3,9:13)],MNPV_mod6_params[c(2,3,9:18)]))
dev.off()


fit_mod76_K <- rstan::stan(file = 'models/Mod76_SxMxT_K.stan',
                           data = df_SxMxT_all,
                           iter = niter,
                           chains = nchains,
                           warmup = nwarmup,
                           refresh=0, control = list(adapt_delta = 0.95,
                                                     max_treedepth = 12)) 

mod76_df <- data.frame(summary(fit_mod76_K)$summary[c(MNPV_mod7_params,SNPV_mod6_params),
                                                   c('mean','sd','se_mean','25%','75%','2.5%','97.5%')])
mod76_df$params <- rownames(mod76_df)

# write_csv(mod76_df,"transmission_fitting/output/mod76_params.csv")
# 
# pdf("transmission_fitting/figures/mod76_K.pdf",height = 8, width = 16)
# plot_summary(mod = fit_mod76_K, nu_bars = c(SNPV_mod6_params[c(1,4:8)],MNPV_mod7_params[c(1,2,4:8,14:18)]),
#              Cs = c(SNPV_mod6_params[c(2,3,9:18)],MNPV_mod7_params[c(3,9:13)]))
# dev.off()

# save(fit_mod65_K, fit_mod56_K,fit_mod67_K,fit_mod76_K,
#      fit_mod63_K,fit_mod36_K,
#      fit_mod75_K,fit_mod57_K,
#      fit_mod73_K,fit_mod37_K,
#      fit_mod35_K,fit_mod53_K,
#      fit_mod77_K,fit_mod66_K,fit_mod55_K,fit_mod33_K,fit_mod0_K,fit_mod1_K, file = "output/models_C2_2021.RData")
# 
load("transmission_fitting/output/models_C2_2021_updated.RData")

elpd_diffs <- loo_compare(loo(fit_mod0_K),
            loo(fit_mod1_K),
            loo(fit_mod33_K),loo(fit_mod35_K),loo(fit_mod36_K),loo(fit_mod37_K),
            loo(fit_mod53_K),loo(fit_mod63_K),loo(fit_mod73_K),
            loo(fit_mod55_K),loo(fit_mod56_K),loo(fit_mod57_K),
            loo(fit_mod65_K),loo(fit_mod75_K),
            loo(fit_mod66_K),loo(fit_mod67_K), loo(fit_mod76_K),
            loo(fit_mod77_K))


elpd_df <- data.frame(elpd_diffs)
elpd_df$mod_num <- rownames(elpd_df)

mod_list <- data.frame(mod_names = c('mod0','mod1',
              'mod33','mod35','mod36','mod37',
              'mod53','mod63','mod73',
              'mod55','mod56','mod57',
              'mod65','mod75',
              'mod66','mod67','mod76',
              'mod77'),
              mod_num = paste0('model',1:18),
              MNPV_mod = c(0,1,
                           3,3,3,3,
                           5,6,7,
                           5,5,5,
                           6,7,
                           6,6,7,
                           7),
              SNPV_mod = c(0,1,
                           3,5,6,7,
                           3,3,3,
                           5,6,7,
                           5,5,
                           6,7,6,
                           7))


elpd_df <- merge(mod_list,elpd_df)

elpd_df <- elpd_df %>% arrange(desc(elpd_diff))

write_csv(elpd_df, "transmission_fitting/output/elpd_diff_df.csv")
elpd_df <- read_csv("transmission_fitting/output/elpd_diff_df.csv")

elpd_df_small <- elpd_df %>% select(mod_names,looic,se_looic)

###################

# mod0_all <- read_csv("transmission_fitting/output/mod0_params.csv")
# mod1_all <- read_csv("transmission_fitting/output/mod1_params.csv")
# mod33_all <- read_csv("transmission_fitting/output/mod33_params.csv")
# mod35_all <- read_csv("transmission_fitting/output/mod35_params.csv")
# mod36_all <- read_csv("transmission_fitting/output/mod36_params.csv")
# mod37_all <- read_csv("transmission_fitting/output/mod37_params.csv")
# mod53_all <- read_csv("transmission_fitting/output/mod53_params.csv")
# mod55_all <- read_csv("transmission_fitting/output/mod55_params.csv")
# mod56_all <- read_csv("transmission_fitting/output/mod56_params.csv")
# mod57_all <- read_csv("transmission_fitting/output/mod57_params.csv")
# mod63_all <- read_csv("transmission_fitting/output/mod63_params.csv")
# mod65_all <- read_csv("transmission_fitting/output/mod65_params.csv")
# mod66_all <- read_csv("transmission_fitting/output/mod66_params.csv")
# mod67_all <- read_csv("transmission_fitting/output/mod67_params.csv")
# mod73_all <- read_csv("transmission_fitting/output/mod73_params.csv")
# mod75_all <- read_csv("transmission_fitting/output/mod75_params.csv")
# mod76_all <- read_csv("transmission_fitting/output/mod76_params.csv")
# mod77_all <- read_csv("transmission_fitting/output/mod77_params.csv")
# 
# mod0_all <- mod0_all %>% mutate(mod_names = 'mod0')
# mod1_all <- mod1_all %>% mutate(mod_names = 'mod1')
# mod33_all <- mod33_all %>% mutate(mod_names = 'mod33')
# mod35_all <- mod35_all %>% mutate(mod_names = 'mod35')
# mod36_all <- mod36_all %>% mutate(mod_names = 'mod36')
# mod37_all <- mod37_all %>% mutate(mod_names = 'mod37')
# mod53_all <- mod53_all %>% mutate(mod_names = 'mod53')
# mod55_all <- mod55_all %>% mutate(mod_names = 'mod55')
# mod56_all <- mod56_all %>% mutate(mod_names = 'mod56')
# mod57_all <- mod57_all %>% mutate(mod_names = 'mod57')
# mod63_all <- mod63_all %>% mutate(mod_names = 'mod63')
# mod65_all <- mod65_all %>% mutate(mod_names = 'mod65')
# mod66_all <- mod66_all %>% mutate(mod_names = 'mod66')
# mod67_all <- mod67_all %>% mutate(mod_names = 'mod67')
# mod73_all <- mod73_all %>% mutate(mod_names = 'mod73')
# mod75_all <- mod75_all %>% mutate(mod_names = 'mod75')
# mod76_all <- mod76_all %>% mutate(mod_names = 'mod76')
# mod77_all <- mod77_all %>% mutate(mod_names = 'mod77')
# 
# mod67_all %>% filter(params %in% c("muC_M[1]",'muC_M[2]','muC_S'))
# mod33_all %>% filter(params %in% c("muC_M",'muC_S'))
# mod1_all %>% filter(params %in% c("muC"))
# 
# all_models_params <- rbind(mod0_all,mod1_all,mod33_all,mod35_all,mod36_all,mod37_all,mod53_all,mod55_all,mod56_all,mod57_all,
#       mod63_all,mod65_all,mod66_all,mod67_all,mod73_all,mod75_all,mod76_all,mod77_all)

model_desc <- data.frame(model = c(0,1,3,5,6,7), desc = c('no isolates',
                                                          'no morphotype',
                                                          'nu and C vary by morph',
                                                          'nu and C vary by morph + tree',
                                                          'nu varies by morph, C by morph + tree',
                                                          'nu varies by morph + tree, C by morph'))

#all_models_params <- merge(all_models_params,model_desc)

write_csv(all_models_params,"transmission_fitting/output/all_model_params.csv")

all_models_params <- read_csv("transmission_fitting/output/all_model_params.csv")

###################
#
# Calculating SSEs
#
###################

rho <- 0.032
xvals <- seq(0,300,0.1)

input_data_all <- input_data %>% mutate(FS = -log(1-(virus/alive))) %>% rename(morph = capsid) %>% 
  mutate(tree = ifelse(tree == "DO",'Douglas-fir','Grand fir'))

non_na_params <- all_models_params %>% filter(level == 'grand', tree %in% c("DO","GR")) %>%
  mutate(tree = ifelse(tree == "DO",'Douglas-fir','Grand fir')) #
na_params_DO <- all_models_params %>% filter(level == 'grand', is.na(tree), !is.na(morph)) %>% mutate(tree = 'Douglas-fir')
na_params_GR <- all_models_params %>% filter(level == 'grand', is.na(tree), !is.na(morph)) %>% mutate(tree = 'Grand fir')

no_morph_SD <- all_models_params %>% filter(level == 'grand', is.na(tree), is.na(morph)) %>% mutate(tree = 'Douglas-fir', morph = "SNPV")
no_morph_MD <- all_models_params %>% filter(level == 'grand', is.na(tree), is.na(morph)) %>% mutate(tree = 'Douglas-fir', morph = "MNPV")
no_morph_SG <- all_models_params %>% filter(level == 'grand', is.na(tree), is.na(morph)) %>% mutate(tree = 'Grand fir', morph = "SNPV")
no_morph_MG <- all_models_params %>% filter(level == 'grand', is.na(tree), is.na(morph)) %>% mutate(tree = 'Grand fir', morph = "MNPV")

all_params2 <- rbind(non_na_params,na_params_DO,na_params_GR,no_morph_MD,no_morph_MG,no_morph_SD,no_morph_SG)

all_params2 %>% 
  select(mean,level,morph,tree,param,mod_names) %>%  #%>% pivot_wider(names_from = 'param',values_from = 'mean') %>%
  mutate(n = 1) %>% 
  group_by(morph,tree) %>% count(n)

all_params2 <- all_params2 %>% 
  select(mean,level,morph,tree,param,mod_names) %>% pivot_wider(names_from = 'param',values_from = 'mean') 

all_params2 %>% filter(mod_names == 'mod67')

input_data_all2 <- merge(input_data_all,all_params2)
input_data_all2 <- input_data_all2 %>% mutate(pred = get_logC(P0,nu,C,7))

input_data_all2 %>% ggplot() + aes(x = P0, y = FS-pred, color = tree) + geom_point() + theme_classic() +
  facet_wrap(~mod_names)

all_sses <- input_data_all2 %>% mutate(se= (FS-pred)^2) %>% group_by(mod_names) %>% 
  summarize(sse = sum(se)) %>% arrange(sse)

elpd_df2 <- merge(elpd_df,all_sses)
elpd_df2 <- elpd_df2 %>% arrange(looic)

write_csv(elpd_df2, "output/elpd_diff_sse_df.csv")

###################
#
# PLOTTING OUTPUT
#
###################

input_data <- read_csv("transmission_fitting/data/STAN_input_data.csv")

input_data <- input_data %>% filter(strain %ni% c('H07','MEW'), alive>0) %>% 
  select(-isolate) 

input_data <- input_data %>% mutate(FS = -log(1-(virus/alive)))

input_data[input_data$tree == "DO",]$tree <- "Douglas-fir"
input_data[input_data$tree == "GR",]$tree <- "Grand fir"

l1 <-50
l2 <- 150
l3 <- 200

input_data$bin <- 0
input_data[input_data$P0 <=l1,]$bin <- 1
input_data[input_data$P0 >l1 &input_data$P0 <=l2,]$bin <- 2
input_data[input_data$P0 >l2 &input_data$P0 <= l3,]$bin <- 3
input_data[input_data$P0 >l3,]$bin <- 4

mean_vals2 <- input_data %>% mutate(FS = -log(1-(virus/alive))) %>%
  group_by(bin,capsid,tree) %>%
  summarize(mean_P0 = mean(P0,na.rm=TRUE),
            mean_FS = mean(FS,na.rm=TRUE),
            sd_P0 = sd(P0,na.rm=TRUE),
            sd_FS = sd(FS,na.rm=TRUE),
            n = length(P0)) %>%
  mutate(se_P0 = sd_P0/sqrt(n),
         se_FS = sd_FS/sqrt(n)) %>% 
  mutate(morph = ifelse(capsid == 'MNPV','Multi-capsid morphotype','Single-capsid morphotype'))

mean_vals2 %>% ggplot() + aes(x = mean_P0, y = mean_FS) +
  geom_point(color = 'red') +
  geom_errorbar(aes(ymin = mean_FS - se_FS,ymax = mean_FS + se_FS), color = 'red') +
  geom_errorbar(aes(xmin = mean_P0 - se_P0,xmax = mean_P0 + se_P0), color = 'red') +
  #geom_point(data = input_data, aes(x = P0, y = FS), color = 'grey55') + 
  theme_classic(base_size = 15) + facet_grid(tree~capsid, scales = 'free')


snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

rho <- 0.032
xvals <- seq(0,300,0.1)

all_grand <- all_models_params %>% filter(level == 'grand')

best_mod <- all_grand %>% filter(mod_names == 'mod67')

get_values <- function(mod_df){
  
  mnpv_vals <- mod_df %>% filter(morph == 'MNPV')
  snpv_vals <- mod_df %>% filter(morph == 'SNPV')
  
  mnpv_vals_C <- mnpv_vals %>% filter(param == "C")
  snpv_vals_C <- snpv_vals %>% filter(param == "C")  
  
  mnpv_vals_nu <- mnpv_vals %>% filter(param == "nu")
  snpv_vals_nu <- snpv_vals %>% filter(param == "nu")  
  
  if (is.na(mnpv_vals_C$tree[1]) == TRUE){
    mnpv_vals_C_DO_mean = mnpv_vals_C %>% pull(mean)
    mnpv_vals_C_GR_mean = mnpv_vals_C %>% pull(mean)
    
    mnpv_vals_C_DO_b1 = mnpv_vals_C %>% pull(`X2.5.`)
    mnpv_vals_C_DO_b2 = mnpv_vals_C %>% pull(`X97.5.`)
    
    mnpv_vals_C_GR_b1 = mnpv_vals_C %>% pull(`X2.5.`)
    mnpv_vals_C_GR_b2 = mnpv_vals_C %>% pull(`X97.5.`)
    
  } else{
    mnpv_vals_C_DO_mean = mnpv_vals_C %>% filter(tree == "DO") %>% pull(mean)
    mnpv_vals_C_GR_mean = mnpv_vals_C %>% filter(tree == "GR")%>% pull(mean)
    
    mnpv_vals_C_DO_b1 = mnpv_vals_C %>% filter(tree == "DO") %>% pull(`X2.5.`)
    mnpv_vals_C_DO_b2 = mnpv_vals_C %>% filter(tree == "DO") %>% pull(`X97.5.`)
    
    mnpv_vals_C_GR_b1 = mnpv_vals_C %>% filter(tree == "GR") %>% pull(`X2.5.`)
    mnpv_vals_C_GR_b2 = mnpv_vals_C %>% filter(tree == "GR") %>% pull(`X97.5.`)
  }
  
  if (is.na(mnpv_vals_nu$tree[1]) == TRUE){
    mnpv_vals_nu_DO_mean = mnpv_vals_nu %>% pull(mean)
    mnpv_vals_nu_GR_mean = mnpv_vals_nu %>% pull(mean)
    
    mnpv_vals_nu_DO_b1 = mnpv_vals_nu %>% pull(`X2.5.`)
    mnpv_vals_nu_DO_b2 = mnpv_vals_nu %>% pull(`X97.5.`)
    
    mnpv_vals_nu_GR_b1 = mnpv_vals_nu %>% pull(`X2.5.`)
    mnpv_vals_nu_GR_b2 = mnpv_vals_nu %>% pull(`X97.5.`)
    
  } else{
    mnpv_vals_nu_DO_mean = mnpv_vals_nu %>% filter(tree == "DO") %>% pull(mean)
    mnpv_vals_nu_GR_mean = mnpv_vals_nu %>% filter(tree == "GR")%>% pull(mean)
    
    mnpv_vals_nu_DO_b1 = mnpv_vals_nu %>% filter(tree == "DO") %>% pull(`X2.5.`)
    mnpv_vals_nu_DO_b2 = mnpv_vals_nu %>% filter(tree == "DO") %>% pull(`X97.5.`)
    
    mnpv_vals_nu_GR_b1 = mnpv_vals_nu %>% filter(tree == "GR") %>% pull(`X2.5.`)
    mnpv_vals_nu_GR_b2 = mnpv_vals_nu %>% filter(tree == "GR") %>% pull(`X97.5.`)
  }
  
  if (is.na(snpv_vals_C$tree[1]) == TRUE){
    snpv_vals_C_DO_mean = snpv_vals_C %>% pull(mean)
    snpv_vals_C_GR_mean = snpv_vals_C %>% pull(mean)
    
    snpv_vals_C_DO_b1 = snpv_vals_C %>% pull(`X2.5.`)
    snpv_vals_C_DO_b2 = snpv_vals_C %>% pull(`X97.5.`)
    
    snpv_vals_C_GR_b1 = snpv_vals_C %>% pull(`X2.5.`)
    snpv_vals_C_GR_b2 = snpv_vals_C %>% pull(`X97.5.`)
    
  } else{
    snpv_vals_C_DO_mean = snpv_vals_C %>% filter(tree == "DO") %>% pull(mean)
    snpv_vals_C_GR_mean = snpv_vals_C %>% filter(tree == "GR")%>% pull(mean)
    
    snpv_vals_C_DO_b1 = snpv_vals_C %>% filter(tree == "DO") %>% pull(`X2.5.`)
    snpv_vals_C_DO_b2 = snpv_vals_C %>% filter(tree == "DO") %>% pull(`X97.5.`)
    
    snpv_vals_C_GR_b1 = snpv_vals_C %>% filter(tree == "GR") %>% pull(`X2.5.`)
    snpv_vals_C_GR_b2 = snpv_vals_C %>% filter(tree == "GR") %>% pull(`X97.5.`)
  }
  
  if (is.na(snpv_vals_nu$tree[1]) == TRUE){
    snpv_vals_nu_DO_mean = snpv_vals_nu %>% pull(mean)
    snpv_vals_nu_GR_mean = snpv_vals_nu %>% pull(mean)
    
    snpv_vals_nu_DO_b1 = snpv_vals_nu %>% pull(`X2.5.`)
    snpv_vals_nu_DO_b2 = snpv_vals_nu %>% pull(`X97.5.`)
    
    snpv_vals_nu_GR_b1 = snpv_vals_nu %>% pull(`X2.5.`)
    snpv_vals_nu_GR_b2 = snpv_vals_nu %>% pull(`X97.5.`)
    
  } else{
    snpv_vals_nu_DO_mean = snpv_vals_nu %>% filter(tree == "DO") %>% pull(mean)
    snpv_vals_nu_GR_mean = snpv_vals_nu %>% filter(tree == "GR")%>% pull(mean)
    
    snpv_vals_nu_DO_b1 = snpv_vals_nu %>% filter(tree == "DO") %>% pull(`X2.5.`)
    snpv_vals_nu_DO_b2 = snpv_vals_nu %>% filter(tree == "DO") %>% pull(`X97.5.`)
    
    snpv_vals_nu_GR_b1 = snpv_vals_nu %>% filter(tree == "GR") %>% pull(`X2.5.`)
    snpv_vals_nu_GR_b2 = snpv_vals_nu %>% filter(tree == "GR") %>% pull(`X97.5.`)
  }
  
  snpv_DO_1 <- data.frame(P0 = xvals,
                          FS = get_logC(xvals,nu_bar = snpv_vals_nu_DO_mean, C = snpv_vals_C_DO_mean, time = 7),
                          tree = 'Douglas-fir', morph = "Single-capsid morphotype",
                          q1 = get_logC(xvals,nu_bar = snpv_vals_nu_DO_b1, C = snpv_vals_C_DO_b2, time = 7),
                          q2 = get_logC(xvals,nu_bar = snpv_vals_nu_DO_b2, C = snpv_vals_C_DO_b1, time = 7))
  
  snpv_GR_1 <- data.frame(P0 = xvals,
                          FS = get_logC(xvals,nu_bar = snpv_vals_nu_GR_mean, C = snpv_vals_C_GR_mean, time = 7),
                          tree = 'Grand fir', morph = "Single-capsid morphotype",
                          q1 = get_logC(xvals,nu_bar = snpv_vals_nu_GR_b1, C = snpv_vals_C_GR_b2, time = 7),
                          q2 = get_logC(xvals,nu_bar = snpv_vals_nu_GR_b2, C = snpv_vals_C_GR_b1, time = 7))
  
  mnpv_DO_1 <- data.frame(P0 = xvals,
                          FS = get_logC(xvals,nu_bar = mnpv_vals_nu_DO_mean, C = mnpv_vals_C_DO_mean, time = 7),
                          tree = 'Douglas-fir', morph = "Multi-capsid morphotype",
                          q1 = get_logC(xvals,nu_bar = mnpv_vals_nu_DO_b1, C = mnpv_vals_C_DO_b2, time = 7),
                          q2 = get_logC(xvals,nu_bar = mnpv_vals_nu_DO_b2, C = mnpv_vals_C_DO_b1, time = 7))
  
  mnpv_GR_1 <- data.frame(P0 = xvals,
                          FS = get_logC(xvals,nu_bar = mnpv_vals_nu_GR_mean, C = mnpv_vals_C_GR_mean, time = 7),
                          tree = 'Grand fir', morph = "Multi-capsid morphotype",
                          q1 = get_logC(xvals,nu_bar = mnpv_vals_nu_GR_b1, C = mnpv_vals_C_GR_b2, time = 7),
                          q2 = get_logC(xvals,nu_bar = mnpv_vals_nu_GR_b2, C = mnpv_vals_C_GR_b1, time = 7))

  all_lines <- rbind(snpv_DO_1,snpv_GR_1,
        mnpv_DO_1,mnpv_GR_1)
    
  
  return(all_lines)
}

mod67_lines <- get_values(best_mod)

pdf("figures/transmission_best_model.pdf",height = 6, width = 10)
mod67_lines %>% ggplot() + geom_ribbon(aes(x = P0,ymin = q1, ymax = q2, fill = tree), alpha = 0.4) + 
  geom_line(aes(x = P0, y = FS, color = tree)) + 
  facet_grid(morph~tree) + 
  scale_color_manual("", values = c('Douglas-fir' = 'green3', 'Grand fir' = 'darkgreen'))+
  scale_fill_manual("", values = c('Douglas-fir' = 'green3', 'Grand fir' = 'darkgreen')) +
  geom_point(data = mean_vals2, aes(x = mean_P0, y = mean_FS, color = tree), size = 3) + 
  geom_errorbar(data = mean_vals2, aes(x = mean_P0, ymin = mean_FS - se_FS,ymax = mean_FS + se_FS, color = tree), size = 1.5) +
  geom_errorbar(data = mean_vals2, aes(y = mean_FS, xmin = mean_P0 - se_P0,xmax = mean_P0 + se_P0, color = tree), size = 1.5) + 
  theme_bw(base_size = 15) + 
  ylab(expression(paste("-ln",bgroup("(", frac(S[t],S[0]),")"))))

dev.off()
#########################
#########################
#########################
#########################
#########################

iso_params <- all_models_params %>% filter(mod_names == 'mod67',level == 'isolate')

iso_lines <- c()
for(i in 1:5){
  MNPV_nu_DO <- iso_params %>% filter(morph == "MNPV", isolate == i, param == 'nu') 
  MNPV_nu_GR <- iso_params %>% filter(morph == "MNPV", isolate == i, param == 'nu') 
  
  MNPV_C_DO <- iso_params %>% filter(morph == "MNPV", isolate == i, param == 'C', tree == "DO") 
  MNPV_C_GR <- iso_params %>% filter(morph == "MNPV", isolate == i, param == 'C', tree == "GR") 
  
  MNPV_FS_DO <- get_logC(xvals,nu_bar = MNPV_nu_DO$mean, C = MNPV_C_DO$mean, time = 7)
  MNPV_b1_DO <- get_logC(xvals,nu_bar = MNPV_nu_DO$`X2.5.`, C = MNPV_C_DO$`X97.5.`, time = 7)
  MNPV_b2_DO <- get_logC(xvals,nu_bar = MNPV_nu_DO$`X97.5.`, C = MNPV_C_DO$`X2.5.`, time = 7)
  
  MNPV_FS_GR <- get_logC(xvals,nu_bar = MNPV_nu_GR$mean, C = MNPV_C_GR$mean, time = 7)
  MNPV_b1_GR <- get_logC(xvals,nu_bar = MNPV_nu_GR$`X2.5.`, C = MNPV_C_GR$`X97.5.`, time = 7)
  MNPV_b2_GR <- get_logC(xvals,nu_bar = MNPV_nu_GR$`X97.5.`, C = MNPV_C_GR$`X2.5.`, time = 7)
  
  mnpv_df_DO <- data.frame(morph ="Multi-capsid morphotype",isolate = i, tree = 'DO', FS = MNPV_FS_DO, b1 = MNPV_b1_DO, b2 = MNPV_b2_DO, P0 = xvals)
  mnpv_df_GR <- data.frame(morph ="Multi-capsid morphotype",isolate = i, tree = 'GR', FS = MNPV_FS_GR, b1 = MNPV_b1_GR, b2 = MNPV_b2_GR, P0 = xvals)
  
  
  SNPV_C_DO <- iso_params %>% filter(morph == "SNPV", isolate == i + 5, param == 'C') 
  SNPV_C_GR <- iso_params %>% filter(morph == "SNPV", isolate == i + 5, param == 'C') 
  
  SNPV_nu_DO <- iso_params %>% filter(morph == "SNPV", isolate == i + 5, param == 'nu', tree == "DO") 
  SNPV_nu_GR <- iso_params %>% filter(morph == "SNPV", isolate == i + 5, param == 'nu', tree == "GR") 
  
  SNPV_FS_DO <- get_logC(xvals,nu_bar = SNPV_nu_DO$mean, C = SNPV_C_DO$mean, time = 7)
  SNPV_b1_DO <- get_logC(xvals,nu_bar = SNPV_nu_DO$`X2.5.`, C = SNPV_C_DO$`X97.5.`, time = 7)
  SNPV_b2_DO <- get_logC(xvals,nu_bar = SNPV_nu_DO$`X97.5.`, C = SNPV_C_DO$`X2.5.`, time = 7)
  
  SNPV_FS_GR <- get_logC(xvals,nu_bar = SNPV_nu_GR$mean, C = SNPV_C_GR$mean, time = 7)
  SNPV_b1_GR <- get_logC(xvals,nu_bar = SNPV_nu_GR$`X2.5.`, C = SNPV_C_GR$`X97.5.`, time = 7)
  SNPV_b2_GR <- get_logC(xvals,nu_bar = SNPV_nu_GR$`X97.5.`, C = SNPV_C_GR$`X2.5.`, time = 7)
  
  snpv_df_DO <- data.frame(morph ="Single-capsid morphotype",isolate = i + 5, tree = 'DO', FS = SNPV_FS_DO, b1 = SNPV_b1_DO, b2 = SNPV_b2_DO, P0 = xvals)
  snpv_df_GR <- data.frame(morph ="Single-capsid morphotype",isolate = i + 5, tree = 'GR', FS = SNPV_FS_GR, b1 = SNPV_b1_GR, b2 = SNPV_b2_GR, P0 = xvals)
  
  temp <- rbind(snpv_df_DO,snpv_df_GR, mnpv_df_GR,mnpv_df_DO)
  iso_lines <- rbind(iso_lines,temp)
  
}

iso_lines[iso_lines$tree == "DO",]$tree <- "Douglas-fir"
iso_lines[iso_lines$tree == "GR",]$tree <- "Grand fir"

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'
iso_info <- input_data2 %>% group_by(strain,isolate) %>% count(isolate) %>% select(-n) %>% arrange(isolate)

iso_lines <- merge(iso_lines,iso_info)
iso_lines$strain <- factor(iso_lines$strain,levels = iso_info$strain)


iso_data <- input_data
iso_data <- iso_data %>% rename(morph = capsid)
iso_data[iso_data$tree == "DO",]$tree <- "Douglas-fir"
iso_data[iso_data$tree == "GR",]$tree <- "Grand fir"
iso_data[iso_data$morph == "SNPV",]$morph <- "Single-capsid morphotype"
iso_data[iso_data$morph == "MNPV",]$morph <- "Multi-capsid morphotype"

iso_data$bin <- 0
iso_data[iso_data$P0 <=100,]$bin <- 1
iso_data[iso_data$P0 >100,]$bin <- 2

iso_mean <- iso_data %>% mutate(FS = -log(1-(virus/alive)),
                                FI = virus/alive) %>%
  group_by(bin,morph,strain,tree) %>%
  summarize(mean_P0 = mean(P0,na.rm=TRUE),
            mean_FS = mean(FS,na.rm=TRUE),
            sd_P0 = sd(P0,na.rm=TRUE),
            sd_FS = sd(FS,na.rm=TRUE),
            mean_FI = mean(FI,na.rm=TRUE),
            sd_FI = sd(FI, na.rm = TRUE),
            n = length(P0)) %>%
  mutate(se_P0 = sd_P0/sqrt(n),
         se_FS = sd_FS/sqrt(n))

iso_mean$strain <- factor(iso_mean$strain,levels = iso_info$strain)

pdf("figures/iso_lines95.pdf",height = 8, width = 16)
iso_lines %>% ggplot() + geom_ribbon(aes(x = P0,ymin = b1, ymax = b2, fill = morph), alpha = 0.4) + 
  geom_line(aes(x = P0, y = FS, color = morph))+
  theme_classic(base_size = 15) + 
  facet_grid(tree~strain) +
  scale_color_manual("", values = c('Single-capsid morphotype' = snpv_col, 'Multi-capsid morphotype' = mnpv_col))  +
  scale_fill_manual("", values = c('Single-capsid morphotype' = snpv_col, 'Multi-capsid morphotype' = mnpv_col))  +
  geom_point(data = iso_mean,aes(x = mean_P0, y = mean_FS, color = morph), size = 2) +
  geom_errorbar(data = iso_mean, aes(x = mean_P0, ymin = mean_FS - se_FS,ymax = mean_FS + se_FS, color = morph)) +
  geom_errorbar(data = iso_mean, aes(y = mean_FS, xmin = mean_P0 - se_P0,xmax = mean_P0 + se_P0, color = morph)) +
  theme_bw(base_size = 15) + 
  theme(legend.position = 'top') +
  ylab(expression(paste("-ln",bgroup("(", frac(S[t],S[0]),")")))) +
  xlab(expression(paste("Initial pathogen density- cadavers per m^2"))) +
  scale_x_continuous(breaks = c(0,150,300))
dev.off()


