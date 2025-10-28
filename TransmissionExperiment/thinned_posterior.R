library(tidyverse)
require(rstan)
require(loo)
library(bayesplot)
library(gridExtra)
library(posterior)

get_logC <- function(x,nu_bar,C,time){
  
  log_val = (1/C^2)*log(1 + nu_bar*C^2*time*x*rho)
  
  return(log_val)
}

load("transmission_fitting/output/best_model_2021.RData")

input_data <- read_csv("transmission_fitting/data/STAN_input_data.csv")

rho <- 0.032
xvals <- seq(0,300,0.1)

input_data$bin <- 0
input_data[input_data$P0 <=100,]$bin <- 1
input_data[input_data$P0 >100 &input_data$P0 <=200,]$bin <- 2
input_data[input_data$P0 >200,]$bin <- 3

mean_vals2 <- input_data %>% mutate(FS = -log(1-(virus/alive))) %>%
  group_by(bin,capsid,tree) %>%
  summarize(mean_P0 = mean(P0,na.rm=TRUE),
            mean_FS = mean(FS,na.rm=TRUE),
            sd_P0 = sd(P0,na.rm=TRUE),
            sd_FS = sd(FS,na.rm=TRUE),
            n = length(P0)) %>%
  mutate(se_P0 = sd_P0/sqrt(n),
         se_FS = sd_FS/sqrt(n))
mean_vals2[mean_vals2$tree == "DO",]$tree <- "Douglas-fir"
mean_vals2[mean_vals2$tree == "GR",]$tree <- "Grand fir"

draws_mat <- as_draws_matrix(fit_mod67_K)

draws_sub_M6 <- subset_draws(draws_mat, c('mu_M','muC_M[1]',"muC_M[2]"))
draws_sub_M6 <- rename_variables(draws_sub_M6, C_DO = `muC_M[1]`, C_GR = `muC_M[2]`)

summarise_draws(draws_sub_M6)

thin_M6 <- thin_draws(draws_sub_M6, thin = 335)
M6_thin_df <- data.frame(thin_M6)

M6_thin_df <- M6_thin_df %>% rename(nu_DO = mu_M) %>% mutate(nu_GR = nu_DO) %>% mutate(morph = "MNPV", draw = 1:225)

draws_sub_S7 <- subset_draws(draws_mat, c('mu_S[1]','mu_S[2]','muC_S'))
draws_sub_S7 <- rename_variables(draws_sub_S7, nu_DO = `mu_S[1]`, nu_GR = `mu_S[2]`)

summarise_draws(draws_sub_S7)

thin_S7 <- thin_draws(draws_sub_S7, thin = 335)
S7_thin_df <- data.frame(thin_S7)

S7_thin_df <- S7_thin_df %>% rename(C_DO = muC_S) %>% mutate(C_GR = C_DO) %>% mutate(morph = "SNPV", draw = 1:225)

all_thinned <- rbind(M6_thin_df,S7_thin_df)

write_csv(all_thinned,"data/thinned_posteriors.csv")

################
#
# Analyze output
#
################

all_thinned <- read_csv("transmission_fitting/data/thinned_posteriors.csv")

thinned_results <- read_csv("transmission_fitting/data/thinned_posterior_simulations.csv")
dens <- read_csv("transmission_fitting/data/mean_simulations.csv")
dens <- dens %>% mutate(tree_sp = ifelse(tree == "DO","Douglas-fir",'Grand fir')) %>% 
  mutate(morph_name = ifelse(morph == "SNPV",'Single-capsid morphotype','Multi-capsid morphotype')) %>% filter(host_dens <= 1000)

thinned_results <- thinned_results %>% mutate(tree_sp = ifelse(tree == "DO","Douglas-fir",'Grand fir'))%>% 
  mutate(morph_name = ifelse(morph == "SNPV",'Single-capsid morphotype','Multi-capsid morphotype')) 

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

pdf("figures/posterior_draws_hostdens.pdf",height = 6, width = 10)
thinned_results %>% ggplot() + aes(x = host_dens, y = frac_I, group = draw) + geom_line(color = 'grey75') +
  theme_classic() +facet_grid(morph~tree_sp) +
  ylab("Fraction Infected") + xlab("Initial Host Density") +
  geom_line(data = dens, aes(x =host_dens, y = frac_I, color = morph_name), group = NA, size = 2) +
  scale_color_manual("", values = c('Single-capsid morphotype' = snpv_col, 'Multi-capsid morphotype' = mnpv_col)) +
  coord_cartesian(xlim = c(0,1000)) 
dev.off()

thin_quant <- thinned_results %>% group_by(host_dens,tree,tree_sp,morph_name) %>% 
  summarize(q05 = quantile(frac_I,probs = c(0.25)),
            q95 = quantile(frac_I,probs = c(0.75)),
            mean = mean(frac_I))


plt_ribbon <- thin_quant %>% ggplot() + aes(x = host_dens, y = mean, ymin = q05, ymax = q95, fill = morph_name) +
  geom_ribbon(color = NA,alpha = 0.25) +
  theme_classic() +facet_grid(~tree_sp) +
  ylab("Fraction Infected") + xlab("Initial Host Density") +
  coord_cartesian(xlim = c(0,1000)) +
  scale_fill_manual("", values = c('Single-capsid morphotype' = snpv_col, 'Multi-capsid morphotype' = mnpv_col))

pdf('figures/posterior_envelopes50.pdf',height = 6, width = 10)
plt_ribbon + geom_line(data = dens, aes(x =host_dens, y = frac_I, color = morph_name, group = morph_name),
                       size = 2, inherit.aes =FALSE) +
  scale_color_manual("", values = c('Single-capsid morphotype' = snpv_col, 'Multi-capsid morphotype' = mnpv_col))
dev.off()

################
################
################

all_params <- read_csv("transmission_fitting/output/all_model_params.csv")

snpv_best_nu <- all_params %>% filter(morph == "SNPV", mod_names =='mod67', level == 'grand', param == 'nu')
mnpv_best_nu <- all_params %>% filter(morph == "MNPV", mod_names =='mod67', level == 'grand', param == 'nu')
mnpv_best_nu_1 <- mnpv_best_nu %>% mutate(tree = "DO")
mnpv_best_nu_2 <- mnpv_best_nu %>% mutate(tree = "GR")

snpv_best_C <- all_params %>% filter(morph == "SNPV", mod_names =='mod67', level == 'grand', param == 'C')
mnpv_best_C <- all_params %>% filter(morph == "MNPV", mod_names =='mod67', level == 'grand', param == 'C')
snpv_best_C_1 <- snpv_best_C %>% mutate(tree = "DO")
snpv_best_C_2 <- snpv_best_C %>% mutate(tree = "GR")

best_param <- rbind(snpv_best_nu, mnpv_best_nu_1,mnpv_best_nu_2, mnpv_best_C, snpv_best_C_1, snpv_best_C_2)

####################
####################
####################
####################


plt_ribbon <- thin_quant %>% ggplot() + aes(x = host_dens, y = mean, ymin = q05, ymax = q95, fill = morph_name) +
  geom_ribbon(color = NA,alpha = 0.25) +
  theme_classic(base_size = 15) +facet_grid(~tree_sp) +
  ylab("Fraction Infected") + xlab("Initial Host Density") +
  coord_cartesian(xlim = c(0,1000)) +
  scale_fill_manual("", values = c('Single-capsid morphotype' = snpv_col, 'Multi-capsid morphotype' = mnpv_col))+
  theme(legend.position = 'none')

plt1 <- plt_ribbon + geom_line(data = dens, aes(x =host_dens, y = frac_I, color = morph_name, group = morph_name),
                               size = 2, inherit.aes =FALSE) +
  scale_color_manual("", values = c('Single-capsid morphotype' = snpv_col, 'Multi-capsid morphotype' = mnpv_col)) 


times_val = 350

plt2 <- best_param %>% 
  mutate(tree_no = ifelse(tree == "DO",1,2)) %>% 
  mutate(param_no = ifelse(param == "C",1,2)) %>% 
  mutate(tree_name = ifelse(tree == "DO",'Douglas-fir','Grand fir')) %>%
  mutate(morph_name = ifelse(morph == "SNPV",'Single-capsid morphotype','Multi-capsid morphotype')) %>%
  mutate(mean = ifelse(param == "C",mean,mean*times_val)) %>%
  mutate(`X75.` = ifelse(param == "C",`X75.`,`X75.` *times_val)) %>% 
  mutate(`X25.` = ifelse(param == "C",`X25.`,`X25.` *times_val)) %>% 
  mutate(nudge_no = ifelse(morph == 'MNPV',-0.1,0.1)) %>% 
  ggplot() + aes(x = param_no + nudge_no, y = mean, color = morph_name, group = interaction(morph,tree)) +
  geom_point(size = 4,alpha=1) +
  geom_errorbar(aes(ymin = `X75.`, ymax = `X25.`), width = 0.2,size = 1.5) +
  theme_classic(base_size = 15) + 
  scale_color_manual("", values = c('Single-capsid morphotype' = snpv_col, 'Multi-capsid morphotype' = mnpv_col)) +
  scale_x_continuous(breaks = c(1,2), labels = c("C", expression(bar(nu))), limits = c(0.5,2.5)) +
  #scale_y_continuous(brekas = c(0.005,0.001,0.015, limits = c()))
  theme(axis.title.x = element_blank(),
        legend.position = c(0.15,0.95),
        legend.background = element_rect(fill = NA)) +
  ylab(expression(bar(nu))) + 
  facet_wrap(~tree_name) +
  scale_y_continuous("C estimate", sec.axis = sec_axis(~./times_val, name=c(expression(bar(nu)~'estimate'))),
                     limits = c(0.75,5)) 

pdf("transmission_fitting/figures/density_params.pdf", height = 8, width = 9.5)
grid.arrange(plt2,plt1)
dev.off()

  