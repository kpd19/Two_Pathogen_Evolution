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

load("data/models_C_2021_rho_Apr6.RData")

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

draws_mat_M6 <- as_draws_matrix(fit_mod6_MNPV_K)

draws_sub_M6 <- subset_draws(draws_mat_M6, c('mu','muC[1]',"muC[2]"))
aws_sub_M6 <- rename_variables(draws_sub_M6, C_DO = `muC[1]`, C_GR = `muC[2]`)

summarise_draws(draws_sub_M6)

thin_M6 <- thin_draws(draws_sub_M6, thin = 150)
M6_thin_df <- data.frame(thin_M6)

M6_thin_df <- M6_thin_df %>% rename(nu_DO = mu) %>% mutate(nu_GR = nu_DO) %>% mutate(morph = "MNPV", draw = 1:225)

draws_mat_S7 <- as_draws_matrix(fit_mod7_SNPV_K)

draws_sub_S7 <- subset_draws(draws_mat_S7, c('mu[1]','mu[2]','muC'))
draws_sub_S7 <- rename_variables(draws_sub_S7, nu_DO = `mu[1]`, nu_GR = `mu[2]`)

summarise_draws(draws_sub_S7)

thin_S7 <- thin_draws(draws_sub_S7, thin = 150)
S7_thin_df <- data.frame(thin_S7)

S7_thin_df <- S7_thin_df %>% rename(C_DO = muC) %>% mutate(C_GR = C_DO) %>% mutate(morph = "SNPV", draw = 1:225)

all_thinned <- rbind(M6_thin_df,S7_thin_df)

write_csv(all_thinned,"data/thinned_posteriors.csv")

##############
#
# Analyze output
#
##############

all_thinned <- read_csv("data/thinned_posteriors.csv")

thinned_results <- read_csv("data/posteriors_df.csv")
dens <- read_csv("_data/test_df.csv")
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

write_csv(thinned_results, "data/thinned_posterior_FI.csv")

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

  
  