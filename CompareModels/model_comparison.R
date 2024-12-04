library(tidyverse)
library(gridExtra)

ll_data <- read_csv("data/morphotype_dist_data.csv")

t5s_all <- read_csv("data/t5s_stoch_all.csv")
nt4s_all <- read_csv( "data/nt4s_fix2_all.csv")
t6s_all <- read_csv( "data/t6s_stoch_all.csv")
ht6s_all <- read_csv( "data/ht6s_stoch_all.csv")
t4s_all <- read_csv( "data/t4s_fix2_all.csv")
t6s_nomorph <- read_csv("data/t6s_noC_all.csv")
t4s_noevo <- read_csv("data/t4s_noevo_stoch_all.csv")

nt4s_all$C <- "M"
t6s_all$C <- "M + T"
ht6s_all$C <- "M"
t4s_all$C <- "M + T"
t4s_noevo$C <- "M + T"
t6s_nomorph$C <- "-"

nt4s_all$s <- "M"
t6s_all$s <- "M + T"
ht6s_all$s <- "M + T"
t4s_all$s <- "M"
t4s_noevo$s <- "-"
t6s_nomorph$s <- "M + T"

t6s_all$p <- 6
ht6s_all$p <- 6
nt4s_all$p <- 4
t4s_all$p <- 4
t4s_noevo$p <- 4
t6s_nomorph$p <- 6

t6s_all$mod_name <- "Selection mosaic for C and s"
ht6s_all$mod_name <- "Selection mosaic for s"
nt4s_all$mod_name <- "No selection mosaic"
t4s_all$mod_name <- "Selection mosaic for C"
t4s_noevo$mod_name <- "Selection mosaic for C, no evolution"
t6s_nomorph$mod_name <- "Selection mosaic for s, C does not vary"

t6s_nomorph <- t6s_nomorph %>% arrange(desc(ll))

sht6 <- ht6s_all %>% select(ll,sigma,pset,C,s,p,mod_name) %>% drop_na(ll)
st6 <- t6s_all %>% select(ll,sigma,pset,C,s,p,mod_name) %>% drop_na(ll)
snt4 <- nt4s_all %>% select(ll,sigma,pset,C,s,p,mod_name) %>% drop_na(ll)
st4 <- t4s_all %>% select(ll,sigma,pset,C,s,p,mod_name) %>% drop_na(ll)
st4ne <- t4s_noevo %>% select(ll,sigma,pset,C,s,p,mod_name) %>% drop_na(ll)
nm6 <- t6s_nomorph %>% select(ll,sigma,pset,C,s,p,mod_name) %>% drop_na(ll)

best_t6s <- st6 %>% arrange(desc(ll)) %>% head(1)
best_ht6s <- sht6 %>% arrange(desc(ll)) %>% head(1)
best_nt4s <- snt4 %>% arrange(desc(ll)) %>% head(1)
best_t4s <- st4 %>% arrange(desc(ll)) %>% head(1)
best_t4ne <- st4ne %>% arrange(desc(ll)) %>% head(1)
best_nm6 <- nm6 %>% arrange(desc(ll)) %>% head(1)

stoch_all <- rbind(snt4,sht6,st6,st4,best_t4ne,nm6)

best_mods2 <- stoch_all %>% arrange(desc(ll)) %>% group_by(mod_name,C,s) %>% filter(ll == max(ll,na.rm=TRUE)) 

aic_table3 <- best_mods2 %>% mutate(AIC = 2*p - 2*ll) %>% ungroup() %>% mutate(dAIC = AIC - min(AIC)) %>% select(-pset)

aic_table3 <- aic_table3 %>% mutate(ll = round(ll,2), AIC = round(AIC,2), dAIC = round(dAIC,2)) %>% rename(mod = mod_name)

write_csv(aic_table3,"data/aic_table.csv")

aic_table3 <- read_csv("data/aic_table.csv")

#################
#
# Percent MNPV
#
#################

t5s_pmnpv <- read_csv("data/pmnpv_t5s_fixnoext.csv")
nt4s_pmnpv <- read_csv("data/pmnpv_nt4s_fix2noext.csv")
ht6s_pmnpv <- read_csv("data/pmnpv_ht6s_fix2noext.csv")
t6s_pmnpv <- read_csv("data/pmnpv_t6s_fix2noext.csv")
t4s_pmnpv <- read_csv("data/pmnpv_t4s_fixnoext.csv")
t6s_nm_pmnpv <- read_csv("data/pmnpv_t6s_nomorph.csv")
t4s_pmnpv_ne <- read_csv("data/pmnpv_t4s_noevo.csv")

ht6s_pmnpv$C <- "M"
t6s_pmnpv$C <- "M + T"
nt4s_pmnpv$C <- "M"
t4s_pmnpv$C <- "M + T"
t4s_pmnpv_ne$C <- "M + T"
t6s_nm_pmnpv$C <- "-"

ht6s_pmnpv$s <- "M + T"
t6s_pmnpv$s <- "M + T"
nt4s_pmnpv$s <- "M"
t4s_pmnpv$s <- "M"
t4s_pmnpv_ne$s <- "-"
t6s_nm_pmnpv$s <- "M + T"

t6s_pmnpv <- t6s_pmnpv %>% select(pdoug,sigma,mean_MNPV,sd,q1,q9,pset,C,s)
ht6s_pmnpv <- ht6s_pmnpv %>% select(pdoug,sigma,mean_MNPV,sd,q1,q9,pset,C,s)
nt4s_pmnpv <- nt4s_pmnpv %>% select(pdoug,sigma,mean_MNPV,sd,q1,q9,pset,C,s)
t4s_pmnpv <- t4s_pmnpv %>% select(pdoug,sigma,mean_MNPV,sd,q1,q9,pset,C,s)
t4s_pmnpv_ne <- t4s_pmnpv_ne %>% select(pdoug,sigma,mean_MNPV,sd,q1,q9,pset,C,s)
t6s_nm_pmnpv <- t6s_nm_pmnpv %>% select(pdoug,sigma,mean_MNPV,sd,q1,q9,pset,C,s)

field_data_round <- ll_data %>% mutate(round_pd = round(Douglas_fir*9)/9) %>% group_by(round_pd) %>% 
  summarize(sum_MNPV = sum(MNPV), sum_tot = sum(total),
            pMNPV = mean(MNPV/total),
            sd_pMNPV = sd(MNPV/total),
            len = length(MNPV),
            pdoug = mean(Douglas_fir),
            sd_pd = sd(Douglas_fir),
            trees= mean(n_trees),
            sd_tree = sd(n_trees)) %>% 
  mutate(se_pMNPV = sd_pMNPV/sqrt(len),
         se_pd = sd_pd/sqrt(len))

size_pick <- 6

plt1 <- t6s_pmnpv %>% filter(sigma == 0.5,pset == 7) %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  geom_text(data = best_t6s,aes(x = 20, y = 100, label = paste0("Log-likelihod = ", round(ll,1))), size = size_pick) +
  #geom_text(data = subset(stoch_df,sigma ==0.5),aes(x = 20, y = 70, label = round(ll,1))) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("Selection mosaic for C and s") + 
  ylim(0,100)

plt6 <- t6s_nm_pmnpv %>% filter(sigma == 0.5,pset == 2) %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  geom_text(data = best_nm6,aes(x = 20, y = 100, label = paste0("Log-likelihod = ", round(ll,1))), size = size_pick) +
  #geom_text(data = subset(stoch_df,sigma ==0.5),aes(x = 20, y = 70, label = round(ll,1))) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("Selection mosaic for s, C does not vary") + 
  ylim(0,100)

plt2 <- ht6s_pmnpv %>% filter(sigma == 0.5,pset == 2) %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  geom_text(data = best_ht6s,aes(x = 20, y = 100, label = paste0("Log-likelihod = ", round(ll,1))), size = size_pick) +
  #geom_text(data = subset(stoch_df,sigma ==0.5),aes(x = 20, y = 70, label = round(ll,1))) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("Selection mosaic for s") + 
  ylim(0,100)

plt3 <- t4s_pmnpv %>% filter(sigma == 0.5,pset == 3) %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  geom_text(data = best_t4s,aes(x = 20, y = 100, label = paste0("Log-likelihod = ", round(ll,1))), size = size_pick) +
  #geom_text(data = subset(stoch_df,sigma ==0.5),aes(x = 20, y = 70, label = round(ll,1))) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("Selection mosaic for C") + 
  ylim(0,100)
nt_best <- nt4s_pmnpv %>% filter(sigma == 0.5, pset == 16) %>% select(-pdoug)
nt_best <- merge(nt_best, data.frame(pdoug = 1:36))

plt4 <- nt_best %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  geom_text(data = best_nt4s,aes(x = 20, y = 100, label = paste0("Log-likelihod = ", round(ll,1))), size = size_pick) +
  #geom_text(data = subset(stoch_df,sigma ==0.5),aes(x = 20, y = 70, label = round(ll,1))) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("No selection mosaic") + 
  ylim(0,100)

plt5 <- t4s_pmnpv_ne %>% filter(sigma == 0.5,pset == 2) %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  geom_text(data = best_t4ne,aes(x = 20, y = 100, label = paste0("Log-likelihod = ", round(ll,1))), size = size_pick) +
  #geom_text(data = subset(stoch_df,sigma ==0.5),aes(x = 20, y = 70, label = round(ll,1))) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("Selection mosaic for C, no evolution") + 
  ylim(0,100)

pdf("figures/pmnpv_top_models.pdf",height = 10, width = 20)
grid.arrange(plt1,plt6,plt2,plt4,plt3,plt5,nrow = 2)
#grid.arrange(plt1,plt3,nrow = 1)
dev.off()

pdf("figures/two_mods.pdf",height = 5, width = 12)
grid.arrange(plt1,plt4,nrow = 1)
#grid.arrange(plt1,plt3,nrow = 1)
dev.off()

num_pick <- 5

t6_psets <- t6s_all %>% select(sS_DO,sS_GR, sM_DO, sM_GR, phi, rho, ll,sigma,pset,C,s,p) %>% drop_na(ll) %>% arrange(desc(ll)) %>% head(num_pick) %>% mutate(id = 1:num_pick) 
ht6_psets <- ht6s_all %>% select(sS_DO,sS_GR, sM_DO, phi, rho, ll,sigma,pset,C,s,p) %>% drop_na(ll) %>% arrange(desc(ll)) %>% head(num_pick) %>% mutate(id = 1:num_pick) 
nt4_psets <- nt4s_all %>% select(sS,sM, phi, rho, ll,sigma,pset,C,s,p) %>% drop_na(ll) %>% arrange(desc(ll)) %>% head(num_pick) %>% mutate(id = 1:num_pick)
t4_psets <- t4s_all %>% select(sS_DO,sM_DO, phi, rho, ll,sigma,pset,C,s,p) %>% drop_na(ll) %>% arrange(desc(ll)) %>% head(num_pick) %>% mutate(id = 1:num_pick)
t4ne_psets <- t4s_noevo %>% select(phiS,phiM,rho,lambda,ll,sigma,pset) %>% drop_na(ll) %>% arrange(desc(ll)) %>% head(num_pick) %>% mutate(id = 1:num_pick) 
t6nm_psets <- t6s_nomorph %>% select(sS_DO,sS_GR, sM_DO, sM_GR, phi, rho,ll,sigma,pset) %>% drop_na(ll) %>% arrange(desc(ll)) %>% head(num_pick) %>% mutate(id = 1:num_pick) 

ht6_int <- c()

for(i in 1:length(ht6_psets$pset)){
  temp <- ht6s_pmnpv %>% filter(pset == ht6_psets$pset[i], sigma == ht6_psets$sigma[i])
  temp$id <- ht6_psets$id[i]

  ht6_int <- rbind(ht6_int,temp)
}

t6_int <- c()

for(i in 1:length(t6_psets$pset)){
  temp <- t6s_pmnpv %>% filter(pset == t6_psets$pset[i], sigma == t6_psets$sigma[i])
  temp$id <- t6_psets$id[i]
  
  t6_int <- rbind(t6_int,temp)
}



nt4_int <- c()

for(i in 1:length(nt4_psets$pset)){
  temp <- nt4s_pmnpv %>% filter(pset == nt4_psets$pset[i], sigma == nt4_psets$sigma[i])
  temp$id <- nt4_psets$id[i]
  
  nt4_int <- rbind(nt4_int,temp)
}


nt4_int <- nt4_int %>% select(-pdoug)
nt4_int <- merge(nt4_int, data.frame(pdoug = 1:35))

t4_int <- c()

for(i in 1:length(t4_psets$pset)){
  temp <- t4s_pmnpv %>% filter(pset == t4_psets$pset[i], sigma == t4_psets$sigma[i])
  temp$id <- t4_psets$id[i]
  
  t4_int <- rbind(t4_int,temp)
}



t4ne_int <- c()

for(i in 1:length(t4ne_psets$pset)){
  temp <- t4s_pmnpv_ne %>% filter(pset == t4ne_psets$pset[i], sigma == t4ne_psets$sigma[i])
  temp$id <- t4ne_psets$id[i]
  
  t4ne_int <- rbind(t4ne_int,temp)
}

t6nm_int <- c()

for(i in 1:length(t6nm_psets$pset)){
  temp <- t6s_nm_pmnpv %>% filter(pset == t6nm_psets$pset[i], sigma == t6nm_psets$sigma[i])
  temp$id <- t6nm_psets$id[i]
  
  t6nm_int <- rbind(t6nm_int,temp)
}


plt1 <- ht6_int %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  geom_text(data = ht6_psets, aes(x = 20, y = 90, label = paste0("LL = ",round(ll,1)))) +
  #geom_text(data = ht6_psets, aes(x = 20, y = 80, label = paste0('sigma = ', sigma))) +
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("C: Selection mosaic for s") + 
  ylim(0,100) + 
  facet_wrap(~id,nrow = 1) + 
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

plt2 <- t6_int %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  geom_text(data = t6_psets, aes(x = 20, y = 90, label = paste0("LL = ",round(ll,1)))) +
  #geom_text(data = t6_psets, aes(x = 20, y = 80, label = paste0('sigma = ', sigma))) +
  xlab("% Douglas-fir") +ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("A: Selection mosaic for C and s") + 
  ylim(0,100) + 
  facet_wrap(~id,nrow = 1) + 
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

plt3 <- nt4_int %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100,  group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  geom_text(data = nt4_psets, aes(x = 20, y = 90, label = paste0("LL = ",round(ll,1)))) +
  #geom_text(data = nt4_psets, aes(x = 20, y = 80, label = paste0('sigma = ', sigma))) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("D: No selection mosaic") + 
  ylim(0,100) + 
  facet_wrap(~id,nrow=1) + 
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

plt4 <- t4_int %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  geom_text(data = t4_psets, aes(x = 20, y = 90, label = paste0("LL = ", round(ll,1)))) +
  #geom_text(data = t4_psets, aes(x = 20, y = 80, label = paste0('sigma = ', sigma))) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("E: Selection mosaic for C") + 
  ylim(0,100) + 
  facet_wrap(~id,nrow=1) + 
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

plt5 <- t4ne_int %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  geom_text(data = t4ne_psets, aes(x = 20, y = 90, label = paste0("LL = ",round(ll,1)))) +
  #geom_text(data = t4ne_psets, aes(x = 20, y = 80, label = paste0('sigma = ', sigma))) +
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("F: Selection mosaic for C, no evolution") + 
  ylim(0,100) + 
  facet_wrap(~id,nrow = 1) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

plt6 <- t6nm_int %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  geom_text(data = t6nm_psets, aes(x = 20, y = 90, label = paste0('LL = ',round(ll,1)))) +
  #geom_text(data = t6nm_psets, aes(x = 20, y = 80, label = paste0('sigma = ', sigma))) +
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  theme(legend.position = 'none') + 
  ggtitle("B: Selection mosaic for s, C does not vary") + 
  ylim(0,100) + 
  facet_wrap(~id,nrow = 1) + 
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

pdf("figures/six_top5.pdf",height = 18, width = 16)
grid.arrange(plt2,plt6,plt1,plt3,plt4, plt5,nrow = 6)
#grid.arrange(plt2,plt3,nrow=3)
dev.off()

