library(tidyverse)
library(gridExtra)

ll_data <- read_csv("data/morphotype_dist_data.csv")

id_noext <- ll_data %>% mutate(n = 1) %>%  group_by(n_trees) %>% mutate(csum = cumsum(n)) %>% filter(csum ==1) %>% pull(id)

field_data_round <- ll_data %>% mutate(round_pd = round(Douglas_fir*9)/9) %>% group_by(round_pd) %>% 
  summarize(sum_MNPV = sum(MNPV), sum_tot = sum(total),
            pMNPV = mean(MNPV/total),
            sd_pMNPV = sd(MNPV/total),
            len = length(MNPV),
            pdoug = mean(Douglas_fir),
            sd_pd = sd(Douglas_fir)) %>% 
  mutate(se_pMNPV = sd_pMNPV/sqrt(len),
         se_pd = sd_pd/sqrt(len))

ll_data2 <- ll_data %>% select(id,n_trees,MNPV,total) %>% rename(pdoug = n_trees)

t5_stoch <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/linesearch/t5s_noext_fix2_sig1.csv")

best_t5s <- t5_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS_DO,sS_GR,sM,phi,rho,round) %>% head(30)

write_csv(best_t5s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t5s_t30_best_fix2.csv")

t5_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-300)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_t5s <- t5_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_t5s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t5s_t30_tried_fix2.csv")
tried_t5s <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t5s_t30_tried_fix2.csv")

df_low <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t5s_fix2_lowstoch_tried1.csv")
df_no <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t5s_fix2_nostoch_tried1.csv")
df_high <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t5s_fix2_stoch_tried1.csv")
df_higher <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t5s_fix2_high_tried1.csv")
df_mid1 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t5s_fix2_mid1_tried1.csv")
df_mid2 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t5s_fix2_mid2_tried1.csv")

ll_low <- df_low %>% mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_no <- df_no %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_stoch <- df_high %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_higher <- df_higher %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid1 <- df_mid1 %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid2 <- df_mid2 %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

t5s_pset_info <- tried_t5s %>% rename(sM_DO = sM) %>% mutate(pset = 1:30) %>% select(sS_DO,sS_GR,sM_DO,phi,rho,pset)

t5s_all <- rbind(ll_no, ll_low, ll_stoch, ll_higher, ll_mid1,ll_mid2)

t5s_all <- merge(t5s_all,t5s_pset_info)

t5s_all <- t5s_all %>% arrange(desc(ll)) 

write_csv(t5s_all, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/t5s_stoch_all.csv")

t5s_all  %>% ggplot() +
  aes(x = sigma, y = ll, group = pset, color = pset) + geom_line() + theme_classic(base_size = 15) + 
  ggtitle("Line search with sigma = 0.5")+ 
  theme(legend.position ='none')

q1 <- 0.025
q2 <- 0.975

pmnpv0 <- df_no %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv005<- df_low %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv05 <- df_high %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv025 <- df_mid1 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv075 <- df_mid2 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv1 <- df_higher %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv <- rbind(pmnpv0,pmnpv005,pmnpv05,pmnpv025,pmnpv075,pmnpv1)
pmnpv <- pmnpv %>% rename(sM = sM_DO) %>% ungroup() %>% select(-sM_GR)

pmnpv <- merge(pmnpv,t5s_pset_info)

write_csv(pmnpv,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/pmnpv_t5s_fixnoext.csv")
#pmnpv <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix/summary/pnpv_t5_fixnoext.csv")


pdf("_plots/_eeid23/t5s_fix_noext_over_data.pdf",height = 10, width = 20)
pmnpv %>% filter(pset == 15, sigma == 0.5) %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, fill = as.factor(sigma), group = sigma)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = pdoug*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = pdoug*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = pdoug*100 - se_pd*100, xmax = pdoug*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  #geom_text(data = subset(stoch_df,sigma ==0),aes(x = 20, y = 90, label = round(ll,1))) +
  #geom_text(data = subset(stoch_df,sigma ==0.5),aes(x = 20, y = 70, label = round(ll,1))) + 
  xlab("% Douglas-fir") + ylab("% MNPV") + 
  facet_grid(sigma~pset)
dev.off()

int <- pmnpv %>% filter(pset == 15,pdoug == 32, sigma == 0)

probs <- df_no %>% filter(pdoug == 32, sS_DO == int$sS_DO)

pdf("_plots/_extremes/linesearch_fix_comp.pdf",height = 4, width = 10)
grid.arrange(plt1,plt2,nrow = 1)
dev.off()

nt4_stoch <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/linesearch/nt4s_noext_fix2_sig1.csv")

best_nt4s <- nt4_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS,sM,phi,rho,round) %>% head(30)

write_csv(best_nt4s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/nt4s_t30_best_fix2.csv")

nt4_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-300)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_nt4s <- nt4_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_nt4s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/nt4s_t30_tried_fix2.csv")
tried_nt4s <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/nt4s_t30_tried_fix2.csv")

df_low <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/nt4s_fix2_lowstoch_tried1.csv")
df_no <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/nt4s_fix2_nostoch_tried1.csv")
df_high <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/nt4s_fix2_stoch_tried1.csv")
df_higher <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/nt4s_fix2_high_tried1.csv")
df_mid1 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/nt4s_fix2_mid1_tried1.csv")
df_mid2 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/nt4s_fix2_mid2_tried1.csv")
df_higher2 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/nt4s_fix2_high2_tried1.csv")
#df_higher3 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/nt4s_fix2_higher2_tried1.csv")


ll_low <- df_low %>% mutate() %>% group_by(id,sS,sM,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS,sM,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_no <- df_no %>%  mutate() %>% group_by(id,sS,sM,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS,sM,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_stoch <- df_high %>%  mutate() %>% group_by(id,sS,sM,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS,sM,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_higher <- df_higher %>% mutate() %>% group_by(id,sS,sM,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS,sM,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid1 <- df_mid1 %>% mutate() %>% group_by(id,sS,sM,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS,sM,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid2 <- df_mid2 %>% mutate() %>% group_by(id,sS,sM,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS,sM,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_higher2 <- df_higher2 %>% group_by(id,sS,sM,phi,rho,sigma) %>%
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS,sM,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

# ll_higher3 <- df_higher3 %>% group_by(id,sS,sM,phi,rho,sigma) %>%
#   summarize(meanL = mean(LL,na.rm=TRUE)) %>%
#   group_by(sS,sM,phi,rho,sigma) %>%
#   summarize(ll = sum(log(meanL)))

tried_nt4s$pset <- 1:30
nt4s_pset_info <- tried_nt4s %>% select(sS,sM,phi,rho,pset)

nt4s_all <- rbind(ll_no, ll_low, ll_stoch, ll_higher, ll_mid1,ll_mid2)

#nt4s_pset_info <- nt4s_all %>% mutate(n = 1) %>% group_by(sS,sM,phi,rho) %>% summarize(sum_n = sum(n)) %>% ungroup() %>% mutate(pset = 1:30) %>% 
#  select(-sum_n)

nt4s_all <- merge(nt4s_all,nt4s_pset_info)

#write_csv(nt4s_all, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/nt4s_fix2_all.csv")
nt4s_all <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/nt4s_fix2_all.csv")
ll_higher2 <- merge(ll_higher2,nt4s_pset_info)

nt4s_all <- rbind(nt4s_all,ll_higher2)
write_csv(nt4s_all, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/nt4s_fix2_all.csv")

nt4s_all %>% arrange(desc(ll)) %>% head(10)

pdf("_plots/_fix2/nt4s_stoch.pdf",height = 8, width = 10)
nt4s_all %>% ggplot() +
  aes(x = sigma, y = ll, group = pset, color = pset) + geom_line() + theme_classic(base_size = 15) + 
  theme(legend.position ='none') + 
  ggtitle("no tree model")
dev.off()

nt4s_all %>% filter(ll >= -525) %>% ggplot() +
  aes(x = as.factor(phi), y = as.factor(sigma), fill =ll, color = ll) + geom_tile() + theme_classic(base_size = 15) + 
  ggtitle("no tree model") + 
  scale_color_viridis_c(option = 'turbo') + 
  scale_fill_viridis_c(option = 'turbo')


q1 <- 0.025
q2 <- 0.975

pmnpv0 <- df_no %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS,sM,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv005<- df_low %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS,sM,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv05 <- df_high %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS,sM,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv025 <- df_mid1 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS,sM,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv075 <- df_mid2 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS,sM,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv1 <- df_higher %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS,sM,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv125 <- df_higher2 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS,sM,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv <- rbind(pmnpv0,pmnpv005,pmnpv05,pmnpv025,pmnpv075,pmnpv1)#pmnpv125)

pmnpv <- merge(pmnpv,nt4s_pset_info)

write_csv(pmnpv,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/pmnpv_nt4s_fix2noext.csv")

pmnpv <- pmnpv %>% select(-pdoug)
pmnpv_nt <- merge(pmnpv,data.frame(pdoug = 1:36))

pdf("_plots/_fix2/nt4s_fix_noext_over_data.pdf",height = 10, width = 20)
pmnpv_nt %>%  filter(pset %in% c(16,2,10,25)) %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = pdoug*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = pdoug*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = pdoug*100 - se_pd*100, xmax = pdoug*100 + se_pd*100), size = 1) +
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  #geom_text(data = subset(stoch_df,sigma ==0),aes(x = 20, y = 90, label = round(ll,1))) +
  #geom_text(data = subset(stoch_df,sigma ==0.5),aes(x = 20, y = 70, label = round(ll,1))) + 
  xlab("% Douglas-fir") + ylab("% MNPV") + 
  facet_grid(pset~sigma)
dev.off()

sig1 <- nt4s_all %>% filter(sigma == 0) %>% arrange(desc(ll)) %>% head(5)
sig2 <- nt4s_all %>% filter(sigma > 0) %>% arrange(desc(ll)) %>% head(5)

sig1 <- rbind(sig1,sig2)
write_csv(sig1,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix/best/nt4s_sigexp1.csv")
sig_exp <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix/stoch/sig_exp1.csv")

sig_ll <- sig_exp %>% group_by(id,sS,sM,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS,sM,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

sig_ll <- sig_ll %>% arrange(desc(ll)) %>% rename(ll_8k = ll)
sig1 <- sig1 %>% rename(ll_2k = ll)

sig_ll <- merge(sig1,sig_ll)
sig_ll %>% select(-pset) %>% arrange(desc(ll_8k))


write_csv(sig_ll, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix/stoch/sigtest_stoch_all.csv")


rlzns <- seq(500,8000,500)

rlzn_test <- c()
for(i in 1:length(rlzns)){
  temp <- sig_exp %>% filter(rep <= rlzns[i]) %>%
    group_by(id,sS,sM,phi,rho,sigma) %>% 
    summarize(meanL = mean(LL,na.rm=TRUE)) %>%
    group_by(sS,sM,phi,rho,sigma) %>%
    summarize(ll = sum(log(meanL)))  
  temp$rep <- rlzns[i]
  
  rlzn_test <- rbind(rlzn_test,temp)
  print(rlzns[i])  
}

mini <- sig_ll %>% select(-ll_2k) %>% select(-ll_8k)

rlzn_test2 <- merge(rlzn_test,mini)

rlzn_test2 %>% ggplot() + aes(x = rep, y = ll, color = sigma, group = interaction(pset,sigma)) +
  geom_line() + theme_classic(base_size = 15) +
  scale_color_viridis_c(option = 'turbo')

q1 <- 0.025
q2 <- 0.975

pmnpv_sig <- sig_exp %>% filter(rep <= 2000) %>% filter(id %in% 1) %>%
  group_by(pdoug,sS,sM,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv_sig <- merge(pmnpv_sig,mini)

write_csv(pmnpv_sig,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix/summary/pnpv_sig_test_fixnoext.csv")


###################
###################
# T6 STOCH
###################
###################

t6_stoch <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/linesearch/t6s_noext_fix2_sig1.csv")

best_t6s <- t6_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,round) %>% head(30)

write_csv(best_t6s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t6s_noext_fix2_sig1.csv")

t6_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-300)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_t6s <- t6_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
#write_csv(tried_t6s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t6s_t30_tried_fix2.csv")
tried_t6s <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t6s_t30_tried_fix2.csv")

df_low <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t6s_fix2_lowstoch_tried1.csv")
df_no <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t6s_fix2_nostoch_tried1.csv")
df_high <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t6s_fix2_stoch_tried1.csv")
df_higher <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t6s_fix2_high_tried1.csv")
df_mid1 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t6s_fix2_mid1_tried1.csv")
df_mid2 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t6s_fix2_mid2_tried1.csv")

ll_low <- df_low %>% mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_no <- df_no %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_stoch <- df_high %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_higher <- df_higher %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid1 <- df_mid1 %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid2 <- df_mid2 %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

tried_t6s$pset <- 1:30
t6s_pset_info <- tried_t6s %>% select(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,pset)

t6s_all <- rbind(ll_no, ll_low, ll_stoch, ll_higher, ll_mid1,ll_mid2)

t6s_all <- merge(t6s_all,t6s_pset_info)

t6s_all <- t6s_all %>% arrange(desc(ll))
t6s_best <- t6s_all %>% arrange(desc(ll)) %>% head(1)

write_csv(t6s_all, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/t6s_stoch_all.csv")
write_csv(t6s_best, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/t6s_stoch_best1.csv")

pdf("_plots/_extremes/t6s_stoch.pdf",height = 8, width = 10)
plt2 <- t6s_all %>% ggplot() +
  aes(x = sigma, y = ll, group = pset, color = pset) + geom_line() + theme_classic(base_size = 15) + 
  ggtitle("Line search with sigma = 0.5")+ 
  theme(legend.position ='none')
plt2
dev.off()

q1 <- 0.025
q2 <- 0.975

pmnpv0 <- df_no %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv005<- df_low %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv05 <- df_high %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv025 <- df_mid1 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv075 <- df_mid2 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv1 <- df_higher %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv <- rbind(pmnpv0,pmnpv005,pmnpv05,pmnpv025,pmnpv075,pmnpv1)

pmnpv <- merge(pmnpv,t6s_pset_info)

write_csv(pmnpv,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/pmnpv_t6s_fix2noext.csv")
#pmnpv <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/pmnpv_t6s_fix2noext.csv")

pdf("_plots/_extremes/nt4s_stoch.pdf",height = 8, width = 10)

t6_psets <- t6s_all %>% select(sS_DO,sS_GR, sM_DO, phi, rho, ll,sigma,pset) %>% drop_na(ll) %>% arrange(desc(ll)) %>% head(10) %>% mutate(id = 1:10) 

t6_int <- c()

for(i in 1:length(t6_psets$pset)){
  temp <- pmnpv %>% filter(pset == t6_psets$pset[i], sigma == t6_psets$sigma[i])
  temp$id <- t6_psets$id[i]
  
  t6_int <- rbind(t6_int,temp)
}

pdf("_plots/_extremes/t6s_top10.pdf",height = 10, width = 20)
t6_int %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = pdoug*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = pdoug*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = pdoug*100 - se_pd*100, xmax = pdoug*100 + se_pd*100), size = 1) +
  geom_text(data = t6_psets, aes(x = 20, y = 90, label = round(ll,1))) +
  geom_text(data = t6_psets, aes(x = 20, y = 80, label = paste0('sigma = ', sigma))) +
  xlab("% Douglas-fir") + ylab("% MNPV") +
  theme(legend.position = 'none') + 
  ggtitle("Model with host tree effects") + 
  ylim(0,100) + 
  facet_wrap(~id,nrow = 2)
dev.off()

###################
###################
# T4 STOCH
###################
###################

t4_stoch <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/linesearch/t4s_noext_fix2_sig1.csv")

best_t4s <- t4_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS,sM,phi,rho,round) %>% head(30)

write_csv(best_t4s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t4s_t30_best_fix2.csv")

t4_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-400)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_t4s <- t4_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_t4s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t4s_t30_tried_fix2.csv")
tried_t4s <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/t4s_t30_tried_fix2.csv")


df_low <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t4s_fix2_lowstoch.csv")
df_no <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t4s_fix2_nostoch.csv")
df_high <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t4s_fix2_stoch.csv")
df_higher <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t4s_fix2_high.csv")
df_mid1 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t4s_fix2_mid1.csv")
df_mid2 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/t4s_fix2_mid2.csv")


ll_low <- df_low %>% mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_no <- df_no %>%  mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_stoch <- df_high %>%  mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_higher <- df_higher %>% mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid1 <- df_mid1 %>% mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid2 <- df_mid2 %>% mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

tried_t4s$pset <- 1:30
t4s_pset_info <- tried_t4s %>% rename(sS_DO = sS, sM_DO = sM) %>% mutate(sM_GR = sM_DO, sS_GR = sS_DO) %>%
  select(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,pset)

t4s_all <- rbind(ll_no, ll_low, ll_stoch, ll_higher, ll_mid1,ll_mid2)

t4s_all <- merge(t4s_all,t4s_pset_info)

write_csv(t4s_all, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/t4s_fix2_all.csv")

t4s_all %>% arrange(desc(ll)) %>% head(10)

pdf("_plots/_fix2/t4s_stoch.pdf",height = 8, width = 10)
t4s_all %>% ggplot() +
  aes(x = sigma, y = ll, group = pset, color = pset) + geom_line() + theme_classic(base_size = 15) + 
  theme(legend.position ='none') + 
  ggtitle("no tree model")
dev.off()

t4s_all %>% filter(ll >= -700) %>% ggplot() +
  aes(x = as.factor(phi), y = as.factor(sigma), fill =ll, color = ll) + geom_tile() + theme_classic(base_size = 15) + 
  ggtitle("no tree model") + 
  scale_color_viridis_c(option = 'turbo') + 
  scale_fill_viridis_c(option = 'turbo')


q1 <- 0.025
q2 <- 0.975

pmnpv0 <- df_no %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv005<- df_low %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv05 <- df_high %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv025 <- df_mid1 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv075 <- df_mid2 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv1 <- df_higher %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv <- rbind(pmnpv0,pmnpv005,pmnpv05,pmnpv025,pmnpv075,pmnpv1)

pmnpv <- merge(pmnpv,t4s_pset_info)

write_csv(pmnpv,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/pmnpv_t4s_fixnoext.csv")


###################
###################
# T6 STOCH NO HET
###################
###################

ht6_stoch <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/linesearch/ht6s_noext_fix2_sig1.csv")

best_ht6s <- ht6_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,round) %>% head(30)

write_csv(best_ht6s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/ht6s_t30_best_fix2.csv")

ht6_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-300)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_ht6s <- ht6_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_ht6s,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/ht6s_t30_tried_fix2.csv")
tried_ht6s <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/best/ht6s_t30_tried_fix2.csv")

df_low <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/ht6s_fix2_lowstoch.csv")
df_no <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/ht6s_fix2_nostoch.csv")
df_high <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/ht6s_fix2_stoch.csv")
df_higher <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/ht6s_fix2_high.csv")
df_mid1 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/ht6s_fix2_mid1.csv")
df_mid2 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/stoch/ht6s_fix2_mid2.csv")

ll_low <- df_low %>% mutate() %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_no <- df_no %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_stoch <- df_high %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_higher <- df_higher %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid1 <- df_mid1 %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

ll_mid2 <- df_mid2 %>% group_by(id,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>%
  summarize(ll = sum(log(meanL)))

tried_ht6s$pset <- 1:30
ht6s_pset_info <- tried_ht6s %>% select(sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,pset)

ht6s_all <- rbind(ll_no, ll_low, ll_stoch, ll_higher, ll_mid1,ll_mid2)

ht6s_all <- merge(ht6s_all,ht6s_pset_info)

ht6s_all <- ht6s_all %>% arrange(desc(ll))
ht6s_best <- ht6s_all %>% arrange(desc(ll)) %>% head(1)

write_csv(ht6s_all, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/ht6s_stoch_all.csv")
write_csv(ht6s_best, "/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/ht6s_stoch_best1.csv")

pdf("_plots/_extremes/ht6s_stoch.pdf",height = 8, width = 10)
plt2 <- ht6s_all %>% ggplot() +
  aes(x = sigma, y = ll, group = pset, color = pset) + geom_line() + theme_classic(base_size = 15) + 
  ggtitle("Line search with sigma = 0.5")+ 
  theme(legend.position ='none')
plt2
dev.off()

q1 <- 0.025
q2 <- 0.975

pmnpv0 <- df_no %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv005<- df_low %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv05 <- df_high %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv025 <- df_mid1 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv075 <- df_mid2 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv1 <- df_higher %>% filter(id %in% id_noext) %>%
  group_by(pdoug,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv <- rbind(pmnpv0,pmnpv005,pmnpv05,pmnpv025,pmnpv075,pmnpv1)

pmnpv <- merge(pmnpv,ht6s_pset_info)

write_csv(pmnpv,"/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/pmnpv_ht6s_fix2noext.csv")
#pmnpv <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/linesearch_fix2/summary/pmnpv_t6s_fix2noext.csv")

ht6_psets <- ht6s_all %>% select(sS_DO,sS_GR, sM_DO, phi, rho, ll,sigma,pset) %>% drop_na(ll) %>% arrange(desc(ll)) %>% head(10) %>% mutate(id = 1:10) 

ht6_int <- c()

for(i in 1:length(ht6_psets$pset)){
  temp <- pmnpv %>% filter(pset == ht6_psets$pset[i], sigma == ht6_psets$sigma[i])
  temp$id <- ht6_psets$id[i]
  
  ht6_int <- rbind(ht6_int,temp)
}

pdf("_plots/_fix2/ht6s_top10.pdf",height = 10, width = 20)
ht6_int %>% 
  ggplot() + geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
                         fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100)) +
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 - sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  # geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100 + sd*100, group = sigma), color = 'grey55', linetype = 'dashed') + 
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = pdoug*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = pdoug*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = pdoug*100 - se_pd*100, xmax = pdoug*100 + se_pd*100), size = 1) +
  geom_text(data = ht6_psets, aes(x = 20, y = 90, label = round(ll,1))) +
  geom_text(data = ht6_psets, aes(x = 20, y = 80, label = paste0('sigma = ', sigma))) +
  xlab("% Douglas-fir") + ylab("% MNPV") +
  theme(legend.position = 'none') + 
  ggtitle("Model with host tree effects") + 
  ylim(0,100) + 
  facet_wrap(~id,nrow = 2)
dev.off()
