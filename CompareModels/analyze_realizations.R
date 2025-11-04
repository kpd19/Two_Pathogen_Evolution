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

t5_stoch <- read_csv("linesearch/t5s_noext_fix2_sig1.csv")

best_t5s <- t5_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS_DO,sS_GR,sM,phi,rho,round) %>% head(30)

t5_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-300)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_t5s <- t5_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)

write_csv(tried_t5s,"best/t5s_t30.csv")

tried_t5s <- read_csv("best/t5s_t30.csv")

# Files are quite large, available on request
df_low <- read_csv("stoch/t5s_fix2_lowstoch_tried1.csv")
df_no <- read_csv("stoch/t5s_fix2_nostoch_tried1.csv")
df_high <- read_csv("stoch/t5s_fix2_stoch_tried1.csv")
df_higher <- read_csv("stoch/t5s_fix2_high_tried1.csv")
df_mid1 <- read_csv("stoch/t5s_fix2_mid1_tried1.csv")
df_mid2 <- read_csv("stoch/t5s_fix2_mid2_tried1.csv")

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

write_csv(t5s_all, "data/t5s_stoch_all.csv")

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

write_csv(pmnpv,"data/pmnpv_t5s_fixnoext.csv")

nt4_stoch <- read_csv("linesearch/nt4s_noext_fix2_sig1.csv")

best_nt4s <- nt4_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS,sM,phi,rho,round) %>% head(30)

nt4_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-300)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_nt4s <- nt4_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_nt4s,"best/nt4s_t30.csv")
tried_nt4s <- read_csv("best/nt4s_t30.csv")

# Files are quite large, available on request
df_low <- read_csv("stoch/nt4s_fix2_lowstoch_tried1.csv")
df_no <- read_csv("stoch/nt4s_fix2_nostoch_tried1.csv")
df_high <- read_csv("stoch/nt4s_fix2_stoch_tried1.csv")
df_higher <- read_csv("stoch/nt4s_fix2_high_tried1.csv")
df_mid1 <- read_csv("stoch/nt4s_fix2_mid1_tried1.csv")
df_mid2 <- read_csv("stoch/nt4s_fix2_mid2_tried1.csv")
df_higher2 <- read_csv("stoch/nt4s_fix2_high2_tried1.csv")

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

tried_nt4s$pset <- 1:30
nt4s_pset_info <- tried_nt4s %>% select(sS,sM,phi,rho,pset)

nt4s_all <- rbind(ll_no, ll_low, ll_stoch, ll_higher, ll_mid1,ll_mid2, ll_higher2)

nt4s_all <- merge(nt4s_all,nt4s_pset_info)

write_csv(nt4s_all, "data/nt4s_fix2_all.csv")

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

write_csv(pmnpv,"data/pmnpv_nt4s_fix2noext.csv")

###################
###################
# T6 STOCH
###################
###################

t6_stoch <- read_csv("linesearch/t6s_noext_fix2_sig1.csv")

best_t6s <- t6_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,round) %>% head(30)

t6_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-300)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_t6s <- t6_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_t6s,"best/t6s_t30.csv")
tried_t6s <- read_csv("best/t6s_t30.csv")

# Files are quite large, available on request
df_low <- read_csv("stoch/t6s_fix2_lowstoch_tried1.csv")
df_no <- read_csv("stoch/t6s_fix2_nostoch_tried1.csv")
df_high <- read_csv("stoch/t6s_fix2_stoch_tried1.csv")
df_higher <- read_csv("stoch/t6s_fix2_high_tried1.csv")
df_mid1 <- read_csv("stoch/t6s_fix2_mid1_tried1.csv")
df_mid2 <- read_csv("stoch/t6s_fix2_mid2_tried1.csv")

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

write_csv(t6s_all, "data/t6s_stoch_all.csv")
write_csv(t6s_best, "data/t6s_stoch_best1.csv")

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

write_csv(pmnpv,"data/pmnpv_t6s_fix2noext.csv")

###################
###################
# T4 STOCH
###################
###################

t4_stoch <- read_csv("linesearch/t4s_noext_fix2_sig1.csv")

best_t4s <- t4_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS,sM,phi,rho,round) %>% head(30)

t4_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-400)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_t4s <- t4_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_t4s,"best/t4s_t30.csv")
tried_t4s <- read_csv("best/t4s_t30.csv")

# Files are quite large, available on request
df_low <- read_csv("stoch/t4s_fix2_lowstoch.csv")
df_no <- read_csv("stoch/t4s_fix2_nostoch.csv")
df_high <- read_csv("stoch/t4s_fix2_stoch.csv")
df_higher <- read_csv("stoch/t4s_fix2_high.csv")
df_mid1 <- read_csv("stoch/t4s_fix2_mid1.csv")
df_mid2 <- read_csv("stoch/t4s_fix2_mid2.csv")


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

write_csv(t4s_all, "data/t4s_fix2_all.csv")

t4s_all %>% arrange(desc(ll)) %>% head(10)

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

write_csv(pmnpv,"data/pmnpv_t4s_fixnoext.csv")


###################
###################
# T6 STOCH NO HET
###################
###################

ht6_stoch <- read_csv("linesearch/ht6s_noext_fix2_sig1.csv")

best_ht6s <- ht6_stoch %>% filter(type =="best") %>%
  group_by(round) %>% filter(iter == max(iter)) %>% arrange(desc(ll)) %>% 
  select(ll,sS_DO,sS_GR,sM_DO,sM_GR,phi,rho,round) %>% head(30)

ht6_stoch %>% filter(type == "best") %>% 
  ggplot() + aes(x = iter, y = ll, color = round, group = round) +
  geom_line() + theme_classic() + 
  coord_cartesian(ylim = c(-700,-300)) + 
  scale_color_viridis_c(option = 'turbo') + 
  xlab("Iteration") + ylab("Best likelihood tried")

tried_ht6s <- ht6_stoch %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_ht6s,"best/ht6s_t30.csv")
tried_ht6s <- read_csv("best/ht6s_t30.csv")

# Files are quite large, available on request
df_low <- read_csv("stoch/ht6s_fix2_lowstoch.csv")
df_no <- read_csv("stoch/ht6s_fix2_nostoch.csv")
df_high <- read_csv("stoch/ht6s_fix2_stoch.csv")
df_higher <- read_csv("stoch/ht6s_fix2_high.csv")
df_mid1 <- read_csv("stoch/ht6s_fix2_mid1.csv")
df_mid2 <- read_csv("stoch/ht6s_fix2_mid2.csv")

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

write_csv(ht6s_all, "data/ht6s_stoch_all.csv")
write_csv(ht6s_best, "data/ht6s_stoch_best1.csv")

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

write_csv(pmnpv,"data/pmnpv_ht6s_fix2noext.csv")
#pmnpv <- read_csv("data/pmnpv_t6s_fix2noext.csv")


###################
###################
# T4 NO EVOLUTION
###################
###################


t4s_noevo <- read_csv("linesearch/t4sne_linesearch.csv")

tried_t4s <- t4s_noevo %>% filter(type == 'tried') %>% arrange(desc(ll)) %>% head(30)
write_csv(tried_t4s,"best/t4s_noevo_t30.csv")

# Files are quite large, available on request
df_low <- read_csv("data/t4s_noevo_lowstoch.csv")
df_no <- read_csv("data/t4s_noevo_nostoch.csv")
df_high <- read_csv("data/t4s_noevo_stoch.csv")
df_higher <- read_csv("data/t4s_noevo_high.csv")
df_mid1 <- read_csv("data/t4s_noevo_mid1.csv")
df_mid2 <- read_csv("data/t4s_noevo_mid2.csv")
df_higher2 <- read_csv("data/t4s_noevo_higher.csv")

ll_low <- df_low %>% group_by(id,phiS,phiM,lambda,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(phiS,phiM,lambda,rho,sigma) %>%
  summarize(ll = sum(log(meanL))) %>% arrange(desc(ll))

ll_no <- df_no %>% group_by(id,phiS,phiM,lambda,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(phiS,phiM,lambda,rho,sigma) %>%
  summarize(ll = sum(log(meanL))) %>% arrange(desc(ll))

ll_stoch <- df_high %>% group_by(id,phiS,phiM,lambda,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(phiS,phiM,lambda,rho,sigma) %>%
  summarize(ll = sum(log(meanL))) %>% arrange(desc(ll))

ll_higher <- df_higher %>% group_by(id,phiS,phiM,lambda,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(phiS,phiM,lambda,rho,sigma) %>%
  summarize(ll = sum(log(meanL))) %>% arrange(desc(ll))

ll_mid1 <- df_mid1 %>% group_by(id,phiS,phiM,lambda,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(phiS,phiM,lambda,rho,sigma) %>%
  summarize(ll = sum(log(meanL))) %>% arrange(desc(ll))

ll_mid2 <- df_mid2 %>% group_by(id,phiS,phiM,lambda,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(phiS,phiM,lambda,rho,sigma) %>%
  summarize(ll = sum(log(meanL))) %>% arrange(desc(ll))

ll_higher2 <- df_higher2 %>% group_by(id,phiS,phiM,lambda,rho,sigma) %>% 
  summarize(meanL = mean(LL,na.rm=TRUE)) %>%
  group_by(phiS,phiM,lambda,rho,sigma) %>%
  summarize(ll = sum(log(meanL))) %>% arrange(desc(ll))

t4s_all <- rbind(ll_no, ll_low, ll_stoch, ll_higher, ll_mid1,ll_mid2,ll_higher2)

t4s_all <- merge(t4s_all,tried_t4s_info)

t4s_all <- t4s_all %>% arrange(desc(ll))

write_csv(t4s_all, "data/t4s_noevo_stoch_all.csv")

q1 <- 0.025
q2 <- 0.975

pmnpv0 <- df_no %>% filter(id %in% id_noext) %>%
  group_by(pdoug,phiS,phiM,rho,lambda,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv005<- df_low %>% filter(id %in% id_noext) %>%
  group_by(pdoug,phiS,phiM,rho,lambda,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv05 <- df_high %>% filter(id %in% id_noext) %>%
  group_by(pdoug,phiS,phiM,rho,lambda,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv025 <- df_mid1 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,phiS,phiM,rho,lambda,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv075 <- df_mid2 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,phiS,phiM,rho,lambda,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv1 <- df_higher %>% filter(id %in% id_noext) %>%
  group_by(pdoug,phiS,phiM,rho,lambda,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv125 <- df_higher2 %>% filter(id %in% id_noext) %>%
  group_by(pdoug,phiS,phiM,rho,lambda,sigma) %>% 
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE ),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pmnpv <- rbind(pmnpv0,pmnpv005,pmnpv05,pmnpv025,pmnpv075,pmnpv1,pmnpv125)


pmnpv <- merge(pmnpv,tried_t4s_info)

write_csv(pmnpv,"data/pmnpv_t4s_noevo.csv")