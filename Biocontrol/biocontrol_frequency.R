library(tidyverse)
library(gridExtra)
library(gg3D)

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

SNPV <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/biocontrol/biocontrol_SNPV2.csv")
MIX <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/biocontrol/biocontrol_MIX.csv")

freq_MNPV <- read_csv("Biocontrol/data/biocontrol_probability_MNPV.csv")
freq_SNPV <- read_csv("Biocontrol/data/biocontrol_probability_SNPV.csv")
freq_MIX <- read_csv("Biocontrol/data/biocontrol_probability_MIX.csv")

freq_MNPV_long <- freq_MNPV %>% pivot_longer(cols = c('S','nu1','nu2'),
                                                   names_to = 'Species',values_to = "Column1") 

freq_SNPV_long <- freq_SNPV %>% pivot_longer(cols = c('S','nu1','nu2'),
                                                   names_to = 'Species',values_to = "Column1") 

freq_MIX_long <- freq_MIX %>% pivot_longer(cols = c('S','nu1','nu2'),
                                             names_to = 'Species',values_to = "Column1") 

thresh_pick <- 0.5
biocontrol_pick <- 150

legend_labels <- c(expression(bar(nu)["SNPV"]),expression(bar(nu)["MNPV"]))

pd_multi <- freq_MNPV_long %>% filter(thresh == thresh_pick, biocontrol == biocontrol_pick,
                               Species %in% c("nu2",'nu1'), Year >= 100, Year <= 300, rep == 1, prob !=0) %>% 
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                     '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                     '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd_num = recode(pdoug, '2' = 1, '7' = 2,'13' = 3,
                         '18' = 4, '24' = 5, 
                         '30' = 6, '35' = 7)) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir',
                                    '80% Douglas-fir', '95% Douglas-fir'))) %>% 
  ggplot() + aes(x = Year,y = Column1, group = Species, color = Species) + geom_path() + theme_classic(base_size = 10) +
  scale_y_log10() + 
  facet_grid(prob~pd, labeller = label_bquote(rows = p == .(prob))) + 
  geom_vline(aes(xintercept = 150),linetype = 'dashed', color = 'grey55') +
  scale_color_manual("", values = c(snpv_col, mnpv_col), labels = parse(text = legend_labels)) + 
  ylab("Transmission risk")

pd_single <- freq_SNPV_long %>% filter(thresh == thresh_pick, biocontrol == biocontrol_pick,
                                      Species %in% c("nu2",'nu1'), Year >= 100, Year <= 300, rep == 1, prob !=0) %>% 
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                     '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                     '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd_num = recode(pdoug, '2' = 1, '7' = 2,'13' = 3,
                         '18' = 4, '24' = 5, 
                         '30' = 6, '35' = 7)) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir',
                                    '80% Douglas-fir', '95% Douglas-fir'))) %>% 
  ggplot() + aes(x = Year,y = Column1, group = Species, color = Species) + geom_path() + theme_classic(base_size = 10) +
  scale_y_log10() + 
  facet_grid(prob~pd, labeller = label_bquote(rows = p == .(prob))) + 
  geom_vline(aes(xintercept = 150),linetype = 'dashed', color = 'grey55') +
  scale_color_manual("", values = c(snpv_col, mnpv_col), labels = parse(text = legend_labels)) + 
  ylab("Transmission risk")


pd_mix <- freq_MIX_long %>% filter(thresh == thresh_pick, biocontrol == biocontrol_pick,
                                       Species %in% c("nu2",'nu1'), Year >= 100, Year <= 300, rep == 1, prob !=0) %>% 
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                     '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                     '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd_num = recode(pdoug, '2' = 1, '7' = 2,'13' = 3,
                         '18' = 4, '24' = 5, 
                         '30' = 6, '35' = 7)) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir',
                                    '80% Douglas-fir', '95% Douglas-fir'))) %>% 
  ggplot() + aes(x = Year,y = Column1, group = Species, color = Species) + geom_path() + theme_classic(base_size = 10) +
  scale_y_log10() + 
  facet_grid(prob~pd, labeller = label_bquote(rows = p == .(prob))) + 
  geom_vline(aes(xintercept = 150),linetype = 'dashed', color = 'grey55') +
  scale_color_manual("", values = c(snpv_col, mnpv_col), labels = parse(text = legend_labels)) + 
  ylab("Transmission risk")


pdf("Biocontrol/figures/transmission_risk_prob_MNPV.pdf",height = 9, width = 12)
pd_multi
dev.off()

pdf("Biocontrol/figures/transmission_risk_prob_SNPV.pdf",height = 9, width = 12)
pd_single
dev.off()


pdf("Biocontrol/figures/transmission_risk_prob_MIX.pdf",height = 9, width = 12)
pd_mix
dev.off()

# freq_MIX, freq_SNPV, freq_MNPV
dat <- freq_MIX

plt_snpv <- dat %>% filter(Year >150, prob %in% c(0,0.2,0.4,0.6,0.8,1)) %>% 
  ggplot() + aes(x = pdoug/37*100, y = nu1, group = interaction(pdoug,prob), fill = as.factor(prob)) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(bar(nu)[SNPV])) + 
  xlab("% Douglas-fir") + ggtitle("Single-capsid morphotype transmission risk") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_brewer("Probability of \nusing biocontrol", palette = 'GnBu')

plt_mnpv <- dat %>% filter(Year >150, prob %in% c(0,0.2,0.4,0.6,0.8,1)) %>% 
  ggplot() + aes(x = pdoug/37*100, y = nu2, group = interaction(pdoug,prob), fill = as.factor(prob)) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(bar(nu)[SNPV])) + 
  xlab("% Douglas-fir") + ggtitle("Multi-capsid morphotype transmission risk") + 
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_fill_brewer("Probability of \nusing biocontrol", palette = 'GnBu')

plt_host <- dat %>% filter(Year >150, prob %in% c(0,0.2,0.4,0.6,0.8,1)) %>% 
  ggplot() + aes(x = pdoug/37*100, y = S, group = interaction(pdoug,prob), fill = as.factor(prob)) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(bar(nu)[SNPV])) + 
  xlab("% Douglas-fir") + ggtitle("Host population") + 
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_fill_brewer("Probability of \nusing biocontrol", palette = 'GnBu')

pdf("Biocontrol/figures/host_transmission_probability_MIX.pdf",height = 10, width = 12)
grid.arrange(plt_snpv, plt_mnpv,plt_host,nrow = 3)
dev.off()


# freq_MIX, freq_SNPV, freq_MNPV
dat <- freq_SNPV

freq <- dat %>% filter(prob != 0) %>% group_by(pdoug,prob,rep) %>%
  summarize(sum_pres = sum(TMB_pres), 
            sum_spark = sum(TMB_spark)) %>% 
  mutate(p = sum_pres/sum_spark)

pdf("Biocontrol/figures/biocontrol_freq_SNPV.pdf",height = 10, width = 12)
freq %>% group_by(pdoug,prob) %>%
  summarize(sum_pres = sum(sum_pres)/20, 
            sum_spark = sum(sum_spark)/20) %>% 
  mutate(p = sum_pres/sum_spark) %>% mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                                                        '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                                                        '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd_num = recode(pdoug, '2' = 1, '7' = 2,'13' = 3,
                         '18' = 4, '24' = 5, 
                         '30' = 6, '35' = 7)) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir',
                                    '80% Douglas-fir', '95% Douglas-fir'))) %>% 
  ggplot() + aes(x = pd, y = prob, color = sum_pres, fill = sum_pres) + geom_tile()  +
  geom_text(aes(x = pd, y = prob, label = sum_pres), color = 'orange', size = 10) + 
  theme_classic(base_size = 15) + 
  scale_color_viridis_c("Avg. # \nsprays in 100 years") + 
  scale_fill_viridis_c("Avg. # \nsprays in 100 years") + 
  ylab("Probability of Applying Biocontrol") + 
  xlab("") + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0))
dev.off()

# freq_MIX, freq_SNPV, freq_MNPV
SNPV_control <- freq_SNPV %>% filter(Year >=150) %>% group_by(pdoug,prob,rep) %>%
  mutate(csum_S = cumsum(S)) %>% mutate(type = "100% single-capsid morphotype")
MNPV_control <- freq_MNPV %>% filter(Year >=150) %>% group_by(pdoug,prob,rep) %>%
  mutate(csum_S = cumsum(S)) %>% mutate(type = '100% multi-capsid morphotype')
MIX_control <- freq_MIX %>% filter(Year >=150) %>% group_by(pdoug,prob,rep) %>%
  mutate(csum_S = cumsum(S)) %>% mutate(type = '50% multi-capsid, 50% single-capsid')

all_control <- rbind(SNPV_control,MNPV_control,MIX_control)

pdf("Biocontrol/figures/spray_effectiveness.pdf",height = 6, width = 8)
all_control %>% filter(pdoug == 18) %>%
  group_by(type,Year,pdoug,prob) %>% summarize(mean_csum = mean(csum_S),
                                               sd_csum = sd(csum_S)) %>% 
  filter(Year == 249) %>% 
  ggplot() + aes(x = prob, mean_csum, group = type, color = type) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_csum - sd_csum/sqrt(20),
                    ymax = mean_csum + sd_csum/sqrt(20)), size = 1, width = 0.05) +
  theme_classic() + 
  scale_color_manual("", values = c(mnpv_col,snpv_col,'#1b9e77')) +
  ylab("Cumulative Host Population over 100 years") + 
  xlab("Probability of Conducting Biocontrol") + 
  ggtitle("50% Douglas-fir Forest")+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'top') + 
  facet_wrap(~pdoug, scales = 'free')
dev.off()
