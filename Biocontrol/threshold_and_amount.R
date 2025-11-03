library(tidyverse)
library(gridExtra)
library(gg3D)

dynamics <- read_csv("Biocontrol/data/biocontrol_dynamics.csv")

stats_vals <- dynamics %>% pivot_longer(cols = c(S,nu1,nu2), names_to = "Species", values_to = 'val')%>%
  group_by(pdoug,type,Species) %>% 
  mutate(nn = 1) %>% summarize(median = median(val),
                               q1 = quantile(val,0.25),
                               q2 = quantile(val,0.75),
                               l_wisker = boxplot.stats(val)$stats[1],
                               u_wisker = boxplot.stats(val)$stats[5],
                               n = sum(nn))

plt_snpv <- dynamics %>% mutate(type = factor(type,
                                           levels = c('No biopesticide applied','100% multi-capsid morphotype',
                                                      '75% multi-capsid, 25% single-capsid', '50% multi-capsid, 50% single-capsid',
                                                      '25% multi-capsid, 75% single-capsid',
                                                      "100% single-capsid morphotype"))) %>% 
  ggplot() + aes(x = pdoug/37*100, y = nu1, group = interaction(pdoug,type), fill = type) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(bar(nu)[SNPV])) + 
  scale_fill_manual("", values = c('No biopesticide applied' = 'grey65',
                                   '100% multi-capsid morphotype' = '#08519c',
                                   '75% multi-capsid, 25% single-capsid' = '#3182bd',
                                   '50% multi-capsid, 50% single-capsid' ='#fef0d9',
                                   '25% multi-capsid, 75% single-capsid' = '#fc8d59',
                                   "100% single-capsid morphotype" = '#e34a33'))+ 
  xlab("% Douglas-fir") + ggtitle("Single-capsid morphotype transmission risk") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')

plt_mnpv <- dfall %>% mutate(type = factor(type,
                                           levels = c('No biopesticide applied','100% multi-capsid morphotype',
                                                      '75% multi-capsid, 25% single-capsid', '50% multi-capsid, 50% single-capsid',
                                                      '25% multi-capsid, 75% single-capsid',
                                                      "100% single-capsid morphotype"))) %>% 
  ggplot() + aes(x = pdoug/37*100, y = nu2, group = interaction(pdoug,type), fill = type) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(bar(nu)[MNPV])) + 
  scale_fill_manual("", values = c('No biopesticide applied' = 'grey65',
                                   '100% multi-capsid morphotype' = '#08519c',
                                   '75% multi-capsid, 25% single-capsid' = '#3182bd',
                                   '50% multi-capsid, 50% single-capsid' ='#fef0d9',
                                   '25% multi-capsid, 75% single-capsid' = '#fc8d59',
                                   "100% single-capsid morphotype" = '#e34a33'))+ 
  xlab("% Douglas-fir") + ggtitle("Multi-capsid morphotype transmission risk") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')

plt_host <- dfall %>% mutate(type = factor(type,
                                           levels = c('No biopesticide applied','100% multi-capsid morphotype',
                                                      '75% multi-capsid, 25% single-capsid', '50% multi-capsid, 50% single-capsid',
                                                      '25% multi-capsid, 75% single-capsid',
                                                      "100% single-capsid morphotype"))) %>% 
  ggplot() + aes(x = pdoug/37*100, y = S, group = interaction(pdoug,type), fill = type) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(log[10]~"Host Population Size")) + 
  scale_fill_manual("", values = c('No biopesticide applied' = 'grey65',
                                   '100% multi-capsid morphotype' = '#08519c',
                                   '75% multi-capsid, 25% single-capsid' = '#3182bd',
                                   '50% multi-capsid, 50% single-capsid' ='#fef0d9',
                                   '25% multi-capsid, 75% single-capsid' = '#fc8d59',
                                   "100% single-capsid morphotype" = '#e34a33'))+ 
  xlab("% Douglas-fir") + ggtitle("Host Population") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

pdf("Biocontrol/figures/host_transmission_biocontrol.pdf",height = 8, width = 8)
grid.arrange(plt_snpv, plt_mnpv,plt_host,nrow = 3,heights = c(1,1,1.3))
dev.off()

all_avg <- read_csv('Biocontrol/data/host_density_bda.csv')

avg_vals <- all_avg %>% group_by(thresh,biocontrol,type,period, pdoug) %>% 
  summarize(mean_pop = mean(avg_pop),
            sd_pop = sd(avg_pop)) %>% filter(pdoug == 18)

pdf("Biocontrol/figures/average_density_thresh_amount.pdf",height = 8, width = 12)
all_avg %>% group_by(thresh,biocontrol,type,period, pdoug) %>% 
  summarize(mean_pop = mean(avg_pop),
            sd_pop = sd(avg_pop)) %>% filter(pdoug == 18) %>%
  ggplot() + aes(x = period, y = mean_pop, group = thresh, color = thresh) + geom_line() + 
  geom_point() + 
  geom_segment(aes(x = period, y = mean_pop - sd_pop, yend = mean_pop + sd_pop)) +
  theme_classic(base_size = 10) + facet_grid(biocontrol~type, labeller = label_bquote(rows = .(biocontrol)~'cad/m^2')) +
  scale_color_viridis_c(expression("Threshold (insects/"~m^2~")"), option = 'viridis') +
  scale_x_continuous("Biopesticide use (years)", breaks = c(1,2,3), labels = c("Before\n(50-149)", 'During\n(150-250)', "After\n(251-350)"), limits = c(0.8,3.2)) +
  theme(legend.position = 'top',
        legend.key.size = unit(0.75, 'cm'),
        panel.spacing.x = unit(0.5,'cm')) + 
  ylab("Average insect population density")
dev.off()

num_all <- read_csv('Biocontrol/data/num_applications.csv')

pdf("Biocontrol/figures/amount_applied.pdf",height = 4, width = 12)
num_all %>% ggplot() +
  aes(x = as.factor(thresh), y = as.factor(biocontrol), color = mean_applications, fill = mean_applications) +
  geom_tile() + theme_classic() + 
  geom_text(aes(x = as.factor(thresh), y = as.factor(biocontrol), label = round(mean_applications,1)), color = 'white', fill = NA) +
  facet_wrap(~type,nrow = 1) +
  scale_color_viridis_c("Avg. # of applications", option = 'turbo')+
  scale_fill_viridis_c("Avg. # of applications", option = 'turbo') +
  theme(legend.position = 'top') +
  xlab(expression("Threshold (insects/"~m^2~")")) +
  ylab(expression("Amount applied (cadavers/"~m^2~")"))
dev.off()

