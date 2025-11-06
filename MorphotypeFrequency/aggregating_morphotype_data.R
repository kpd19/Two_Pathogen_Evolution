library(tidyverse)
library(gridExtra)
library(mgcv)

site_inf <- read_csv("MorphotypeFrequency/data/site_infection_data_2019-2020.csv")
ll_data <- read_csv("MorphotypeFrequency/data/literature_dist_data.csv")
pd_data <- read_csv("MorphotypeFrequency/data/forest_composition_wide.csv")
tree_sp <- read_csv("MorphotypeFrequency/tree_polygons/tree_sp_df.csv")

# calculating the infection totals for all the natural sites that weren't included in biopesticide spray projects
site_inf_totals <- site_inf %>% filter(sprayed == FALSE) %>% group_by(state,latitude,longitude,site,year) %>%
  summarize(sum_MNPV = sum(MNPV,na.rm=TRUE), sum_SNPV = sum(SNPV,na.rm=TRUE),
            coinf = sum(coinfections,na.rm=TRUE),sum_total = sum(infected),
            sum_insect = sum(total),
            sum_infected = sum(infected)) %>% rename(MNPV = sum_MNPV,SNPV = sum_SNPV, Site = site) %>%
  mutate(total_iso = MNPV + SNPV,
         total_known = MNPV + SNPV - coinf,
         pSNPV = SNPV/total_iso,
         pMNPV = MNPV/total_iso,
         p_coinf = coinf/total_known) %>% 
  mutate(frac_inf = sum_infected/sum_insect,
         dom = max(pMNPV,pSNPV), ndom = min(pMNPV,pSNPV)) %>% mutate(rat = ndom/dom)

site_inf_totals %>% ungroup() %>% summarize(sum_coinf = sum(coinf),
                                            sum_tk = sum(total_known),
                                            sum_iso = sum(total_iso),
                                            sum_MNPV = sum(MNPV),
                                            sum_SNPV = sum(SNPV))

ll_small <- pd_data %>% rename(Douglas_fir = `Douglas-fir`, Site = site) %>% select(Site,Douglas_fir)
# analysing coinfections

site_inf_totals2 <- merge(site_inf_totals,ll_small)

crit <- 1.96

mod1 <- glm(p_coinf ~ rat, data = site_inf_totals2, family = 'quasibinomial')

x1 <- seq(0,0.84,0.01)
pred1 <- predict.glm(mod1, newdata = data.frame(rat = x1), type = 'response', se.fit = TRUE)
df1 <- data.frame(rat = x1, y = pred1$fit, se = pred1$se.fit)
df1 <- df1 %>% mutate(upper = y + 1.96*se,
                     lower = y - 1.96*se)
p1 <- summary(mod1)$coefficients[2,4]

site_inf_totals2 %>% ungroup() %>% filter(sum_infected >0) %>% 
  ggplot() +
  geom_ribbon(data = df1, aes(x = rat, ymin = lower, ymax = upper),
              fill = 'blue', alpha = 0.25, linetype = 'dashed', color = 'grey55')+
  geom_point(aes(x = sum_total,y = p_coinf), size = 2) + theme_classic(base_size = 15) + 
  xlab("Non-dominant:Dominant") + ylab("Proportion Coinfected") +
  geom_line(data = df1, aes(x = rat, y = y), color = 'blue', size = 1.25)     + 
  annotate(geom = 'text', x = -Inf, y = Inf, label =  paste0("P-val = ", round(p1,6)), hjust = -0.15, vjust = 1)

plt1 <- site_inf_totals2 %>% ungroup() %>% filter(sum_infected >0) %>% 
  ggplot() +
  geom_ribbon(data = df1, aes(x = rat, ymin = lower, ymax = upper),
              fill = 'blue', alpha = 0.25, linetype = 'dashed', color = 'grey55')+
  geom_point(aes(x = rat,y = p_coinf), size = 2) + theme_classic(base_size = 15) + 
  xlab("Non-dominant:Dominant") + ylab("Proportion Coinfected") +
  geom_line(data = df1, aes(x = rat, y = y), color = 'blue', size = 1.25)     + 
  annotate(geom = 'text', x = -Inf, y = Inf, label =  paste0("P-val = ", round(p1,6)), hjust = -0.15, vjust = 1)

mod2 <- glm(p_coinf ~ frac_inf, data = site_inf_totals, family = 'quasibinomial')

x2 <- seq(0.11,0.61,0.01)
pred2 <- predict.glm(mod2, newdata = data.frame(frac_inf = x2), type = 'response', se.fit = TRUE)
df2 <- data.frame(frac_inf = x2, y = pred2$fit, se = pred2$se.fit)
df2 <- df2 %>% mutate(upper = y + 1.96*se,
                      lower = y - 1.96*se)
p2 <- summary(mod2)$coefficients[2,4]

plt2 <- site_inf_totals2 %>% ungroup() %>% filter(sum_infected >0) %>% 
  ggplot() +
  geom_ribbon(data = df2, aes(x = frac_inf, ymin = lower, ymax = upper),
              fill = 'blue', alpha = 0.25, linetype = 'dashed', color = 'grey55')+
  geom_point(aes(x = frac_inf,y = p_coinf), size = 2) + theme_classic(base_size = 15) + 
  xlab("Fraction Infected") + ylab("Proportion Coinfected") +
  geom_line(data = df2, aes(x = frac_inf, y = y), color = 'blue', size = 1.25)    + 
  annotate(geom = 'text', x = -Inf, y = Inf, label =  paste0("P-val = ", round(p2,3)), hjust = -0.15, vjust = 1)

mod3 <- gam(p_coinf ~ s(pMNPV, k = 4), data = site_inf_totals)

x3 <- seq(0,1,0.01)
pred3 <- predict.gam(mod3, newdata = data.frame(pMNPV = x3), type = 'response', se.fit = TRUE)
df3 <- data.frame(pMNPV = x3, y = pred3$fit, se = pred3$se.fit)
df3 <- df3 %>% mutate(upper = y + 1.96*se,
                      lower = y - 1.96*se)
p3 <- summary(mod3)$s.table[,4]

plt3 <- site_inf_totals2 %>% ungroup() %>% filter(sum_infected >0) %>% 
  ggplot() +
  geom_ribbon(data = df3, aes(x = pMNPV, ymin = lower, ymax = upper),
              fill = 'blue', alpha = 0.25, linetype = 'dashed', color = 'grey55')+
  geom_point(aes(x = pMNPV,y = p_coinf), size = 2) + theme_classic(base_size = 15) + 
  xlab("Fraction Infected with\n Multi-capsid morphotype") + ylab("Proportion Coinfected") +
  geom_line(data = df3, aes(x = pMNPV, y = y), color = 'blue', size = 1.25)   + 
  annotate(geom = 'text', x = -Inf, y = Inf, label =  paste0("P-val = ", round(p3,3)), hjust = -0.15, vjust = 1)

mod4 <- gam(p_coinf ~ s(pSNPV, k = 4), data = site_inf_totals)

x4 <- seq(0,1,0.01)
pred4 <- predict.gam(mod4, newdata = data.frame(pSNPV = x4), type = 'response', se.fit = TRUE)
df4 <- data.frame(pSNPV = x4, y = pred4$fit, se = pred4$se.fit)
df4 <- df4 %>% mutate(upper = y + 1.96*se,
                      lower = y - 1.96*se)
p4 <- summary(mod4)$s.table[,4]

plt4 <- site_inf_totals2 %>% ungroup() %>% filter(sum_infected >0) %>% 
  ggplot() +
  geom_ribbon(data = df4, aes(x = pSNPV, ymin = lower, ymax = upper),
              fill = 'blue', alpha = 0.25, linetype = 'dashed', color = 'grey55')+
  geom_point(aes(x = pSNPV,y = p_coinf), size = 2) + theme_classic(base_size = 15) + 
  xlab("Fraction Infected with\n Single-capsid morphotype") + ylab("Proportion Coinfected") +
  geom_line(data = df4, aes(x = pSNPV, y = y), color = 'blue', size = 1.25)    + 
  annotate(geom = 'text', x = -Inf, y = Inf, label =  paste0("P-val = ", round(p4,3)), hjust = -0.15, vjust = 1)

mod5 <- glm(frac_inf ~ rat, data = site_inf_totals2, family = 'quasibinomial')

x5 <- seq(0,0.84,0.01)
pred5 <- predict.glm(mod5, newdata = data.frame(rat = x5), type = 'response', se.fit = TRUE)
df5 <- data.frame(rat = x5, y = pred5$fit, se = pred5$se.fit)
df5 <- df5 %>% mutate(upper = y + 1.96*se,
                      lower = y - 1.96*se)

p5 <- summary(mod5)$coefficients[2,4]

plt5 <- site_inf_totals2 %>% ungroup() %>% filter(sum_infected >0) %>% 
  ggplot() +
  geom_ribbon(data = df5, aes(x = rat, ymin = lower, ymax = upper),
              fill = 'blue', alpha = 0.25, linetype = 'dashed', color = 'grey55')+
  geom_point(aes(y = frac_inf,x = rat)) + theme_classic(base_size = 15) + 
  ylab("Fraction infected") + xlab("Non-dominant:Dominant") +
  geom_line(data = df5, aes(x = rat, y = y), color = 'blue', size = 1.25)   + 
  annotate(geom = 'text', x = -Inf, y = Inf, label =  paste0("P-val = ", round(p5,3)), hjust = -0.15, vjust = 1)

mod6 <- glm(p_coinf ~ Douglas_fir, data = site_inf_totals2, family = 'quasibinomial')

x6 <- seq(0,1,0.01)
pred6 <- predict.glm(mod6, newdata = data.frame(Douglas_fir = x6), type = 'response', se.fit = TRUE)
df6 <- data.frame(Douglas_fir = x6, y = pred6$fit, se = pred6$se.fit)
df6 <- df6 %>% mutate(upper = y + 1.96*se,
                      lower = y - 1.96*se)

p6 <- summary(mod6)$coefficients[2,4]

plt6 <- site_inf_totals2 %>% ungroup() %>% filter(sum_infected >0) %>% 
  ggplot() +
  geom_ribbon(data = df6, aes(x = Douglas_fir, ymin = lower, ymax = upper),
              fill = 'blue', alpha = 0.25, linetype = 'dashed', color = 'grey55')+
  geom_point(aes(x = Douglas_fir,y = p_coinf), size = 2) + theme_classic(base_size = 15) + 
  xlab("Proportion Douglas-fir") + ylab("Proportion Coinfected") +
  geom_line(data = df6, aes(x = Douglas_fir, y = y), color = 'blue', size = 1.25) + 
  annotate(geom = 'text', x = -Inf, y = Inf, label =  paste0("P-val = ", round(p6,3)), hjust = -0.15, vjust = 1)


pdf("morphotype_dist/figures/prop_coinf_ratio.pdf",height = 8, width = 12)
grid.arrange(plt2, plt1,plt5, plt3,plt4, plt6,nrow=2)
dev.off()

site_summary <- site_inf_totals %>% filter(total_iso >0) %>% arrange(Site) %>%
  select(Site,MNPV,SNPV,total_iso,state,latitude,longitude) %>% 
  rename(site = Site, total = total_iso) %>% mutate(elevation = NA, source = 'Dwyer Lab PCR Collections')
site_summary$site_no <- 111:128

ll_data <- rbind(ll_data,site_summary)

write_csv(ll_data, "MorphotypeFrequency/data/morphotye_distribution_data.csv")

head(ll_data)
head(pd_data)

pd_data2 <- pd_data %>% rename(Douglas_fir = `Douglas-fir`) %>% select(site,Douglas_fir)

ll_data2 <- merge(ll_data,pd_data2, all = TRUE)

site1 <- ll_data2 %>% ggplot() + aes(x = MNPV/(SNPV + MNPV)) + geom_histogram(bins = 11) + theme_classic(base_size = 15) + 
  xlab("% Multi-capsid morphotype") + ylab("# Sites")

site2 <- ll_data2 %>% ggplot() + aes(x = SNPV + MNPV) + geom_histogram(bins = 10) + theme_classic(base_size = 15) + 
  xlab("# Virus IDs") + ylab("# Sites") + 
  scale_x_log10()

site3 <- ll_data2 %>% ggplot() + aes(x = Douglas_fir*100) + geom_histogram(bins = 11) + theme_classic(base_size = 15) + 
  xlab("% Douglas-fir") + ylab("# Sites")

site4 <- ll_data2 %>% ggplot() + aes(x = state) + geom_bar() + theme_classic(base_size = 15) + 
  xlab("State") + ylab("# Sites")

pdf("MorphotypeFrequency/figures/site_info.pdf", height = 5, width = 15)
grid.arrange(site1,site2,site3,site4,nrow = 1)
dev.off()

pd_long <- pd_data %>% pivot_longer(colnames(pd_data)[11:38])

pd_long <- merge(pd_long,tree_sp)

pd_long <- pd_long %>% mutate(genus2 = ifelse(genus == 'Misc. gymnosperms', "Misc. gymnosperms", paste0("*",genus,"*")))

pdf("MorphotypeFrequency/figures/dist_tree_sp_by_sp.pdf",height = 10, width = 10)
pd_long %>% ggplot() + aes(x = site_no, y = value, fill = genus2) + geom_bar(stat = 'identity', position = 'stack') + 
  theme_classic(base_size = 9) + 
  facet_wrap(~name,nrow = 7) + 
  theme(legend.position = 'top',
        legend.text = element_markdown()) + 
  scale_fill_discrete("") + 
  ylab("% of tree species") + xlab("Site no.")
dev.off()

prop_dataset <- pd_long %>% group_by(name,genus) %>% summarize(mean_p = mean(value)) %>% arrange(desc(mean_p))

prop_dataset <- prop_dataset %>% 
  arrange(desc(mean_p))

prop_dataset$name <- factor(prop_dataset$name, levels = prop_dataset$name)

prop_dataset <- prop_dataset %>% mutate(genus2 = ifelse(genus == 'Misc. gymnosperms', "Misc. gymnosperms", paste0("*",genus,"*")))

pdf("MorphotypeFrequency/figures/tree_prop.pdf")
prop_dataset %>% ggplot() + aes(x = name, y = mean_p*100, fill = genus2) + geom_bar(stat = 'identity') + 
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("% of total dataset") + xlab("") + 
  theme(legend.position = 'top',
        legend.text = element_markdown()) +
  scale_fill_discrete("")
dev.off()
