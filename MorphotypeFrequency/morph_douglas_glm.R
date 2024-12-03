library(tidyverse)
library(gridExtra)

ll_data <- read_csv("data/data_for_ll.csv")

ll_data <- ll_data %>% mutate(MNPV_perc = MNPV/total)

glm(MNPV_perc ~ Douglas_fir, data = subset(ll_data,total>=0), family = 'binomial',weights = total)

glm_mod <- glm(MNPV_perc ~ Douglas_fir, data = subset(ll_data,total>=0), family = 'binomial',weights = total)
glm_mod <- glm(MNPV_perc ~ Douglas_fir, data = subset(ll_data,total>=0), family = 'binomial',weights = total)
glm0 <- glm(MNPV_perc ~ 1, data = subset(ll_data,total>=0), family = 'binomial',weights = total)
summary(glm_mod)
summary(glm0)

logLik(glm_mod)
logLik(glm0)

AIC(glm0) - AIC(glm_mod)
stats::anova(glm0,glm_mod, test = "Chisq")

xnums <- seq(0,1,0.01)
ypred_S <- predict.glm(glm_mod, list(Douglas_fir = xnums), se.fit = TRUE)$fit
ypred_S_se <- predict.glm(glm_mod, list(Douglas_fir = xnums), se.fit = TRUE)$se.fit

ilink <- family(glm_mod)$linkinv

mod_df <- data.frame(df = xnums, mperc = ypred_S, se_pred = ypred_S_se)

mod_df <- mod_df %>% mutate(fit_resp = ilink(mperc),
                            right_upr = ilink(mperc) + (1.96*se_pred),
                            right_lwr = ilink(mperc) - (1.96*se_pred))

mod_df2 <- mod_df %>% mutate(df = df*100, mperc = mperc*100, fit_resp = fit_resp*100, se_pred = se_pred*100, right_upr = right_upr*100, right_lwr = right_lwr*100)

field_data_round <- ll_data %>% mutate(round_pd = round(Douglas_fir*9)/9) %>% group_by(round_pd) %>% 
  summarize(sum_MNPV = sum(MNPV), sum_tot = sum(total),
            pMNPV = mean(MNPV/total),
            sd_pMNPV = sd(MNPV/total),
            len = length(MNPV),
            pdoug = mean(Douglas_fir),
            sd_pd = sd(Douglas_fir)) %>% 
  mutate(se_pMNPV = sd_pMNPV/sqrt(len),
         se_pd = sd_pd/sqrt(len))

pdf("figures/glm_mod.pdf",height = 5, width = 9)
mod_df2 %>% ggplot() + geom_line(aes(x = df, y = fit_resp), color = 'black', size = 1.2) +
  geom_ribbon(aes(x = df, ymax = right_upr, ymin = right_lwr),
              fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') + theme_classic(base_size = 15) +
  geom_point(data = field_data_round, aes(x = pdoug*100, y = pMNPV*100, size = sum_tot)) +
  geom_errorbar(data = field_data_round, aes(x = pdoug*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = pdoug*100 - se_pd*100, xmax = pdoug*100 + se_pd*100), size = 1) +
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") + 
  theme(legend.position = 'none') + 
  scale_y_continuous(breaks = c(0,20,40,60,80))
dev.off()
