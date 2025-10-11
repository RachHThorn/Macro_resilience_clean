# R Thornley
# 25/09/2025
# Check and plot model residuals for the interaction models

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

###############################################################################
# LOAD MODEL RESULTS; CALC RESIDUALS
###############################################################################

dat <- readRDS("results/models_robust_comms_rough.rds")
dat <- dat %>% filter(Comm_var == "variability")
dat <- dat %>% dplyr::select(mod) 
resids_long <- dat %>% mutate(resid = map(mod, residuals),
                              fitted  = map(mod, fitted)) %>%
  dplyr::select(resid, fitted) %>% unnest(cols = c(resid, fitted))


# Residuals vs Fitted, faceted by demo_trait
resids_long %>% 
  filter(experiment == "DIST") %>%
  filter(time_period == "T0-T1") %>%
  ggplot(aes(x = fitted, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.7) +
  facet_wrap(~ Demo_trait, scales = "free") +
  labs(title = "Residuals vs Fitted (robust models)",
       x = "Fitted values", y = "Residuals") +
  theme_bw()

################################################################################
# MARGINAL EFFECTS
################################################################################


# resids_long %>% 
filter(experiment == "NPK") %>%
  filter(time_period == "T0-T1") %>%
  ggplot(aes(x = fitted, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.7) +
  facet_wrap(~ Demo_trait, scales = "free") +
  labs(title = "Residuals vs Fitted (robust models)",
       x = "Fitted values", y = "Residuals") +
  theme_bw()

