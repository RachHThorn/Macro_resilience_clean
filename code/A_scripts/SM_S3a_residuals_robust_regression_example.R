# R Thornley
# 25/09/2025
# Check and plot model residuals for the robust linear models - 
# one example for the appendices

###############################################################################
# LOAD MODEL RESULTS; CALC RESIDUALS
###############################################################################

dat <- readRDS("results/robust_model_results.rds")
names(dat)
dat$experiment
dat$time_period
dat <- dat %>% filter(experiment == "DIST") %>% filter(model == "Ordbeta") %>% filter(time_period == "T0-T1")
dat <- dat %>% dplyr::select(mod) 
dat <- dat %>% mutate(resid = purrr::map(mod, residuals),
                              fitted  = purrr::map(mod, fitted)) %>%
  dplyr::select(resid, fitted) %>% unnest(cols = c(resid, fitted))

# rename the life history traits for the plot
dat <- 
  dat %>% mutate(
    New_demo_trait = case_when(
      Demo_trait == "T_generation"    ~ "Generation time",
      Demo_trait == "R0"              ~ "Net reproductive rate",
      Demo_trait == "percapita_repro" ~ "Per capita reproduction",
      Demo_trait == "Lmean"           ~ "Mean life expectancy",
      Demo_trait == "age_repro"       ~ "Age at first reproduction",
      Demo_trait == "Lambda"          ~ "Population growth rate",
      Demo_trait == "Reactivity"      ~ "Reactivity",
      Demo_trait == "FirstStepAtt"    ~ "First step attenuation",
      Demo_trait == "Lmax"            ~ "Maximum longevity", # changed this name
      TRUE                            ~ Demo_trait
    ))

# Residuals vs Fitted, faceted by demo_trait
plot <- dat %>% 
  ggplot(aes(x = fitted, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.7) +
  facet_wrap(~ New_demo_trait, scales = "free") +
  labs(title = "Residuals vs Fitted (robust models)",
       x = "Fitted values", y = "Residuals") +
  theme_bw()+
  theme(plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )
plot
ggsave("figures/SM_S3a_residuals_robust_OLS_example.jpeg", plot, width = 9, height = 7)


