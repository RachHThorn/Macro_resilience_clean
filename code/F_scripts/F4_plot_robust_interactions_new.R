# R Thornley
# 09/10/2025
# Project: P1_COMPADRE_DRAGNET
# Script: F4_plot_robust_interactions_new

rm(list = ls())

################################################################################
# Instructions
################################################################################

# 1) LOAD MODELS AND FILTER FOR SIG INTERACTION EFFECTS
# 2) ITERATE OVER ALL MODELS IN THE DATASET
# 3) CREATE NEW PLOTTING VARS
# 4) REFINED FACETED PLOT

################################################################################
# 1) LOAD MODELS AND FILTER FOR SIG INTERACTION EFFECTS
################################################################################

dat <- readRDS("results/robust_model_comm_results_T0_scaled_not_logged.rds") %>%
  mutate(mod_name = paste0(hypoth, "_", Demo_trait, "_", Comm_var, "_", experiment, "_", time_period_cover_change))
names(dat)

# count number of models calculated
dat %>% 
  mutate(mod_name = paste0(hypoth, "_", Demo_trait, "_", Comm_var, "_", experiment, "_", time_period_cover_change)) %>%
  pull(mod_name) %>% length()
# 216 models

# get a list of the models we want
mods_wanted <- dat %>% 
  dplyr::select(time_period_cover_change, experiment, Demo_trait, Comm_var, results, hypoth) %>%
  unnest(results) %>%
  filter(term != "(Intercept)") %>%
  filter(p.value < 0.05 & term == "Demo_value:Comm_value") %>%
  mutate(mod_name = paste0(hypoth, "_", Demo_trait, "_", Comm_var, "_", experiment, "_", time_period_cover_change)) %>%
  ungroup() %>%
  pull(mod_name)

sig_mods <- dat %>% filter(mod_name %in% mods_wanted) %>% arrange(hypoth)
sig_mods
sig_mods %>% group_by(Comm_var) %>% tally() %>% arrange(-n)
list_sig_mods <- sig_mods %>% pull(mod_name)
list_sig_mods

# Filter for the Comm_vars we want to report on (we don't want to reporton some of them)
comm_wanted <- c("variability", "meanSigmas", "site_time_invsimp", "mean_cover")
comm_wanted <- sig_mods %>% filter(Comm_var %in% comm_wanted)
comm_wanted <- comm_wanted$mod_name
comm_wanted
# there are 39 significant models 

# create a significance flag
dat <- dat %>% mutate(sig_flag = if_else(mod_name %in% list_sig_mods, TRUE, FALSE))


# filter the data for the hypotheses and the comm_var you want
# we want to keep the non-sig models for all the times 
# functional diversity
sig_mods %>% filter(Comm_var == "variability")
FD <- dat %>% 
  filter(
    Comm_var == "variability" & hypoth == "H1" & Demo_trait == "R0" |
      Comm_var == "variability" & hypoth == "H3" & Demo_trait == "Reactivity" |
      Comm_var == "variability" & hypoth == "H5" & Demo_trait == "Lmax") %>%
  arrange(hypoth)

# specialisation
sig_mods %>% filter(Comm_var == "meanSigmas")
FS <- dat %>% 
  filter(
    Comm_var == "meanSigmas" & hypoth == "H1" & Demo_trait == "R0" |
      Comm_var == "meanSigmas" & hypoth == "H5" & Demo_trait == "Lmax") %>% 
  arrange(hypoth)

# then species diversity (invsimp)
sig_mods %>% filter(Comm_var == "site_time_invsimp")
SD <- dat %>% 
  filter(Comm_var == "site_time_invsimp" & hypoth == "H5" & Demo_trait == "Lmean")%>% 
  arrange(hypoth)


################################################################################
# 2) ITERATE OVER ALL MODELS IN THE DATASET: see what is looks like
# use the first approach with only the preds in range of data used in the model
################################################################################

# Function to perform the predictions 
get_preds <- function(mod, z_vals, n = 200) {
  
  # Extract model frame and ranges
  mf <- model.frame(mod)
  demo_range <- range(mf$Demo_value, na.rm = TRUE)
  
  # Build prediction grid
  grid <- expand.grid(
    Demo_value = seq(demo_range[1], demo_range[2], length.out = n),
    Comm_value = z_vals
  )
  
  # Get predictions
  pr <- predict(mod, newdata = grid, se.fit = TRUE, interval = "confidence")
  fit <- as.data.frame(pr$fit)
  se  <- pr$se.fit
  
  # Combine into a data frame
  pred_df <- cbind(grid, fit, se = se)
  pred_df$Comm_value <- factor(pred_df$Comm_value, levels = sort(unique(z_vals)))
  
  return(pred_df)
}

# we need to define the z vals - the level of the community variable at which the 
# slopes are calculated

# z vals for FD
# z_vals = c(-2, 0, 2)
# map this function over all the models
FD_pred <- FD %>% mutate(preds = purrr::map(mod, get_preds, z_vals = c(-2, 0, 2))) %>% 
  dplyr::select(time_period_cover_change, experiment, Demo_trait, hypoth, Comm_var, mod_name, preds) %>% 
  unnest(preds)

# z vals for FS
# z_vals = c(-2, 0, 2)
FS_pred <- FS %>% mutate(preds = purrr::map(mod, get_preds, z_vals = c(-2, 0, 2))) %>% 
  dplyr::select(time_period_cover_change, experiment, Demo_trait, hypoth, Comm_var, mod_name, preds) %>% 
  unnest(preds)

# z vals for SD
# z_vals = c(-2, 0, 2)
SD_pred <- SD %>% mutate(preds = purrr::map(mod, get_preds, z_vals = c(-1.5, 0, 2))) %>% 
  dplyr::select(time_period_cover_change, experiment, Demo_trait, hypoth, Comm_var, mod_name, preds, sig_flag) %>% 
  unnest(preds)

################################################################################
# 3) SELECT MODELS TO REPORT ON AND CREATE NEW PLOTTING VARS
################################################################################
# select manually which of the models we want in Fig. 4

list_sig_mods

# create some categorical values 
# functional diversity
unique(FD_pred$mod_name)
FD_plot_dat <- FD_pred %>%
  mutate(Comm_value_new = factor(Comm_value,
                                 levels = sort(unique(Comm_value)),
                                 labels = c("Low", "Medium", "High"))) %>%
  mutate(Plot_var = case_when(
    mod_name == "H1_R0_variability_DIST_T0-T1" ~ "H1\n Functional diversity \nReproductive output",
    mod_name == "H1_R0_variability_DIST_T0-T2" ~ "H1\n Functional diversity \nReproductive output",
    mod_name == "H1_R0_variability_DIST_T0-T3" ~ "H1\n Functional diversity \nReproductive output",
    mod_name == "H3_Reactivity_variability_DIST_T0-T1" ~ "H3\n Functional diversity \nReactivity",
    mod_name == "H3_Reactivity_variability_DIST_T0-T2" ~ "H3\n Functional diversity \nReactivity",
    mod_name == "H3_Reactivity_variability_DIST_T0-T3" ~ "H3\n Functional diversity \nReactivity",
    mod_name == "H5_Lmax_variability_NPK_T0-T1" ~ "H5\n Functional diversity \nMaximum life expectancy",
    mod_name == "H5_Lmax_variability_NPK_T0-T2" ~ "H5\n Functional diversity \nMaximum life expectancy",
    mod_name == "H5_Lmax_variability_NPK_T0-T3" ~ "H5\n Functional diversity \nMaximum life expectancy"))

# Functional specialisation
unique(FS_pred$mod_name)
FS_plot_dat <- FS_pred %>%
  mutate(Comm_value_new = factor(Comm_value,
                                 levels = sort(unique(Comm_value)),
                                 labels = c("Low", "Medium", "High"))) %>%
  mutate(Plot_var = case_when(
    mod_name == "H1_R0_meanSigmas_DIST_T0-T1" ~ "H1\n Functional specialisation \nReproductive output",
    mod_name == "H1_R0_meanSigmas_DIST_T0-T2" ~ "H1\n Functional specialisation \nReproductive output",
    mod_name == "H1_R0_meanSigmas_DIST_T0-T3" ~ "H1\n Functional specialisation \nReproductive output",
    mod_name == "H5_Lmax_meanSigmas_NPK_T0-T1" ~ "H5\n Functional specialisation \nMaximum life expectancy",
    mod_name == "H5_Lmax_meanSigmas_NPK_T0-T2" ~ "H5\n Functional specialisation \nMaximum life expectancy",
    mod_name == "H5_Lmax_meanSigmas_NPK_T0-T3" ~ "H5\n Functional specialisation \nMaximum life expectancy"))

# Species Diversity
unique(SD_pred$mod_name)
SD_plot_dat <- SD_pred %>%
  mutate(Comm_value_new = factor(Comm_value,
                                 levels = sort(unique(Comm_value)),
                                 labels = c("Low", "Medium", "High"))) %>%
  mutate(Plot_var = case_when(
    mod_name == "H5_Lmean_site_time_invsimp_NPK_T0-T1" ~ "H5\n Species diversity \nMean life expectancy",
    mod_name == "H5_Lmean_site_time_invsimp_NPK_T0-T2" ~ "H5\n Species diversity \nMean life expectancy",
    mod_name == "H5_Lmean_site_time_invsimp_NPK_T0-T3" ~ "H5\n Species diversity \nMean life expectancy"))


all_plot_dat <- rbind(FD_plot_dat, FS_plot_dat, SD_plot_dat)
unique(all_plot_dat$mod_name)

################################################################################
# 4) EXTRACT P VALUES AND EFFECTS FROM MODELS
################################################################################

FD_results <- 
  FD %>% 
  dplyr::select(time_period_cover_change, experiment, Demo_trait, hypoth, Comm_var, results) %>%
  unnest(results) %>%
  filter(term == "Demo_value:Comm_value")

FS_results <- 
  FS %>% 
  dplyr::select(time_period_cover_change, experiment, Demo_trait, hypoth, Comm_var, results) %>%
  unnest(results) %>%
  filter(term == "Demo_value:Comm_value")

SD_results <- 
  SD %>% 
  dplyr::select(time_period_cover_change, experiment, Demo_trait, hypoth, Comm_var, results) %>%
  unnest(results) %>%
  filter(term == "Demo_value:Comm_value")

all_results <- rbind(FD_results, FS_results, SD_results)
# tidy some names in all_results for joining 
all_results <-
  all_results %>% 
  mutate(mod_name = paste0(hypoth, "_", Demo_trait, "_", Comm_var, "_", experiment, "_", time_period_cover_change)) %>%
  dplyr::select(mod_name, estimate, std.error, p.value)
all_results

all_results <- all_results %>%
  mutate(sig_label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

################################################################################
# 5) CREATE PLOT DAT
################################################################################

# join to prediction values for plotting 
all_plot_dat_new <- all_plot_dat %>% left_join(all_results)

# check that we have all the models in the data frame after joining
setdiff(unique(all_plot_dat$mod_name), unique(all_results$mod_name))
intersect(names(all), names(all_results))

unique(all_plot_dat_new$mod_name)

################################################################################
# 6) PLOT BY HYPOTHESIS FOR BETTER VISUALISATION
################################################################################

# Load colour palette
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73")

################################################################################
# 7) PLOT BEST MODELS AS FULL PAGE FOR MANUSCRIPT
################################################################################

all_plot <- 
  ggplot(all_plot_dat_new, aes(x = Demo_value, y = fit,
                               colour = Comm_value_new,
                               fill = Comm_value_new,
                               linetype = Comm_value_new)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  theme_bw() +
  facet_grid(Plot_var ~ time_period_cover_change, scales = "free") +
  labs(
    x = "Demographic metric (scaled)",
    y = "Species abundance change at the site level",
    colour = "Community complexity metric (scaled)",
    fill   = "Community complexity metric (scaled)",
    linetype = "Community complexity metric (scaled)") +
  scale_fill_manual(values = okabe_ito, name = "Community complexity metric (scaled)") +
  scale_colour_manual(values = okabe_ito, name = "Community complexity metric (scaled)") +
  scale_linetype_manual(
    values = c("Low" = "dotted", "Medium" = "dashed", "High" = "solid"),
    name = "Community complexity metric (scaled)") +
  guides(fill = "none",
         colour = guide_legend(
           override.aes = list(fill = alpha("grey50", 0.2), size = 1.2)),
         linetype = guide_legend(
           override.aes = list(fill = alpha("grey50", 0.2), size = 1.2))) +
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(colour = "black", size = 8),
    legend.position = "top",
    legend.justification = "center",
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5, 15, 5, 5)
  ) +
  coord_cartesian(clip = "off")

all_plot

###############################################################################
# this code below is now redundant for this version of the manuscript

# # p value label
# p_lab <-  all_plot_dat_new %>%
#   group_by(time_period_cover_change, Plot_var) %>%
#   summarise(p_value = dplyr::first(na.omit(p.value)), .groups = "drop") %>%
#   mutate(
#     p_value = ifelse(p_value < 0.0001, "<0.0001",
#                      formatC(p_value, format = "f", digits = 4)))
# p_lab
# 
# e_lab <- all_plot_dat_new %>%
#   group_by(time_period_cover_change, Plot_var) %>%
#   summarise('Interaction estimate' = dplyr::first(na.omit(estimate)), .groups = "drop")
# e_lab
# 
# # append a single label for interaction estimate per facet
# all_plot_labs <- all_plot +
#   geom_text(
#     data = e_lab,
#     aes(label = paste0("Estimate = ", signif(`Interaction estimate`, 2))),
#     x = -Inf, y = Inf, hjust = -0.05, vjust = 1.7,
#     inherit.aes = FALSE, size = 3, show.legend = FALSE
#   ) +
#   geom_text(
#     data = p_lab,
#     aes(label = paste0("p-value = ", p_value)),
#     x = -Inf, y = Inf, hjust = -0.045, vjust = 4,
#     inherit.aes = FALSE, size = 3, show.legend = FALSE
#   )+
#   xlim(-2.5, 10)

################################################################################

ggsave("figures/Fig_4_predicted_interactions_new_09_11_2025.pdf", plot = all_plot, width = 17, height = 21, units = "cm", dpi = 600)

###############################################################################
# alternative - so all x axes are free
###############################################################################

install.packages("ggh4x")
library(ggh4x)

new_plot <-
  ggplot(all_plot_dat_new, aes(x = Demo_value, y = fit,
                               colour = Comm_value_new,
                               fill = Comm_value_new,
                               linetype = Comm_value_new)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  theme_bw() +
  ggh4x::facet_grid2(
    rows   = vars(Plot_var),
    cols   = vars(time_period_cover_change),
    scales = "free",        # allow different x ranges
    independent = "x"         # make x scales independent per *panel*
  ) +
  labs(
    x = "Demographic metric (scaled)",
    y = "Species abundance change at the site level",
    colour = "Community complexity metric (scaled)",
    fill   = "Community complexity metric (scaled)",
    linetype = "Community complexity metric (scaled)") +
  scale_fill_manual(values = okabe_ito, name = "Community complexity metric (scaled)") +
  scale_colour_manual(values = okabe_ito, name = "Community complexity metric (scaled)") +
  scale_linetype_manual(values = c("Low" = "dotted", "Medium" = "dashed", "High" = "solid"),
                        name = "Community complexity metric (scaled)") +
  guides(fill = "none",
         colour = guide_legend(override.aes = list(fill = alpha("grey50", 0.2), size = 1.2)),
         linetype = guide_legend(override.aes = list(fill = alpha("grey50", 0.2), size = 1.2))) +
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(colour = "black", size = 8),
    legend.position = "top",
    legend.justification = "center",
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5, 15, 5, 5)
  ) +
  coord_cartesian(clip = "off")
new_plot

# lets use this plot with the completely independent axes
ggsave("figures/Fig_4_predicted_interactions_independent_axes_new_09_11_2025.pdf", 
       plot = new_plot, width = 17, height = 21, units = "cm", dpi = 600)
