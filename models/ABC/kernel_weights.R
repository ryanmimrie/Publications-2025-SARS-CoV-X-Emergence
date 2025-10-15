# ==============================================================================
# ===== SARS-CoV-X Emergence: Kernel Weights ===================================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following calculates epanechnikov kernel weights from the output of
# ABC_Wuhan.R and calculate_distances.R.

# File paths:
# This script uses the here() library to set paths dynamically. Path errors may occur
# if this script is opened in an already running Rstudio session. To resolve this,
# close Rstudio and reopen it by double-clicking on this file.

# ----- 0.2. Dependecies -------------------------------------------------------

library(tidyverse)
library(HDInterval)
library(patchwork)
library(here)

# ----- 0.3. Load Data ---------------------------------------------------------

data <- read.csv(here("models", "ABC", "distances_Wuhan.csv"))

# ------------------------------------------------------------------------------
# ----- 1. Kernel Weights ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Calculate Kernel Weights ------------------------------------------

threshold_rmse <- 0.01

epanechnikov_kernel <- function(u) {ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)}

data <- filter(data, mean_rmse <= threshold_rmse) %>%
  mutate(scaled_rmse = mean_rmse / threshold_rmse,
         weight = epanechnikov_kernel(scaled_rmse)) %>%
  filter(weight > 0) %>%
  mutate(weight = weight / sum(weight, na.rm = TRUE))

# ----- 1.2. Plot approximate posteriors ---------------------------------------

priors <- data.frame(parameter = NA, type = NA, mean = NA, sd = NA)[-1,]

priors <- rbind(priors, data.frame(parameter = "R0_wuhan", type = "lognormal", mean = 1, sd = 0.5))
priors <- rbind(priors, data.frame(parameter = "R0_alpha", type = "lognormal", mean = 1, sd = 0.5))
priors <- rbind(priors, data.frame(parameter = "R0_delta", type = "lognormal", mean = 1, sd = 0.5))
priors <- rbind(priors, data.frame(parameter = "R0_omicron", type = "lognormal", mean = 1, sd = 0.5))

priors <- rbind(priors, data.frame(parameter = "recovery_wuhan", type = "logit-normal", mean = -1, sd = 0.6))
priors <- rbind(priors, data.frame(parameter = "recovery_alpha", type = "logit-normal", mean = -1, sd = 0.6))
priors <- rbind(priors, data.frame(parameter = "recovery_delta", type = "logit-normal", mean = -1, sd = 0.6))
priors <- rbind(priors, data.frame(parameter = "recovery_omicron", type = "logit-normal", mean = -1, sd = 0.6))

priors <- rbind(priors, data.frame(parameter = "waning_wuhan", type = "logit-normal", mean = -5, sd = 2))
priors <- rbind(priors, data.frame(parameter = "waning_alpha", type = "logit-normal", mean = -5, sd = 2))
priors <- rbind(priors, data.frame(parameter = "waning_delta", type = "logit-normal", mean = -5, sd = 2))
priors <- rbind(priors, data.frame(parameter = "waning_omicron", type = "logit-normal", mean = -5, sd = 2))

priors_discretised <- data.frame()

for (i in 1:nrow(priors)) {
  p <- priors[i, ]
  
  if (p$type == "lognormal") {
    y_vals <- dlnorm(seq(0, 6, 0.01), meanlog = p$mean, sdlog = p$sd)
    priors_discretised <- rbind(priors_discretised,
                                data.frame(parameter = p$parameter,
                                           x = seq(0, 6, 0.01),
                                           y = y_vals))
    
  } else if (p$type == "logit-normal") {
    x_vals_ln <- seq(0.0001, 0.9999, 0.001)
    logit_x <- qlogis(x_vals_ln)
    y_vals <- dnorm(logit_x, mean = p$mean, sd = p$sd) / (x_vals_ln * (1 - x_vals_ln))
    
    priors_discretised <- rbind(priors_discretised,
                                data.frame(parameter = p$parameter,
                                           x = x_vals_ln,
                                           y = y_vals))
  }
}

data_long <- dplyr::select(data, 1:12, weight) %>%
  pivot_longer(cols = -weight, names_to = "parameter", values_to = "value")

data_long$parameter <- factor(data_long$parameter, levels = colnames(data)[1:12])
priors_discretised$parameter <- factor(priors_discretised$parameter, levels = colnames(data)[1:12])

p1 <- ggplot() +
  geom_density(data = filter(data_long, grepl("R0", parameter)), aes(x = value, weight = weight), color = "blue") +
  geom_line(data = filter(priors_discretised, grepl("R0", parameter)), aes(x = x, y = y), linetype = "dotted") +
  scale_x_continuous(limits = c(0, 6)) +
  facet_wrap(~parameter, nrow = 1, scales = "free_y") +
  theme(aspect.ratio = 1,
        axis.title.x = element_blank())

p2 <- ggplot() +
  geom_density(data = filter(data_long, grepl("recovery", parameter)), aes(x = value, weight = weight), color = "blue") +
  geom_line(data = filter(priors_discretised, grepl("recovery", parameter)), aes(x = x, y = y), linetype = "dotted") +
  scale_x_continuous(limits = c(0, 0.5)) +
  facet_wrap(~parameter, nrow = 1, scales = "free_y") +
  theme(aspect.ratio = 1,
        axis.title.x = element_blank())

p3 <- ggplot() +
  geom_density(data = filter(data_long, grepl("waning", parameter)), aes(x = value, weight = weight), color = "blue") +
  geom_line(data = filter(priors_discretised, grepl("waning", parameter)), aes(x = x, y = y), linetype = "dotted") +
  scale_x_continuous(limits = c(0, 0.03)) +
  facet_wrap(~parameter, nrow = 1, scales = "free_y") +
  theme(aspect.ratio = 1)

p1 / p2 / p3

# ----- 1.3. Summary Statistics ------------------------------------------------

summary_stats <- data_long %>% group_by(parameter) %>%
  summarise(mode = {
      d <- density(value, weights = weight)
      d$x[which.max(d$y)]},
    hdi_vals = list(hdi(value, credMass = 0.95, weights = weight))) %>%
  mutate(hdi_lower = map_dbl(hdi_vals, 1), hdi_upper = map_dbl(hdi_vals, 2)) %>%
  select(-hdi_vals)

summary_stats

# ----- 1.4. Write Posterior ---------------------------------------------------

posterior_Wuhan <- dplyr::select(data, R0_wuhan, recovery_wuhan, weight)

write_csv(posterior_Wuhan, "posterior_Wuhan.csv")






