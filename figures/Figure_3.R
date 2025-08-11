# ==============================================================================
# ===== SARS-CoV-X Emergence Modeling: Figure 2 ================================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script wrangles and plots Figure 2 from Imrie & Bissett et al.,
# (2025) "Post-pandemic changes in population immunity have reduced the
# likelihood of emergence of zoonotic coronaviruses".

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(here)
library(scales)
library(viridis)
library(viridisLite)

# ----- 0.3. Load Data ---------------------------------------------------------

# NTS: Collapse to single csv before publication

dataa <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_2_output_1"), pattern = "\\.csv$", full.names = TRUE), read.csv))
datab <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_2_output_2"), pattern = "\\.csv$", full.names = TRUE), read.csv))
datac <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_2_output_3"), pattern = "\\.csv$", full.names = TRUE), read.csv))

# ------------------------------------------------------------------------------
# ----- 1. Data Wrangling ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Data Wrangling ----------------------------------------------------

data <- dataa
data$emergences <- data$emergences + datab$emergences + datac$emergences
data$trials <- data$trials + datab$trials + datac$trials

data$prob <- data$emergences / data$trials

data_bg <- filter(data, vaccine_coverage == 0) %>%
  group_by(virus) %>%
  summarise(background = mean(prob))

data <- left_join(data, data_bg)

data$change <- pmin(0, data$prob - data$background)

data <- left_join(data, data.frame(virus = 1:4,
                                     strain = c("RaTG13 (bat)", "Rs4084 (bat)", "GX/P1E (pangolin)", "SARS-CoV")))

data$strain <- factor(data$strain, levels = c("SARS-CoV", "Rs4084 (bat)", "GX/P1E (pangolin)", "RaTG13 (bat)"))

data_bg <- left_join(data_bg, data.frame(virus = 1:4,
                                     strain = c("RaTG13 (bat)", "Rs4084 (bat)", "GX/P1E (pangolin)", "SARS-CoV")))

# ------------------------------------------------------------------------------
# ----- 2. Plots ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Figure 2A ---------------------------------------------------------

palette <- magma(12)

palette[12] <- "#fdfdcf" # Soften top colour (25% closer to white)

p1 <- ggplot(data) +
  geom_tile(aes(x = vaccine_timing, y = vaccine_coverage, fill = change)) +
  facet_wrap(~strain) +
  scale_fill_gradientn(name = "Change in\nemergence\nprobability", colours = palette, guide = guide_colorbar(title.position = "right")) +
  scale_x_continuous(name = "Vaccine program start (days from exposure event)", expand = c(0,0)) +
  scale_y_continuous(name = "Vaccine program uptake (%)", expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        text = element_text(color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

ggsave(here("figures", "figure 2A raw.svg"), plot = p1, dpi = 300, width = 6, height = 6)

# ----- 2.2. Figure 2B ---------------------------------------------------------

# These colours are added manually to the svg of Figure 2A to properly align
# between subplots. This plot is therefore quite bare.

ggplot(data_bg) +
  geom_tile(aes(x = 1, y = 1, fill = background)) +
  facet_wrap(~strain) +
  scale_fill_viridis_c(limits = c(0, 0.07))




