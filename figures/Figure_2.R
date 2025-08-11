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

data2a <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_2_output_1"), pattern = "\\.csv$", full.names = TRUE), read.csv))
data2b <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_2_output_2"), pattern = "\\.csv$", full.names = TRUE), read.csv))
data2c <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_2_output_3"), pattern = "\\.csv$", full.names = TRUE), read.csv))

# ------------------------------------------------------------------------------
# ----- 1. Data Wrangling ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Data Wrangling ----------------------------------------------------

data2 <- data2a
data2$emergences <- data2$emergences + data2b$emergences + data2c$emergences
data2$trials <- data2$trials + data2b$trials + data2c$trials

data2$prob <- data2$emergences / data2$trials

data2_bg <- filter(data2, vaccine_coverage == 0) %>%
  group_by(virus) %>%
  summarise(background = mean(prob))

data2 <- left_join(data2, data2_bg)

data2$change <- pmin(0, data2$prob - data2$background)

data2 <- left_join(data2, data.frame(virus = 1:4,
                                     strain = c("RaTG13 (bat)", "Rs4084 (bat)", "GX/P1E (pangolin)", "SARS-CoV")))

data2$strain <- factor(data2$strain, levels = c("SARS-CoV", "Rs4084 (bat)", "GX/P1E (pangolin)", "RaTG13 (bat)"))

data2_bg <- left_join(data2_bg, data.frame(virus = 1:4,
                                     strain = c("RaTG13 (bat)", "Rs4084 (bat)", "GX/P1E (pangolin)", "SARS-CoV")))

# ------------------------------------------------------------------------------
# ----- 2. Plots ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Figure 2A ---------------------------------------------------------

palette <- magma(12)

palette[12] <- "#fdfdcf" # Soften top colour (25% closer to white)

p1 <- ggplot(data2) +
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

ggplot(data2_bg) +
  geom_tile(aes(x = 1, y = 1, fill = background)) +
  facet_wrap(~strain) +
  scale_fill_viridis_c(limits = c(0, 0.07))




