library(tidyverse); library(here)

data2 <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_2_output"), pattern = "\\.csv$", full.names = TRUE), read.csv))

data2 <- left_join(data2, data.frame(virus = 1:4,
                                     strain = c("RaTG13", "Rs4084", "GX/P1E", "SARS-1")))

data2$strain <- factor(data2$strain, levels = c("SARS-1", "Rs4084", "GX/P1E", "RaTG13"))

ggplot(data2) +
  geom_tile(aes(x = vaccine_timing, y = vaccine_coverage, fill = emergences/trials)) +
  facet_wrap(~strain) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1)

ggplot(filter(data2, vaccine_coverage == 1)) +
  geom_line(aes(x = vaccine_timing, y = emergences/trials, color = strain, group = strain)) +
  theme_bw() +
  theme(aspect.ratio = 1)


data3a <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_3_output_1"), pattern = "\\.csv$", full.names = TRUE), read.csv))
data3b <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_3_output_2"), pattern = "\\.csv$", full.names = TRUE), read.csv))
data3c <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_3_output_3"), pattern = "\\.csv$", full.names = TRUE), read.csv))
data3d <- do.call(rbind, lapply(list.files(path = here("models", "figures", "figure_3_output_4"), pattern = "\\.csv$", full.names = TRUE), read.csv))

data3 <- data3a
data3$emergences <- data3$emergences + data3b$emergences + data3c$emergences + data3d$emergences
data3$trials <- data3$trials + data3b$trials + data3c$trials + data3d$trials


ggplot(data3) +
  geom_tile(aes(x = vaccine_timing, y = vaccine_coverage, fill = emergences/trials)) +
  facet_grid(cols = vars(natural_cross_immunity), rows = vars(R0)) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1)




data4a <- read_csv(here("models", "figures", "figure_4_output_1.csv"))
data4b <- read_csv(here("models", "figures", "figure_4_output_2.csv"))

data4 <- data4a
data4$emergences <- data4$emergences + data4b$emergences
data4$trials <- data4$trials + data4b$trials

data4$prob <- data4$emergences / data4$trials

data4_bg <- filter(data4, vaccine_coverage == 0) %>%
  group_by(R0, natural_cross_immunity) %>%
  summarise(background = mean(prob))

data4 <- left_join(data4, data4_bg)

data4$change <- pmax(0, data4$prob - data4$background)

ggplot(data4) +
  geom_tile(aes(x = vaccine_timing, y = vaccine_coverage, fill = change)) +
  facet_grid(cols = vars(natural_cross_immunity), rows = vars(R0)) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1)
