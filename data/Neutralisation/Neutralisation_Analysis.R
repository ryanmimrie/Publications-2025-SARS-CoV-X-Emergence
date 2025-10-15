# ==============================================================================
# ===== SARS-CoV-X Emergence Modeling: Neutralisation Analysis =================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script wrangles, analyses and plots Figure 1 from Imrie & Bissett
# et al., (2025) "Post-pandemic changes in population immunity have reduced the
# likelihood of emergence of zoonotic coronaviruses".

# ----- 0.2. Dependencies ------------------------------------------------------

library(readxl)
library(tidyverse)
library(lme4)
library(merTools)
library(emmeans)
library(gridExtra)
library(here)
library(viridisLite)
library(patchwork)

# ----- 0.3. Load Data Custom Function -----------------------------------------

get_data_from_sheets <- function(sheet_name) {
  
  sheet_data <- read_excel(file_path, sheet = sheet_name, range = "B1:M18")
  
  plate1 <- sheet_data[1:8, 1:12] %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "column", values_to = "plate1_value") %>%
    mutate(row = rep(1:8, each = 12))
  
  plate2 <- sheet_data[10:17, 1:12] %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "column", values_to = "plate2_value") %>%
    mutate(row = rep(1:8, each = 12))
  
  plate1 <- dplyr::select(plate1, plate1_value)
  plate2 <- dplyr::select(plate2, plate2_value)
  
  colnames(plate1) <- "ID"
  colnames(plate2) <- "fluorescence"
  
  combined <- cbind(plate1, plate2)
  
  combined$sheet <- sheet_name
  
  combined$sheet <- str_remove_all(combined$sheet, "Plate ")
  combined <- combined %>%
    separate(sheet, into = c("plate", "virus"), sep = " ", extra = "merge") %>%
    mutate(virus = str_trim(virus, side = "right"))

  return(combined)
}

# ----- 0.4. Load Data ---------------------------------------------------------

file_path <- here("data", "Neutralisation", "OCT2024_ Rs4084_pseudotype fixed template.xlsx")
sheet_names <- excel_sheets(file_path)
data1 <- map_dfr(sheet_names[1:9], get_data_from_sheets)

file_path <- here("data", "Neutralisation", "OCT2024_GXP1E_pseudotype fixed template.xlsx")
sheet_names <- excel_sheets(file_path)
data2 <- map_dfr(sheet_names[1:9], get_data_from_sheets)

file_path <- here("data", "Neutralisation", "OCT2024_RaTG13_pseudotype fixed template.xlsx")
sheet_names <- excel_sheets(file_path)
data3 <- map_dfr(sheet_names[1:9], get_data_from_sheets)

file_path <- here("data", "Neutralisation", "OCT2024_SARS-CoV-1_pseudotype fixed template.xlsx")
sheet_names <- excel_sheets(file_path)
data4 <- map_dfr(sheet_names[1:9], get_data_from_sheets)

data <- rbind(data1, data2, data3, data4) %>% na.omit()

data$ID <- ifelse(data$ID == "38803", "RG38803", data$ID) # Fix typo

# ------------------------------------------------------------------------------
# ----- 1. Data Wrangling ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Calculate Neutralisation % ----------------------------------------

data$fluorescence <- as.numeric(data$fluorescence)

data_NSC <- filter(data, ID == "NSC")

data_NSC <- data_NSC %>% group_by(plate, virus) %>%
  summarise(fluorescence_NSC = mean(fluorescence))

data <- data %>% filter(!ID %in% c("NSC", "-ve", "+ve")) %>%
  group_by(ID, plate, virus) %>%
  summarise(fluorescence_mean = mean(fluorescence))

data <- left_join(data, data_NSC)

data$neutralisation <- 100 - (data$fluorescence_mean / data$fluorescence_NSC) * 100

data$neutralisation <- ifelse(data$neutralisation < 0, 0, data$neutralisation)

# ----- 1.2. Add Sample Metadata -----------------------------------------------

metadata <- read_csv(here("data", "Neutralisation", "2023_SARSCoV2_data.csv"))[,1:3]

data$virus <- ifelse(data$virus == "CoV-1", "SARS-CoV-1",
                     ifelse(data$virus == "pangolin", "GX/P1E (pangolin)",
                            ifelse(data$virus == "RaTG13", "RaTG13 (bat)", "Rs4084 (bat)")))

data$virus <- factor(data$virus, levels = c("SARS-CoV-1", "Rs4084 (bat)", "GX/P1E (pangolin)", "RaTG13 (bat)"))

data <- left_join(data, metadata)

data <- na.omit(data)

data$Doses <- ifelse(data$Doses == 0, "Unvaccinated",
                         ifelse(data$Doses == 1, "Vaccinated: 1 dose", "Vaccinated: 2 doses"))

data$history <- ifelse(grepl("I_W", data$Group), "Recovered: Wuhan",
                       ifelse(grepl("I_A", data$Group), "Recovered: Alpha",
                              ifelse(grepl("I_D", data$Group), "Recovered: Delta", "Naive")))

data$history <- factor(data$history , levels = c("Naive", "Recovered: Wuhan", "Recovered: Alpha", "Recovered: Delta"))

# ------------------------------------------------------------------------------
# ----- 2. Figure 1 ------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 Plots --------------------------------------------------------------

p1 <- ggplot(filter(data, Doses == "Unvaccinated")) +
  geom_point(aes(x = virus, y = neutralisation), position = position_jitter(width = 0.1, height = 0), alpha = 0.25) +
  geom_boxplot(aes(x = virus, y = neutralisation, fill = virus), alpha = 0.75, outlier.shape = NA) +
  facet_grid(cols = vars(history), rows = vars(Doses)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) +
  scale_x_discrete(labels = c("SARS-CoV", "GX/P1E\n(pangolin)", "Rs4084\n(bat)", "RaTG13\n(bat)")) +
  scale_y_continuous(name = "Neutralisation (%)", limits = c(0, 100)) +
  scale_fill_manual(values = mako(5)[2:5], labels = c("SARS-CoV", "GX/P1E (pangolin)", "Rs4084 (bat)", "RaTG13 (bat)")) +
  ggtitle("A.") +
  coord_fixed(ratio = 0.0375)

p2 <- ggplot(filter(data, Doses == "Vaccinated: 1 dose")) +
  geom_point(aes(x = virus, y = neutralisation), position = position_jitter(width = 0.1, height = 0), alpha = 0.25) +
  geom_boxplot(aes(x = virus, y = neutralisation, fill = virus), alpha = 0.75, outlier.shape = NA) +
  facet_grid(cols = vars(history), rows = vars(Doses)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) +
  scale_x_discrete(labels = c("SARS-CoV", "GX/P1E\n(pangolin)", "Rs4084\n(bat)", "RaTG13\n(bat)")) +
  scale_y_continuous(name = "Neutralisation (%)", limits = c(0, 100)) +
  scale_fill_manual(values = mako(5)[2:5], labels = c("SARS-CoV", "GX/P1E (pangolin)", "Rs4084 (bat)", "RaTG13 (bat)")) +
  ggtitle("B.") +
  coord_fixed(ratio = 0.0375)

p3 <- ggplot(filter(data, Doses == "Vaccinated: 2 doses")) +
  geom_point(aes(x = virus, y = neutralisation), position = position_jitter(width = 0.1, height = 0), alpha = 0.25) +
  geom_boxplot(aes(x = virus, y = neutralisation, fill = virus), alpha = 0.75, outlier.shape = NA) +
  facet_grid(cols = vars(history), rows = vars(Doses)) +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) +
  scale_x_discrete(labels = c("SARS-CoV", "GX/P1E\n(pangolin)", "Rs4084\n(bat)", "RaTG13\n(bat)")) +
  scale_y_continuous(name = "Neutralisation (%)", limits = c(0, 100)) +
  scale_fill_manual(values = mako(5)[2:5], labels = c("SARS-CoV", "GX/P1E (pangolin)", "Rs4084 (bat)", "RaTG13 (bat)")) +
  ggtitle("C.") +
  coord_fixed(ratio = 0.0375)

p1 / p2 / p3


