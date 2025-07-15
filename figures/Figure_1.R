library(readxl)
library(tidyverse)
library(lme4)
library(merTools)
library(emmeans)
library(gridExtra)
library(here)

process_sheet <- function(sheet_name) {
  
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

file_path <- here("data", "OCT2024_ Rs4084_pseudotype fixed template.xlsx")
sheet_names <- excel_sheets(file_path)
data1 <- map_dfr(sheet_names[1:9], process_sheet)

file_path <- here("data", "OCT2024_GXP1E_pseudotype fixed template.xlsx")
sheet_names <- excel_sheets(file_path)
data2 <- map_dfr(sheet_names[1:9], process_sheet)

file_path <- here("data", "OCT2024_RaTG13_pseudotype fixed template.xlsx")
sheet_names <- excel_sheets(file_path)
data3 <- map_dfr(sheet_names[1:9], process_sheet)

file_path <- here("data", "OCT2024_SARS-CoV-1_pseudotype fixed template.xlsx")
sheet_names <- excel_sheets(file_path)
data4 <- map_dfr(sheet_names[1:9], process_sheet)

data <- rbind(data1, data2, data3, data4)

data <- na.omit(data)

data$fluorescence <- as.numeric(data$fluorescence)

#rm(list = setdiff(ls(), "data"))

data_NSC <- filter(data, ID == "NSC")

data_NSC <- data_NSC %>% group_by(plate, virus) %>%
  summarise(fluorescence_NSC = mean(fluorescence))

data <- data %>% filter(!ID %in% c("NSC", "-ve", "+ve")) %>%
  group_by(ID, plate, virus) %>%
  summarise(fluorescence_mean = mean(fluorescence))


data <- left_join(data, data_NSC)

data$ID <- ifelse(data$ID == "38803", "RG38803", data$ID)

data$neutralisation <- 100 - (data$fluorescence_mean / data$fluorescence_NSC) * 100

data$neutralisation <- ifelse(data$neutralisation < 0, 0, data$neutralisation)

metadata <- read_csv("2023_SARSCoV2_data.csv")[,1:3]

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
  scale_fill_manual(values = c("#e74c3c", "#f1c40f", "#1abc9c", "#3498db"), labels = c("SARS-CoV", "GX/P1E (pangolin)", "Rs4084 (bat)", "RaTG13 (bat)")) +
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
  scale_fill_manual(values = c("#e74c3c", "#f1c40f", "#1abc9c", "#3498db"), labels = c("SARS-CoV", "GX/P1E (pangolin)", "Rs4084 (bat)", "RaTG13 (bat)")) +
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
  scale_fill_manual(values = c("#e74c3c", "#f1c40f", "#1abc9c", "#3498db"), labels = c("SARS-CoV", "GX/P1E (pangolin)", "Rs4084 (bat)", "RaTG13 (bat)")) +
  ggtitle("C.") +
  coord_fixed(ratio = 0.0375)

data$cat <- ifelse(data$history == "Naive" & data$Doses != "Unvaccinated", "Vaccinated",
                   ifelse(data$history != "Naive" & data$Doses == "Unvaccinated", "Infected",
                          ifelse(data$history == "Naive" & data$Doses == "Unvaccinated", "Naive", "Infected + Vaccinated")))

data2 <- filter(data, cat != "Naive")
      
data2$cat <- factor(data2$cat, levels = c("Vaccinated", "Infected", "Infected + Vaccinated"))

data2 <- left_join(data2, data.frame(virus = as.character(unique(data$virus)),
                                   similarity = c(75.9, 97.4, 76.9, 92.2)))  
              
data2$simcat <- as.character(data2$similarity)

ggplot(data2) +
  geom_boxplot(aes(x = simcat, y = neutralisation))

data2$sim_zeroed <- data2$similarity-75.9


data_simplified <- data

data_simplified$history <- as.character(data_simplified$history)

data_simplified$history <- ifelse(data_simplified$history == "Naive", "Naive", "Recovered")


combined_se <- function(se1, se2, corr) {
  # Calculate the variances
  var1 <- se1^2
  var2 <- se2^2
  
  # Calculate the covariance
  cov <- corr * se1 * se2
  
  # Calculate the combined variance
  combined_var <- var1 + var2 + 2 * cov
  
  # Return the combined standard error
  return(sqrt(combined_var))
}

library(lmerTest)

model <- lmer(neutralisation ~ history * Doses + (1 | virus) + (1 | ID), data = data_simplified)

summary(model)

emmeans(model, ~ history * Doses)

model <- lmer(neutralisation ~
                virus + (1 | history) + (1 | Doses) +
                (1 | history:Doses) + (1 | ID), data = data_simplified)

summary(model)

emmeans(model, ~ virus)





model <- lmer(neutralisation ~ I(sim_zeroed^2) + (1 | history) + (1 | Doses) + (1 | history:Doses) + (1 | ID), data = data2)

summary(model)

emmeans(model, ~ I(sim_zeroed^2))

model_null <- lmer(neutralisation ~ 1 + (1 | history) + (1 | Doses) + (1 | history:Doses) + (1 | ID), data = data2)

anova(model, model_null)




preds <- predictInterval(model, newdata = data2, which = "fixed", n.sims = 1000, level = 0.68)

# Add predictions to data
data2 <- cbind(data2, preds)

# Summarise by similarity
data2_sum <- data2 %>%
  group_by(similarity) %>%
  summarise(neut.mean = mean(fit),
            neut.lower = mean(lwr),
            neut.upper = mean(upr))

ggplot(data2_sum) +
  geom_smooth(aes(x = similarity, y = neut.lower), method = "lm", formula = y ~ poly(x, 2), se = F) + # lower se
  geom_smooth(aes(x = similarity, y = neut.upper), method = "lm", formula = y ~ poly(x, 2), se = F) + # upper se
  geom_smooth(aes(x = similarity, y = neut.mean), method = "lm", formula = y ~ poly(x, 2), se = F) + # mean line
  labs(y = "Neutralisation (%)")


similarity_seq <- seq(75, 100, length.out = 200)

# Fit polynomial models to predict mean, lower, and upper
mean_model <- lm(neut.mean ~ poly(similarity, 2), data = data2_sum)
lower_model <- lm(neut.lower ~ poly(similarity, 2), data = data2_sum)
upper_model <- lm(neut.upper ~ poly(similarity, 2), data = data2_sum)

predictions <- data.frame(
  similarity = similarity_seq,
  neut.mean = predict(mean_model, newdata = data.frame(similarity = similarity_seq)),
  neut.lower = predict(lower_model, newdata = data.frame(similarity = similarity_seq)),
  neut.upper = predict(upper_model, newdata = data.frame(similarity = similarity_seq))
)

predictions$neut.upper <- ifelse(predictions$neut.upper > 100, 100, predictions$neut.upper)

data2_sum <- data2 %>% group_by(similarity) %>%
  summarise(neut.mean = mean(neutralisation))

p4 <- ggplot(predictions) +
  geom_point(data = data2, aes(x = similarity, y = neutralisation), position = "jitter", alpha = 0.1) +
  geom_ribbon(aes(x = similarity, ymin = neut.lower, ymax = neut.upper), alpha = 0.2) +
  geom_line(aes(x = similarity, y = neut.mean), color = "red", size = 0.75) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16)) +
  scale_x_continuous(name = "Spike aa similarity (%)") +
  scale_y_continuous(name = "Neutralisation (%)") +
  coord_fixed(ratio = 1/4) +
  ggtitle("E.")

combined <- grid.arrange(
  arrangeGrob(p1, p2, p3, ncol = 1),
  arrangeGrob(p4, nrow = 1, heights = c(0.2)),
  widths = c(10, 2.2)
)


ggsave("~/Research/Odin/Plots/Figure1_raw.png", combined, dpi = 300, width = 16, height = 8)
ggsave("~/Research/Odin/Plots/Figure1_raw.svg", combined, width = 16, height = 8)


# Naive
# Recovered: Wuhan
# Recovered: Alpha
# Recovered: Delta

# Unvaccinated
# Vaccinated: 1 dose
# Vaccinated: 2 doses

tmp <- filter(data, history == "Recovered: Delta", Doses == "Vaccinated: 2 doses")

pairwise.t.test(tmp$neutralisation, tmp$virus, p.adjust.method = "holm", pool.sd = FALSE)

library(ggtree)
library(treeio)
library(ape)
library(phytools)

setwd("~/Research/Odin/Viruses/data")

tree <- read.tree("al_S_aa.fasta.treefile")

tree$tip.label <- c("GX/P1E (pangolin)", "Wuhan-Hu-1", "RaTG13 (bat)", "Rs4084 (bat)", "SARS-CoV")

rooted_tree <- midpoint.root(tree)

rooted_tree$node.label <-  c("", "1", "0.99", "1")

p5 <- ggtree(rooted_tree) +
  geom_tiplab(size = 4.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust = -0.3, size = 3.5) +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0.9, 5.1)) +
  geom_treescale(x = 0, y = 0.75, width = 0.05, offset = 0.075, fontsize = 3.5) +
  ggtitle("D.") +
  theme(text = element_text(size = 13))

ggsave("~/Research/Odin/Plots/Figure1_phylo.png", p5, dpi = 300, width = 4, height = 3)
ggsave("~/Research/Odin/Plots/Figure1_phylo.svg", p5, width = 4, height = 3)

table(data$virus, data$history, data$Doses)
