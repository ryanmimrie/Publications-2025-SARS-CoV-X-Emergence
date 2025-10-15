# ==============================================================================
# ===== SARS-X Emergence: Emulator Training ====================================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script trains the XGBoost GBRT model emulator using the training
# and test datasets produced in the script Emulator_Sample_Schemes.R

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse); library(here); library(xgboost)

# ----- 0.3. Load Data ---------------------------------------------------------

data_training <- read_csv(here("models", "emulator", "training_1.csv"))
data_training <- rbind(data_training, read_csv(here("models", "emulator", "training_2.csv")))
data_training <- rbind(data_training, read_csv(here("models", "emulator", "training_3.csv")))
data_training <- rbind(data_training, read_csv(here("models", "emulator", "training_4.csv")))
data_training <- rbind(data_training, read_csv(here("models", "emulator", "training_5.csv")))

data_training$prob <- data_training$emergences / data_training$trials

data_test <- read_csv(here("models", "emulator", "test_1.csv"))

data_test$prob <- data_test$emergences / data_test$trials

# ------------------------------------------------------------------------------
# ----- 1. XGBoost -------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Wrangle XGBoost Matrices ------------------------------------------

x_train <- as.matrix(data_training[,1:9])

y_train <- data_training$prob

d_train <- xgb.DMatrix(data = x_train, label = y_train)


x_test <- as.matrix(data_test[,1:9])

y_test <- data_test$prob

d_test <- xgb.DMatrix(data = x_test, label = y_test)

# ----- 1.2. CV Folds ----------------------------------------------------------

xgb_cv <- xgb.cv(data = d_train,
                 objective = "reg:squarederror",
                 nrounds = 10000,
                 nfold = 5,
                 metric = "rmse",
                 eta = 0.01,
                 max_depth = 9,
                 subsample = 0.7,
                 colsample_bytree = 1.0,
                 verbose = 2,
                 nthread = 9,
                 early_stopping_rounds = 50)

xgb_cv$best_iteration # 5237, use as upper cutoff for training

# ----- 1.3. Model Training ----------------------------------------------------

xgb_model <- xgb.train(data = d_train,
                       nrounds = 5237,
                       early_stopping_rounds = 50,
                       watchlist = list(train = d_train, eval = d_test),
                       objective = "reg:squarederror",
                       eta = 0.01,
                       max_depth = 9,
                       subsample = 0.7,
                       nthread = 11,
                       verbose = 2)

save(xgb_model, file = here("models", "emulator", "emulator.Rdata"))

# ------------------------------------------------------------------------------
# ----- 2. Validation ----------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Correlation Plot --------------------------------------------------

data_test$pred <- predict(xgb_model, newdata = x_test, iterationrange = c(1, xgb_model$best_iteration))

ggplot(data_test) +
  geom_point(aes(x = prob, y = pred), alpha = 0.1) +
  geom_smooth(aes(x = prob, y = pred), method = "lm") +
  geom_abline(linetype = "dashed", color = "black") +
  scale_x_continuous(name = "Probability of Emergence (Model)") +
  scale_y_continuous(name = "Probability of Emergence (Emulator)") +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

cor.test(data_test$pred, data_test$prob)
cor.test(data_test$pred, data_test$prob)$estimate^2

# ----- 2.2. Sanity Check - 2D Plane Plot --------------------------------------

data <- expand.grid(R0_Omicron = 4,
                    R0_SARSX = 4,
                    vaccine_coverage = seq(0, 1, length.out = 100),
                    vaccine_timing = seq(-360, 360, length.out = 100),
                    vaccine_cross_immunity = 1,
                    natural_cross_immunity = 0.5,
                    waning_rate_Omicron = 0.01,
                    waning_rate_SARSX = 0.01,
                    waning_rate_vaccine = 0.01)

x_new <- as.matrix(data)

data$prob <- predict(xgb_model, newdata = x_new, iterationrange = c(1, xgb_model$best_iteration))

ggplot(data) +
  geom_tile(aes(x = vaccine_timing, y = vaccine_coverage, fill = prob)) +
  scale_fill_viridis_c(option = "A") +
  theme(aspect.ratio = 1)

