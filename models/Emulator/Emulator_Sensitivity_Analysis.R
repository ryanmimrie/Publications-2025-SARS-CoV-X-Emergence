# ==============================================================================
# ===== SARS-X Emergence: Sobol Sensitivity Analysis ===========================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script runs Sobol sensitivity analysis on 9 parameters, with
# probability of emergence estimates provided from the GBRT emulator.

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse); library(here); library(sensitivity)

# ----- 0.3. Load Data ---------------------------------------------------------

load(here("models", "emulator", "emulator.Rdata"))

# ------------------------------------------------------------------------------
# ----- 1. Sobol Sensitivity Analysis ------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Parameter Ranges --------------------------------------------------

params <- c("R0_Omicron",
            "R0_SARSX",
            "vaccine_coverage",
            "vaccine_timing",
            "vaccine_cross_immunity",
            "natural_cross_immunity",
            "waning_rate_Omicron",
            "waning_rate_SARSX",
            "waning_rate_vaccine")

mins <- c(2, 2, 0, -90, 0, 0, 0, 0, 0)
maxs <- c(6, 6, 1,  90, 1, 1, 0.02, 0.02, 0.02)
ranges <- maxs - mins

# ----- 1.2. Sobol Design ------------------------------------------------------

n <- 200000
d <- length(params)

n * (d + 2) # Total sample size

X1 <- data.frame(matrix(runif(d * n), nrow = n))
X2 <- data.frame(matrix(runif(d * n), nrow = n))

scale_to_range <- function(X) {sweep(X, 2, ranges, "*") + matrix(rep(mins, each = nrow(X)), nrow = nrow(X))}

X1_scaled <- scale_to_range(X1)
X2_scaled <- scale_to_range(X2)

x <- sobol2002(model = NULL, X1_scaled, X2_scaled, nboot = 1000)

data_sobol <- as.data.frame(x$X)
colnames(data_sobol) <- params

data_sobol <- as.matrix(data_sobol)

# ----- 1.3. Get Emulator Predictions ------------------------------------------

y <- predict(xgb_model, newdata = data_sobol, iterationrange = c(1, xgb_model$best_iteration))

# ----- 1.4. Calculate Sobol Indices -------------------------------------------

tell(x,y)

data_soboltable <- data.frame(parameter = params,
                              first_order = x$S$original,
                              first_order_low = x$S$`min. c.i.`,
                              first_order_high = x$S$`max. c.i.`,
                              total_order = x$`T`$original,
                              total_order_low = x$`T`$`min. c.i.`,
                              total_order_high = x$`T`$`max. c.i.`)

data_soboltable <- data_soboltable %>%
  mutate(across(where(is.numeric), ~ round(pmax(., 0), 3))) %>%
  arrange(desc(first_order))

data_soboltable
