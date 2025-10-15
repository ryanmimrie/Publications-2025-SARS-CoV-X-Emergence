# ==============================================================================
# ===== SARS-X Emergence: Emulator Sample Schemes ==============================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script produces the Sobol and Latin Hypercube parameter sample
# schemes used to produce the training and test datasets for model emulation.

# ----- 0.2. Dependencies ------------------------------------------------------

library(randtoolbox); library(tidyverse); library(lhs); library(here)

# ------------------------------------------------------------------------------
# ----- 1. Sobol Samples -------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Parameter Ranges --------------------------------------------------

R0_Omicron              <- c(2, 6)
R0_SARSX                <- c(2, 6)
vaccine_coverage        <- c(0, 1)
vaccine_timing          <- c(-360, 360)
vaccine_cross_immunity  <- c(0, 1)
natural_cross_immunity  <- c(0, 1)
waning_rate_Omicron     <- c(0, 0.02)
waning_rate_SARSX       <- c(0, 0.02)
waning_rate_vaccine     <- c(0, 0.02)

# ----- 1.2. Sobol Design Scheme -----------------------------------------------

sob <- sobol(n = 500000, dim = 9)

design <- data.frame(R0_Omicron             = sob[,1] * diff(R0_Omicron) + R0_Omicron[1],
                     R0_SARSX               = sob[,2] * diff(R0_SARSX) + R0_SARSX[1],
                     vaccine_coverage       = sob[,3] * diff(vaccine_coverage) + vaccine_coverage[1],
                     vaccine_timing         = sob[,4] * diff(vaccine_timing) + vaccine_timing[1],
                     vaccine_cross_immunity = sob[,5] * diff(vaccine_cross_immunity) + vaccine_cross_immunity[1],
                     natural_cross_immunity = sob[,6] * diff(natural_cross_immunity) + natural_cross_immunity[1],
                     waning_rate_Omicron    = sob[,7] * diff(waning_rate_Omicron) + waning_rate_Omicron[1],
                     waning_rate_SARSX      = sob[,8] * diff(waning_rate_SARSX) + waning_rate_SARSX[1],
                     waning_rate_vaccine    = sob[,9] * diff(waning_rate_vaccine) + waning_rate_vaccine[1])

# ----- 1.3. Write Samples in Batches ------------------------------------------

write_csv(design[1:100000,], here("models", "emulator", "sample_1.csv"))
write_csv(design[100001:200000,], here("models", "emulator", "sample_2.csv"))
write_csv(design[200001:300000,], here("models", "emulator", "sample_3.csv"))
write_csv(design[300001:400000,], here("models", "emulator", "sample_4.csv"))
write_csv(design[400001:500000,], here("models", "emulator", "sample_5.csv"))


# ------------------------------------------------------------------------------
# ----- 2. LHC Samples ---------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 LHC Design Scheme --------------------------------------------------

lhs_sample <- randomLHS(n = 5000, k = 9)

test_design <- data.frame(R0_Omicron             = lhs_sample[,1] * diff(R0_Omicron) + R0_Omicron[1],
                          R0_SARSX               = lhs_sample[,2] * diff(R0_SARSX) + R0_SARSX[1],
                          vaccine_coverage       = lhs_sample[,3] * diff(vaccine_coverage) + vaccine_coverage[1],
                          vaccine_timing         = lhs_sample[,4] * diff(vaccine_timing) + vaccine_timing[1],
                          vaccine_cross_immunity = lhs_sample[,5] * diff(vaccine_cross_immunity) + vaccine_cross_immunity[1],
                          natural_cross_immunity = lhs_sample[,6] * diff(natural_cross_immunity) + natural_cross_immunity[1],
                          waning_rate_Omicron    = lhs_sample[,7] * diff(waning_rate_Omicron) + waning_rate_Omicron[1],
                          waning_rate_SARSX      = lhs_sample[,8] * diff(waning_rate_SARSX) + waning_rate_SARSX[1],
                          waning_rate_vaccine    = lhs_sample[,9] * diff(waning_rate_vaccine) + waning_rate_vaccine[1])

# ----- 2.2 Write Samples ------------------------------------------------------

write_csv(test_design, here("models", "emulator", "sample_test_1.csv"))


