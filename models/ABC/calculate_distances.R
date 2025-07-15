# ==============================================================================
# ===== SARS-CoV-X Emergence: ABC Distance Calculation =========================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following calls the python script calculate_distances.py to calculate
# mean RMSE across stochastic replicates of model iterations used in the ABC
# calibation of SARS-CoV-2 R0, recovery rates, and waning rates.

# ----- 0.2. Dependencies ------------------------------------------------------

library(reticulate)
library(here)

source_python(here("models", "ABC", "calculate_distances.py"))

# ------------------------------------------------------------------------------
# ----- 1. Calculate distances -------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Calculate distances by file prefix --------------------------------

process_files(here("models", "ABC", "ABC_samples_Wuhan"),
              file_prefix = "1",
              data = here("models", "ABC", "prevalences.csv"))

# ----- 1.2. Combine distances from parallel runs ------------------------------

dist_files <- list.files(path = here("models", "ABC", "ABC_samples_Wuhan"),  pattern = "^dist", full.names = TRUE)

all_dists <- do.call(rbind, lapply(dist_files, read.csv))

write.csv(all_dists, file = here("models", "ABC", "distances_Wuhan.csv"), row.names = FALSE)



