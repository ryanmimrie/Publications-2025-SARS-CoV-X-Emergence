library(reticulate)
library(here)

# Source the Python script
source_python(here("models", "ABC", "calculate_distances.py"))

process_files(here("models", "ABC", "ABC_samples_Wuhan"),
              file_prefix = "1",
              data = here("models", "ABC", "prevalences.csv"))


