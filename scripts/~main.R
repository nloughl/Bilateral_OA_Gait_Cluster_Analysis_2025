# Main -- runs the scripts in order (not including visualization .rmd files, run those separately) 

rm(list = ls())

# libraries 
library(here)


# List all needed R scripts
script_files <- c(
  "Bilateral_v3d_data2023_preprocessing.R",
  "Bilateral_RC_MWF_merge_2023.R",
  "Bilateral_RC_data2025_preprocessing.R",
  "Bilateral_v3d_data2025_preprocessing.R",
  "Bilateral_merge_RC_v3d_2025.R",
  "Bilateral_one_leg_2023.R",
  "Bilateral_one_leg_2025.R",
  "Bilateral_merge_2023_2025.R",
  "Bilateral_PCA.r",
  "Bilateral_MDS.r",
  "Bilateral_HC_Ward.r" 
)

# Run each script in order
log_file <- "error_log.txt"

for (script in script_files) {
  scripts_dir <- here::here("scripts")  # Set the path to the scripts folder
  script_path <- file.path(scripts_dir, script)
  
  if (file.exists(script_path)) {
    cat("Running:", script_path, "\n")
    
    tryCatch(
      source(script_path),
      error = function(e) {
        msg <- paste(Sys.time(), "- Error in", script_path, ":", e$message, "\n")
        cat(msg, file = log_file, append = TRUE)
      }
    )
  } else {
    warning("File not found:", script_path)
  }
}


