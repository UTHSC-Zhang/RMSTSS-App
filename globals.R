# Specify the directory where your R scripts are located
script_directory <- "R/"

# Get a list of all .R files in the directory
r_files <- list.files(
  path = script_directory, 
  pattern = "\\.R$", 
  full.names = TRUE
)

# Loop through the file list and source each one
for (file in r_files) {
  source(file)
}

# Optional: Print a message to confirm
cat("All R scripts in the '", script_directory, "' directory have been sourced.\n", sep = "")
