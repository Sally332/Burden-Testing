
# Familial Variant Filtering Script
# Description: Filters pre-annotated family-based variant tables to retain rare, likely pathogenic variants for burden testing.

# Load necessary libraries
suppressMessages({
  library(dplyr)
  library(readr)
})

# ---------------------------- USER INPUTS ---------------------------- #

# Set project directory
base.dir    <- "path/to/project"        # <- EDIT: Replace with your project path
input.dir   <- file.path(base.dir, "data/families")
output.dir  <- file.path(base.dir, "results/familial_filtered")

# Define filenames of annotated variant tables (one per family)
file.names <- c(
  "FAMILY123_nonsyn_low_outgroup_count.txt",
  "FAMILY456_nonsyn_low_outgroup_count.txt",
  "FAMILY789_nonsyn_low_outgroup_count.txt"
)

# ---------------------------- SCRIPT START ---------------------------- #

# Create output directory if it doesn't exist
if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

# Extract family names from file names
family.names <- sub("_.*", "", file.names)

# Loop over files and apply filtering
for (i in seq_along(file.names)) {
  file.path.input <- file.path(input.dir, file.names[i])
  file.path.output <- file.path(output.dir, paste0(family.names[i], "_filtered.tsv"))

  # Read data
  dat <- read_tsv(file.path.input, show_col_types = FALSE)

  # Check for expected columns
  expected_cols <- c("Gene", "Consequence", "CADD", "gnomAD_AF", "Outgroup_Count")
  missing_cols <- setdiff(expected_cols, names(dat))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns in", file.names[i], ":", paste(missing_cols, collapse = ", ")))
    next
  }

  # Apply filtering criteria (adjustable)
  dat.filtered <- dat %>%
    filter(
      Consequence %in% c("missense_variant", "stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"),
      is.na(gnomAD_AF) | gnomAD_AF < 0.01,
      is.na(CADD) | CADD > 15,
      Outgroup_Count == 0
    )

  # Write filtered result
  write_tsv(dat.filtered, file.path.output)

  cat("Filtered", file.names[i], "->", file.path.output, "
")
}
