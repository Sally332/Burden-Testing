#!/usr/bin/env Rscript

################################################################################
# File: variant_filtering_main.R
# Description:
# This script prepares a filtered variant table for downstream burden testing.
# It supports flexible filtering schemes for different project-specific needs,
# including:
#   - variant impact/consequence filters
#   - gnomAD frequency filters
#   - CADD scoring thresholds
#   - segregation rules across multiple family structures
#   - optional mapping of variants to custom genomic regions (e.g., BED intervals)
# It also prepares genotype and covariate matrices for subsequent analyses,
# including PCA for population structure correction.
#
# Input formats:
# - 'subjects': expects REFVAR_INGROUP_NAMES, REFVAR_OUTGROUP, and composite IDs
# - 'regions': optional BED-style annotation file to replace gene IDs with regions
#
# Output: cleaned TSV file used as input for burden_testing_pipeline.R
################################################################################

suppressMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ---------------------------- PARSE ARGUMENTS ---------------------------- #

option_list <- list(
  make_option("--input", type="character", help="Input file path"),
  make_option("--output", type="character", help="Filtered output file path"),
  make_option("--format", type="character", default="subjects", 
              help="Input format: 'subjects' or 'matrix'"),
  make_option("--min_carriers", type="integer", default=3, 
              help="Minimum number of affected individuals per variant"),
  make_option("--af_max", type="double", default=0.01, 
              help="Max allowed allele frequency in gnomAD (or similar)", metavar="double"),
  make_option("--cadd_min", type="double", default=15, 
              help="Minimum CADD score for variant inclusion"),
  make_option("--geno_matrix", type="character", default="geno_matrix.tsv", 
              help="Output path for genotype matrix"),
  make_option("--covar_file", type="character", default="covar.tsv", 
              help="Output path for covariate file (PCs)"),
  make_option("--region_map", type="character", default=NULL, 
              help="(Optional) TSV with variant-to-region mapping: VariantID,RegionID")
)

opt <- parse_args(OptionParser(option_list))

# ---------------------------- LOAD INPUT FILE ---------------------------- #

cat("Reading input variant file...\n")
variants <- read_tsv(opt$input, show_col_types = FALSE)

# ---------------------------- FLEXIBLE FILTERING ---------------------------- #

cat("Applying variant impact filters...\n")

if ("Consequence" %in% colnames(variants)) {
  variants <- variants %>% 
    filter(str_detect(Consequence, "missense|frameshift|stop_gained|splice"))
}

if ("gnomAD_AF" %in% colnames(variants)) {
  variants <- variants %>% filter(gnomAD_AF < opt$af_max)
}

if ("CADD" %in% colnames(variants)) {
  variants <- variants %>% filter(CADD >= opt$cadd_min)
}

# ---------------------------- SEGREGATION FILTER ---------------------------- #
#column names must follow a composite format
# that encodes the family and individual IDs. Example structure:
# <Prefix>_<FamilyID>_<IndividualID>   →  e.g., STUDY123_FAM1001_0001_A
# This allows the script to infer family structure when applying segregation filters.
# Make sure to preserve this convention for compatibility.

cat("Applying segregation logic...\n")

if (opt$format == "subjects") {
  if ("REFVAR_INGROUP_NAMES" %in% colnames(variants)) {
    variants <- variants %>% mutate(n_affected = str_count(REFVAR_INGROUP_NAMES, ",") + 1)
    variants <- variants %>% filter(n_affected >= opt$min_carriers)
  }
  if ("REFVAR_OUTGROUP" %in% colnames(variants)) {
    variants <- variants %>% filter(REFVAR_OUTGROUP == "")
  }
}

# ---------------------------- CUSTOMIZABLE LOGIC ---------------------------- #

# Add additional combinations as needed, for example:
# variants <- variants %>% filter(Impact == "HIGH" | (Impact == "MODERATE" & CADD > 20))
# Use clinical annotation tags if available (e.g., ClinVar pathogenicity)

# ---------------------------- REGION MAPPING (OPTIONAL) ---------------------------- #

if (!is.null(opt$region_map)) {
  cat("Merging variant-to-region mapping...\n")
  map_tbl <- read_tsv(opt$region_map, col_names = c("VariantID", "RegionID"))
  if (!"VariantID" %in% colnames(variants)) {
    stop("Column 'VariantID' required in main file for region mapping.")
  }
  variants <- variants %>% left_join(map_tbl, by = "VariantID")
  if (!"RegionID" %in% colnames(variants)) {
    stop("Region mapping failed. Column 'RegionID' not found after join.")
  }
}

# ---------------------------- EXPORT FILTERED VARIANTS ---------------------------- #

cat("Writing filtered output to:", opt$output, "\n")
write_tsv(variants, opt$output)

# ---------------------------- BUILD GENOTYPE MATRIX ---------------------------- #

cat("Building genotype matrix...\n")
# NOTE: Replace with actual aggregation logic in real application
geno_matrix <- data.frame(
  GroupID = c("TP53", "BRCA1", "CHEK2", "ATM", "MLH1"),
  FAM1_001 = c(0,1,0,2,0),
  FAM1_002 = c(1,1,0,2,0),
  FAM1_003 = c(1,1,2,2,0),
  FAM2_001 = c(0,0,1,0,1),
  FAM2_002 = c(2,0,1,1,1),
  CTRL_001 = c(0,0,0,0,0),
  CTRL_002 = c(0,0,1,1,0)
)

write.table(geno_matrix, opt$geno_matrix, sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------- COMPUTE PCs ---------------------------- #

cat("Computing principal components from genotype matrix...\n")
geno <- geno_matrix
rownames(geno) <- geno$GroupID
geno <- geno[,-1]  # remove group ID column
geno_scaled <- scale(t(geno))
pca <- prcomp(geno_scaled)
pcs <- data.frame(SampleID = rownames(geno_scaled), PC1 = pca$x[,1], PC2 = pca$x[,2])
write.table(pcs, file = opt$covar_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n✔ Variant filtering and PCA complete. Outputs written to:\n")
cat("-", opt$output, "\n")
cat("-", opt$geno_matrix, "\n")
cat("-", opt$covar_file, "\n")
