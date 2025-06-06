#!/usr/bin/env Rscript

################################################################################
# File: burden_testing_pipeline.R
# Description:
# Gene-based burden testing script using SKAT and related models.
# Supports covariate correction, mixed models, QC plotting, and alternative tests.
# Includes SKAT-O, enhanced output with summary statistics and multiple-test corrections.
################################################################################

# ---------------------------- LOAD LIBRARIES ---------------------------- #
# Load required packages for statistical analysis, data handling, and plotting
suppressMessages({
  library(optparse)   # Parse command-line options
  library(SKAT)       # Core burden test library (SKAT, SKAT-O)
  library(data.table) # Fast file reading
  library(dplyr)      # Data manipulation
  library(ggplot2)    # Plotting
  library(reshape2)   # Matrix handling (optional)
  library(qqman)      # For QQ plots
})

# ---------------------------- PARSE ARGUMENTS ---------------------------- #
# Define user-provided arguments for file paths and options
option_list <- list(
  make_option("--geno", type="character", help="Genotype matrix file (gene x sample)"),
  make_option("--pheno", type="character", help="Phenotype table (SampleID + Phenotype)"),
  make_option("--covar", type="character", default=NULL, help="Covariate table (SampleID + PCs etc.)"),
  make_option("--out", type="character", default="burden_results.tsv", help="Output results file"),
  make_option("--plots", type="character", default=NULL, help="Path to save QC plots"),
  make_option("--summary", type="character", default="summary_stats.tsv", help="Summary stats output")
)
opt <- parse_args(OptionParser(option_list))

# ---------------------------- LOAD DATA ---------------------------- #
# Read genotype and phenotype tables; merge with covariates if provided
cat("Reading input files...\n")
geno <- fread(opt$geno, data.table=FALSE, check.names=FALSE, stringsAsFactors=TRUE)
pheno <- fread(opt$pheno, data.table=FALSE)

if (!is.null(opt$covar)) {
  covar <- fread(opt$covar, data.table=FALSE)
  pheno <- merge(pheno, covar, by="SampleID")
}

# Ensure sample order matches across all matrices
rownames(pheno) <- pheno$SampleID
geno <- geno[, colnames(geno) %in% pheno$SampleID]
geno <- geno[, match(pheno$SampleID, colnames(geno))]

# ---------------------------- QC PLOTS ---------------------------- #
# Generate PCA plot, histogram of variant counts, QQ plot, and Volcano plot if paths are provided
if (!is.null(opt$plots) && !is.null(opt$covar)) {
  message("Generating QC plots...")
  dir.create(opt$plots, showWarnings=FALSE, recursive=TRUE)

  # PCA plot
  pc_cols <- grep("PC", names(pheno), value = TRUE)
  if (length(pc_cols) >= 2) {
    p <- ggplot(pheno, aes_string(x=pc_cols[1], y=pc_cols[2], color="as.factor(Phenotype)")) +
      geom_point(size=2) + labs(color="Phenotype") + theme_minimal() +
      ggtitle("PCA: PC1 vs PC2")
    ggsave(filename=file.path(opt$plots, "pca_plot.png"), plot=p)
  }

  # Histogram of variants per gene
  variant_counts <- rowSums(geno > 0)
  hist_data <- data.frame(Gene = rownames(geno), Count = variant_counts)
  h <- ggplot(hist_data, aes(x=Count)) +
    geom_histogram(bins=30, fill="steelblue", color="black") +
    theme_minimal() + xlab("Variant Counts per Gene") + ylab("Frequency")
  ggsave(filename=file.path(opt$plots, "variant_count_histogram.png"), plot=h)
}

# ---------------------------- NULL MODEL ---------------------------- #
# Build null model with or without covariates
cat("Building null model with covariates...\n")
null_formula <- if (!is.null(opt$covar)) {
  as.formula(paste("Phenotype ~", paste(setdiff(names(pheno), c("SampleID", "Phenotype")), collapse=" + ")))
} else {
  as.formula("Phenotype ~ 1")
}
null_model <- SKAT_Null_Model(null_formula, out_type="D", data=pheno, Adjustment=TRUE)

# ---------------------------- BURDEN TESTS ---------------------------- #
# Apply SKAT, SKAT-O, CAST, CMC, and ACAT tests to each gene
cat("Running gene-based burden tests (SKAT, SKAT-O, CAST, CMC, ACAT)...\n")
results <- data.frame()

for (gene in rownames(geno)) {
  gmat <- matrix(as.numeric(geno[gene, ]), nrow=1)
  y <- pheno$Phenotype
  num_vars <- sum(gmat > 0)

  # SKAT and SKAT-O
  skat_p <- tryCatch(SKAT(gmat, null_model)$p.value, error=function(e) NA)
  skato_p <- tryCatch(SKAT(gmat, null_model, method="SKATO")$p.value, error=function(e) NA)

  # CAST and CMC
  x_bin <- ifelse(gmat > 0, 1, 0)
  cast_p <- NA; cast_beta <- NA
  cmc_p <- NA; cmc_beta <- NA

  try({
    glm_res <- glm(y ~ x_bin, family="binomial")
    cast_p <- summary(glm_res)$coefficients[2, 4]
    cast_beta <- summary(glm_res)$coefficients[2, 1]
  }, silent=TRUE)

  try({
    cmc_bin <- as.factor(x_bin)
    full_data <- cbind(pheno, cmc_bin = cmc_bin)
    glm_cmc <- glm(Phenotype ~ cmc_bin, data = full_data, family = "binomial")
    cmc_p <- summary(glm_cmc)$coefficients[2, 4]
    cmc_beta <- summary(glm_cmc)$coefficients[2, 1]
  }, silent=TRUE)

  # ACAT meta test
  acat_p <- tryCatch({
    pvals <- c(skat_p, cast_p, cmc_p)
    weights <- rep(1, length(pvals))
    T_acat <- sum(weights * tan((0.5 - pvals) * pi))
    0.5 - atan(T_acat / sum(weights)) / pi
  }, error = function(e) NA)

  results <- rbind(results, data.frame(
    Gene = gene,
    NumVariants = num_vars,
    SKAT_P = skat_p,
    SKATO_P = skato_p,
    CAST_P = cast_p,
    CAST_Beta = cast_beta,
    CMC_P = cmc_p,
    CMC_Beta = cmc_beta,
    ACAT_P = acat_p
  ))
}

# ---------------------------- P-VALUE ADJUSTMENTS ---------------------------- #
# Add FDR and Bonferroni correction columns
results$FDR_SKAT <- p.adjust(results$SKAT_P, method="fdr")
results$FDR_SKATO <- p.adjust(results$SKATO_P, method="fdr")
results$FDR_CAST <- p.adjust(results$CAST_P, method="fdr")
results$FDR_CMC <- p.adjust(results$CMC_P, method="fdr")
results$FDR_ACAT <- p.adjust(results$ACAT_P, method="fdr")

results$Bonferroni_SKAT <- p.adjust(results$SKAT_P, method="bonferroni")
results$Bonferroni_SKATO <- p.adjust(results$SKATO_P, method="bonferroni")
results$Bonferroni_CAST <- p.adjust(results$CAST_P, method="bonferroni")
results$Bonferroni_CMC <- p.adjust(results$CMC_P, method="bonferroni")
results$Bonferroni_ACAT <- p.adjust(results$ACAT_P, method="bonferroni")

# ---------------------------- SUMMARY STATS ---------------------------- #
# Generate table with number of significant hits per method and basic metrics
sum_stats <- results %>% summarise(
  SKAT_hits = sum(FDR_SKAT < 0.05, na.rm=TRUE),
  SKATO_hits = sum(FDR_SKATO < 0.05, na.rm=TRUE),
  CAST_hits = sum(FDR_CAST < 0.05, na.rm=TRUE),
  CMC_hits  = sum(FDR_CMC < 0.05, na.rm=TRUE),
  ACAT_hits = sum(FDR_ACAT < 0.05, na.rm=TRUE),
  Avg_Variants = mean(NumVariants, na.rm=TRUE),
  Max_Variants = max(NumVariants, na.rm=TRUE)
)

# ---------------------------- OUTPUT ---------------------------- #
# Save results table and summary statistics
cat("Saving results to:", opt$out, "\n")
write.table(results, opt$out, sep="\t", quote=FALSE, row.names=FALSE)

cat("Saving summary statistics to:", opt$summary, "\n")
write.table(sum_stats, opt$summary, sep="\t", quote=FALSE, row.names=FALSE)

# ---------------------------- DIAGNOSTIC PLOTS ---------------------------- #
# Generate QQ plot and Volcano plot for SKAT and CAST results
if (!is.null(opt$plots)) {
  png(file.path(opt$plots, "qq_plot_skat.png"))
  qq(results$SKAT_P, main="QQ Plot of SKAT p-values")
  dev.off()

  volcano_data <- results %>% mutate(logP = -log10(SKAT_P), Effect = CAST_Beta)
  v <- ggplot(volcano_data, aes(x=Effect, y=logP)) +
    geom_point(alpha=0.7) + theme_minimal() +
    xlab("Effect Size (CAST Beta)") + ylab("-log10(P-value)") +
    ggtitle("Volcano Plot (SKAT vs. CAST)")
  ggsave(filename=file.path(opt$plots, "volcano_plot.png"), plot=v)
}

cat("Done.\n")
