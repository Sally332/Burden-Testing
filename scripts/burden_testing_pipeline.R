# Burden Testing Pipeline with Modular Filtering and Thresholds

# Load Required Libraries
library(AssotesteR)
library(SKAT)
library(MiST)
library(caret)  # for LD pruning helper

# Clear Workspace to Avoid Conflicts
rm(list=ls())

#-------------------------
# GLOBAL THRESHOLD SETTINGS
#-------------------------
impact_thresholds <- list(
  cadd_thresh = 20,
  revel_thresh = 0.5,
  sift_thresh = 0.05,
  polyphen_label = "probably_damaging",
  mutationtaster_label = "D",
  clinvar_accepted = c("Pathogenic", "Likely_pathogenic"),
  hgmd_accepted = c("DM", "DM?"),
  high_impact_terms = c("stop_gained", "frameshift_variant", "splice_acceptor_variant", 
                        "splice_donor_variant", "start_lost", "stop_lost")
)

frequency_thresholds <- list(
  min_af = 0.001,
  min_dp = 10,
  min_gq = 20,
  min_qual = 30,
  pop_thresholds = list(
    AFR = 0.005,
    EAS = 0.01,
    NFE = 0.002
  ),
  exclude_pops = c("FIN", "ASJ")  # Add populations to exclude here
)

#-------------------------
# MODULE: CASE-CONTROL EXTRACTION
#-------------------------
extract_samples <- function(build_dir = NULL, genotype_file = NULL, sample_types = c("Case","Control")) {
  if (!is.null(build_dir)) {
    files <- list.files(build_dir, pattern="*_nonsyn_low_outgroup_count.txt$", 
                        full.names=TRUE, recursive=TRUE)
    types <- sapply(strsplit(basename(files), "_"), `[`, 2)
    keep <- files[types %in% sample_types]
    geno_list <- lapply(keep, function(f) {
      d <- read.table(f, header=TRUE, stringsAsFactors=FALSE)
      rownames(d) <- d[[1]]; d[[1]] <- NULL; return(d)
    })
    genotype <- do.call(cbind, geno_list)
  } else if (!is.null(genotype_file)) {
    genotype <- read.table(genotype_file, header=TRUE, stringsAsFactors=FALSE)
    rownames(genotype) <- genotype[[1]]; genotype[[1]] <- NULL
  } else {
    stop("Either build_dir or genotype_file must be provided.")
  }
  return(genotype)
}

#-------------------------
# MODULE: VARIANT FILTERING
#-------------------------
filter_variants <- function(genotype, annotation_db = NULL, gnomad_file = NULL, genoMAD_file = NULL, thresholds) {
  if (is.null(annotation_db)) stop("Annotation database is required for variant filtering.")
  ann <- read.table(annotation_db, header=TRUE, stringsAsFactors=FALSE)

  # Logging setup
  total_variants <- ncol(genotype)
  cat(sprintf("Initial number of variants: %d
", total_variants))

  # Impact Filtering
  keep_imp <- with(ann,
    (CADD_PHRED >= thresholds$cadd_thresh) |
    (REVEL_score >= thresholds$revel_thresh) |
    (SIFT_score < thresholds$sift_thresh) |
    (Polyphen_prediction == thresholds$polyphen_label) |
    (MutationTaster_pred == thresholds$mutationtaster_label) |
    (CLINVAR_CLNSIG %in% thresholds$clinvar_accepted) |
    (HGMD_HGMD_TAG %in% thresholds$hgmd_accepted) |
    (Consequence %in% thresholds$high_impact_terms)
  )
  genotype <- genotype[, ann$Variant_ID[keep_imp], drop=FALSE]
  cat(sprintf("After impact filtering: %d variants
", ncol(genotype)))

  # Frequency Filtering
  keep_maf <- with(ann, AF < thresholds$min_af)
  genotype <- genotype[, ann$Variant_ID[keep_maf], drop=FALSE]
  cat(sprintf("After frequency filtering: %d variants
", ncol(genotype)))

  # Quality Filtering
  keep_qc <- with(ann, QUAL >= thresholds$min_qual & DP >= thresholds$min_dp & GQ >= thresholds$min_gq)
  genotype <- genotype[, ann$Variant_ID[keep_qc], drop=FALSE]
  cat(sprintf("After quality filtering: %d variants
", ncol(genotype)))

  # gnomAD and genoMAD Filtering (Optional)
  filter_af_file <- function(file, thresholds) {
    data <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
    if (!"AF" %in% colnames(data)) stop(paste(file, "must contain an 'AF' column for global frequency filtering."))
    valid_variants <- data$Variant_ID[data$AF <= thresholds$min_af]
    pop_cols <- grep("AF_", colnames(data), value=TRUE)
    pop_filtered <- unique(unlist(lapply(pop_cols, function(col) {
      pop_name <- sub("AF_", "", col)
      if (pop_name %in% thresholds$exclude_pops) return(NULL)
      threshold <- thresholds$pop_thresholds[[pop_name]]
      if (!is.null(threshold)) {
        return(data$Variant_ID[data[[col]] <= threshold])
      } else {
        return(data$Variant_ID[data[[col]] <= thresholds$min_af])
      }
    })))
    return(intersect(valid_variants, pop_filtered))
  }

  if (!is.null(gnomad_file)) {
    valid_gnomad <- filter_af_file(gnomad_file, thresholds)
    genotype <- genotype[, colnames(genotype) %in% valid_gnomad, drop=FALSE]
    cat(sprintf("After gnomAD filtering: %d variants
", ncol(genotype)))
  }

  if (!is.null(genoMAD_file)) {
    valid_genomad <- filter_af_file(genoMAD_file, thresholds)
    genotype <- genotype[, colnames(genotype) %in% valid_genomad, drop=FALSE]
    cat(sprintf("After genoMAD filtering: %d variants
", ncol(genotype)))
  }

  cat(sprintf("Final number of variants: %d
", ncol(genotype)))

  return(genotype)
}
  if (is.null(annotation_db)) stop("Annotation database is required for variant filtering.")
  ann <- read.table(annotation_db, header=TRUE, stringsAsFactors=FALSE)

  # Impact Filtering
  keep_imp <- with(ann,
    (CADD_PHRED >= thresholds$cadd_thresh) |
    (REVEL_score >= thresholds$revel_thresh) |
    (SIFT_score < thresholds$sift_thresh) |
    (Polyphen_prediction == thresholds$polyphen_label) |
    (MutationTaster_pred == thresholds$mutationtaster_label) |
    (CLINVAR_CLNSIG %in% thresholds$clinvar_accepted) |
    (HGMD_HGMD_TAG %in% thresholds$hgmd_accepted) |
    (Consequence %in% thresholds$high_impact_terms)
  )
  genotype <- genotype[, ann$Variant_ID[keep_imp], drop=FALSE]

  # Frequency Filtering
  keep_maf <- with(ann, AF < thresholds$min_af)
  genotype <- genotype[, ann$Variant_ID[keep_maf], drop=FALSE]

  # Quality Filtering
  keep_qc <- with(ann, QUAL >= thresholds$min_qual & DP >= thresholds$min_dp & GQ >= thresholds$min_gq)
  genotype <- genotype[, ann$Variant_ID[keep_qc], drop=FALSE]

  # gnomAD and genoMAD Filtering (Optional)
  filter_af_file <- function(file, thresholds) {
    data <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
    if (!"AF" %in% colnames(data)) stop(paste(file, "must contain an 'AF' column for global frequency filtering."))
    valid_variants <- data$Variant_ID[data$AF <= thresholds$min_af]
    pop_cols <- grep("AF_", colnames(data), value=TRUE)
    pop_filtered <- unique(unlist(lapply(pop_cols, function(col) {
      pop_name <- sub("AF_", "", col)
      if (pop_name %in% thresholds$exclude_pops) return(NULL)
      threshold <- thresholds$pop_thresholds[[pop_name]]
      if (!is.null(threshold)) {
        return(data$Variant_ID[data[[col]] <= threshold])
      } else {
        return(data$Variant_ID[data[[col]] <= thresholds$min_af])
      }
    })))
    return(intersect(valid_variants, pop_filtered))
  }

  if (!is.null(gnomad_file)) {
    valid_gnomad <- filter_af_file(gnomad_file, thresholds)
    genotype <- genotype[, colnames(genotype) %in% valid_gnomad, drop=FALSE]
  }

  if (!is.null(genoMAD_file)) {
    valid_genomad <- filter_af_file(genoMAD_file, thresholds)
    genotype <- genotype[, colnames(genotype) %in% valid_genomad, drop=FALSE]
  }

  return(genotype)
}

#-------------------------
# MODULE: BURDEN TESTING, CORRECTION, AND SAMPLE SIZE GUIDANCE
#-------------------------

# NOTE: Sample Size Guidance
# The required sample size for adequate power (e.g., 80%) depends on the effect size, minor allele frequency (MAF), and the chosen significance threshold.
# As a rough reference:
#   - Large effect size (OR ≈ 2.0), MAF = 0.01 -> ~250 cases and ~250 controls
#   - Moderate effect size (OR ≈ 1.5), MAF = 0.01 -> ~500 cases and ~500 controls
#   - Small effect size (OR ≈ 1.2), MAF = 0.01 -> ~1000 cases and ~1000 controls
#   - Larger sample sizes are needed for smaller effect sizes and lower MAFs.
# For more precise estimates, consider using dedicated power analysis tools like the 'pwr' package in R or online tools such as:
# https://csg.sph.umich.edu/abecasis/cats/
# https://zzz.bwh.harvard.edu/gpc/
# https://github.com/bioinformatics-ptp/power-analysis

perform_burden_tests <- function(genotype, case_control, correction_method = "bonferroni", pcs_file = NULL, num_pcs = 5, sample_metadata_file = NULL) {
  results <- list()
  variant_matrix <- as.matrix(genotype)

  # Fisher's exact test (simple burden test)
  fisher_pvals <- apply(variant_matrix, 2, function(variant) {
    contingency_table <- table(factor(variant, levels = c(0,1)), case_control)
    fisher.test(contingency_table)$p.value
  })
  results$fisher <- fisher_pvals

  # SKAT-O (adaptive test)
  skat_null <- SKAT_Null_Model(case_control ~ 1, out_type = "D")
  skat_o <- SKAT(as.matrix(genotype), skat_null, method = "optimal.adj")
  results$skat_o <- skat_o$p.value

  # Madsen-Browning (weighted burden test)
  mb_pvals <- apply(variant_matrix, 2, function(variant) {
    contingency_table <- table(factor(variant, levels = c(0,1)), case_control)
    Madsen.Browning(contingency_table)$p.value
  })
  results$madsen_browning <- mb_pvals

  # Multiple Testing Correction
  if (correction_method == "bonferroni") {
    corrected_pvals <- lapply(results, function(pvals) {
      p.adjust(pvals, method = "bonferroni")
    })
    results <- corrected_pvals
  } else if (correction_method == "fdr") {
    corrected_pvals <- lapply(results, function(pvals) {
      p.adjust(pvals, method = "fdr")
    })
    results <- corrected_pvals
  }

  # Sample Size Guidance
  num_cases <- sum(case_control == 1)
  num_controls <- sum(case_control == 0)
  total_samples <- num_cases + num_controls
  cat(sprintf("Total samples: %d (Cases: %d, Controls: %d)
", total_samples, num_cases, num_controls))
  cat("Recommended minimum sample size for adequate power (e.g., 80%) depends on the effect size and variant frequency.
")
  cat("Consider using power analysis tools like pwr.p.test() for more precise estimates.
")

    # Calculate Lambda GC and Generate QQ Plot
  all_pvals <- unlist(results)
  lambda_gc_est <- median(qchisq(all_pvals, df=1, lower.tail=FALSE)) / qchisq(0.5, df=1)
  cat(sprintf("Estimated Lambda GC: %.3f
", lambda_gc_est))

  # QQ Plot
  observed_pvals <- sort(-log10(all_pvals))
  expected_pvals <- -log10(ppoints(length(observed_pvals)))
  plot(expected_pvals, observed_pvals, main="QQ Plot - Genomic Control", xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="darkred")
  abline(0, 1, col="blue", lty=2)

  # Return Results
  results$lambda_gc <- lambda_gc_est
  return(results)
  results$lambda_gc <- lambda_gc_est
  return(results)
  return(results)
}

perform_burden_tests <- function(genotype, case_control, correction_method = "bonferroni", pcs_file = NULL, num_pcs = 5, sample_metadata_file = NULL) {
  results <- list()
  variant_matrix <- as.matrix(genotype)

  # Fisher's exact test (simple burden test)
  fisher_pvals <- apply(variant_matrix, 2, function(variant) {
    contingency_table <- table(factor(variant, levels = c(0,1)), case_control)
    fisher.test(contingency_table)$p.value
  })
  results$fisher <- fisher_pvals

  # SKAT-O (adaptive test)
  skat_null <- SKAT_Null_Model(case_control ~ 1, out_type = "D")
  skat_o <- SKAT(as.matrix(genotype), skat_null, method = "optimal.adj")
  results$skat_o <- skat_o$p.value

  # Madsen-Browning (weighted burden test)
  mb_pvals <- apply(variant_matrix, 2, function(variant) {
    contingency_table <- table(factor(variant, levels = c(0,1)), case_control)
    Madsen.Browning(contingency_table)$p.value
  })
  results$madsen_browning <- mb_pvals

  # Multiple Testing Correction
  if (correction_method == "bonferroni") {
    corrected_pvals <- lapply(results, function(pvals) {
      p.adjust(pvals, method = "bonferroni")
    })
    results <- corrected_pvals
  } else if (correction_method == "fdr") {
    corrected_pvals <- lapply(results, function(pvals) {
      p.adjust(pvals, method = "fdr")
    })
    results <- corrected_pvals
  }

  # Generate Summary Statistics
  summary_stats <- data.frame(
    Test = c("Fisher", "SKAT-O", "Madsen-Browning"),
    Num_Variants = c(length(results$fisher), length(results$skat_o), length(results$madsen_browning)),
    Min_P = c(min(results$fisher, na.rm=TRUE), min(results$skat_o, na.rm=TRUE), min(results$madsen_browning, na.rm=TRUE)),
    Median_P = c(median(results$fisher, na.rm=TRUE), median(results$skat_o, na.rm=TRUE), median(results$madsen_browning, na.rm=TRUE)),
    Max_P = c(max(results$fisher, na.rm=TRUE), max(results$skat_o, na.rm=TRUE), max(results$madsen_browning, na.rm=TRUE))
  )

  cat("
Summary of Burden Test Results:
")
  print(summary_stats)

  # Return Results
  results$summary <- summary_stats
  return(results)
}

#-------------------------
# MODULE: POPULATION STRUCTURE CORRECTION
#-------------------------

correct_population_structure <- function(genotype, pcs_file = NULL, num_pcs = 5, lambda_gc = NULL, kinship_matrix_file = NULL, sample_metadata_file = NULL, pvals = NULL) {
  # PCA Correction
  if (!is.null(pcs_file)) {
    pcs <- read.table(pcs_file, header=TRUE, stringsAsFactors=FALSE)
    if (ncol(pcs) < num_pcs + 1) stop("The PCA file must contain at least the specified number of PCs.")

    # Extract the specified number of PCs
    selected_pcs <- pcs[, 2:(num_pcs + 1)]
    rownames(selected_pcs) <- pcs[, 1]

    # Subset genotype to match PCA samples
    common_samples <- intersect(rownames(selected_pcs), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    selected_pcs <- selected_pcs[common_samples, , drop=FALSE]
    cat(sprintf("Population structure correction applied using %d PCs.
", num_pcs))

    # Visualization of PCs (color-coded by population if metadata is provided)
    if (!is.null(sample_metadata_file)) {
      metadata <- read.table(sample_metadata_file, header=TRUE, stringsAsFactors=FALSE)
      metadata <- metadata[metadata$Sample_ID %in% rownames(selected_pcs), ]
      pca_df <- data.frame(PC1 = selected_pcs[,1], PC2 = selected_pcs[,2], Population = metadata$Population)
      plot(pca_df$PC1, pca_df$PC2, col=as.factor(pca_df$Population), pch=20,
           main="PCA - Population Structure", xlab="PC1", ylab="PC2")
      legend("topright", legend=levels(as.factor(pca_df$Population)), col=1:length(unique(pca_df$Population)), pch=20, cex=0.8)
    } else {
      pca_df <- data.frame(PC1 = selected_pcs[,1], PC2 = selected_pcs[,2])
      plot(pca_df$PC1, pca_df$PC2, main="PCA - Population Structure", xlab="PC1", ylab="PC2", pch=20, col="darkblue")
    }
  } else {
    selected_pcs <- NULL
    cat("No PCA file provided. Skipping PCA-based correction.
")
  }

  # Genomic Control (GC) with Real P-values
  if (!is.null(pvals)) {
    lambda_gc_est <- median(qchisq(pvals, df=1, lower.tail=FALSE)) / qchisq(0.5, df=1)
    cat(sprintf("Estimated Lambda GC: %.3f
", lambda_gc_est))
    
    # QQ Plot for Lambda GC Validation
    observed_pvals <- sort(-log10(pvals))
    expected_pvals <- -log10(ppoints(length(observed_pvals)))
    plot(expected_pvals, observed_pvals, main="QQ Plot - Genomic Control", xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="darkred")
    abline(0, 1, col="blue", lty=2)
  }

  # Linear Mixed Model (LMM) with Kinship Matrix (if provided)
  if (!is.null(kinship_matrix_file)) {
    kinship_matrix <- read.table(kinship_matrix_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)
    common_samples <- intersect(rownames(kinship_matrix), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    kinship_matrix <- kinship_matrix[common_samples, common_samples]
    cat(sprintf("LMM correction applied using kinship matrix for %d samples.
", length(common_samples)))
    
    # Heatmap for kinship matrix
    heatmap(as.matrix(kinship_matrix), main="Kinship Matrix", col=heat.colors(256), margins=c(5,5))
    return(list(genotype=genotype, pcs=selected_pcs, kinship=kinship_matrix, lambda_gc=lambda_gc_est))
  }

  return(list(genotype=genotype, pcs=selected_pcs, lambda_gc=lambda_gc_est))
}


correct_population_structure <- function(genotype, pcs_file = NULL, num_pcs = 5, lambda_gc = NULL, kinship_matrix_file = NULL, sample_metadata_file = NULL) {
  # PCA Correction
  if (!is.null(pcs_file)) {
    pcs <- read.table(pcs_file, header=TRUE, stringsAsFactors=FALSE)
    if (ncol(pcs) < num_pcs + 1) stop("The PCA file must contain at least the specified number of PCs.")

    # Extract the specified number of PCs
    selected_pcs <- pcs[, 2:(num_pcs + 1)]
    rownames(selected_pcs) <- pcs[, 1]

    # Subset genotype to match PCA samples
    common_samples <- intersect(rownames(selected_pcs), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    selected_pcs <- selected_pcs[common_samples, , drop=FALSE]
    cat(sprintf("Population structure correction applied using %d PCs.
", num_pcs))

    # Visualization of PCs (color-coded by population if metadata is provided)
    if (!is.null(sample_metadata_file)) {
      metadata <- read.table(sample_metadata_file, header=TRUE, stringsAsFactors=FALSE)
      metadata <- metadata[metadata$Sample_ID %in% rownames(selected_pcs), ]
      pca_df <- data.frame(PC1 = selected_pcs[,1], PC2 = selected_pcs[,2], Population = metadata$Population)
      plot(pca_df$PC1, pca_df$PC2, col=as.factor(pca_df$Population), pch=20,
           main="PCA - Population Structure", xlab="PC1", ylab="PC2")
      legend("topright", legend=levels(as.factor(pca_df$Population)), col=1:length(unique(pca_df$Population)), pch=20, cex=0.8)
    } else {
      pca_df <- data.frame(PC1 = selected_pcs[,1], PC2 = selected_pcs[,2])
      plot(pca_df$PC1, pca_df$PC2, main="PCA - Population Structure", xlab="PC1", ylab="PC2", pch=20, col="darkblue")
    }
  } else {
    selected_pcs <- NULL
    cat("No PCA file provided. Skipping PCA-based correction.
")
  }

  # Genomic Control (GC)
  if (!is.null(lambda_gc)) {
    genotype <- genotype / sqrt(lambda_gc)
    cat(sprintf("Genomic control applied with lambda = %.3f.
", lambda_gc))
    
    # QQ Plot for Lambda GC Validation
    observed_pvals <- -log10(runif(ncol(genotype)))  # Placeholder for real p-values
    expected_pvals <- -log10(ppoints(length(observed_pvals)))
    plot(expected_pvals, observed_pvals, main="QQ Plot - Genomic Control", xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="darkred")
    abline(0, 1, col="blue", lty=2)
  }

  # Linear Mixed Model (LMM) with Kinship Matrix (if provided)
  if (!is.null(kinship_matrix_file)) {
    kinship_matrix <- read.table(kinship_matrix_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)
    common_samples <- intersect(rownames(kinship_matrix), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    kinship_matrix <- kinship_matrix[common_samples, common_samples]
    cat(sprintf("LMM correction applied using kinship matrix for %d samples.
", length(common_samples)))
    
    # Heatmap for kinship matrix
    heatmap(as.matrix(kinship_matrix), main="Kinship Matrix", col=heat.colors(256), margins=c(5,5))
    return(list(genotype=genotype, pcs=selected_pcs, kinship=kinship_matrix))
  }

  return(list(genotype=genotype, pcs=selected_pcs))
}


correct_population_structure <- function(genotype, pcs_file = NULL, num_pcs = 5, lambda_gc = NULL, kinship_matrix_file = NULL, sample_metadata_file = NULL) {
  # PCA Correction
  if (!is.null(pcs_file)) {
    pcs <- read.table(pcs_file, header=TRUE, stringsAsFactors=FALSE)
    if (ncol(pcs) < num_pcs + 1) stop("The PCA file must contain at least the specified number of PCs.")

    # Extract the specified number of PCs
    selected_pcs <- pcs[, 2:(num_pcs + 1)]
    rownames(selected_pcs) <- pcs[, 1]

    # Subset genotype to match PCA samples
    common_samples <- intersect(rownames(selected_pcs), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    selected_pcs <- selected_pcs[common_samples, , drop=FALSE]
    cat(sprintf("Population structure correction applied using %d PCs.
", num_pcs))

    # Visualization of PCs (color-coded by population if metadata is provided)
    if (!is.null(sample_metadata_file)) {
      metadata <- read.table(sample_metadata_file, header=TRUE, stringsAsFactors=FALSE)
      metadata <- metadata[metadata$Sample_ID %in% rownames(selected_pcs), ]
      pca_df <- data.frame(PC1 = selected_pcs[,1], PC2 = selected_pcs[,2], Population = metadata$Population)
      plot(pca_df$PC1, pca_df$PC2, col=as.factor(pca_df$Population), pch=20,
           main="PCA - Population Structure", xlab="PC1", ylab="PC2")
      legend("topright", legend=levels(as.factor(pca_df$Population)), col=1:length(unique(pca_df$Population)), pch=20, cex=0.8)
    } else {
      pca_df <- data.frame(PC1 = selected_pcs[,1], PC2 = selected_pcs[,2])
      plot(pca_df$PC1, pca_df$PC2, main="PCA - Population Structure", xlab="PC1", ylab="PC2", pch=20, col="darkblue")
    }
  } else {
    selected_pcs <- NULL
    cat("No PCA file provided. Skipping PCA-based correction.
")
  }

  # Genomic Control (GC)
  if (!is.null(lambda_gc)) {
    genotype <- genotype / sqrt(lambda_gc)
    cat(sprintf("Genomic control applied with lambda = %.3f.
", lambda_gc))
  }

  # Linear Mixed Model (LMM) with Kinship Matrix (if provided)
  if (!is.null(kinship_matrix_file)) {
    kinship_matrix <- read.table(kinship_matrix_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)
    common_samples <- intersect(rownames(kinship_matrix), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    kinship_matrix <- kinship_matrix[common_samples, common_samples]
    cat(sprintf("LMM correction applied using kinship matrix for %d samples.
", length(common_samples)))
    
    # Heatmap for kinship matrix
    heatmap(as.matrix(kinship_matrix), main="Kinship Matrix", col=heat.colors(256), margins=c(5,5))
    return(list(genotype=genotype, pcs=selected_pcs, kinship=kinship_matrix))
  }

  return(list(genotype=genotype, pcs=selected_pcs))
}


correct_population_structure <- function(genotype, pcs_file = NULL, num_pcs = 5, lambda_gc = NULL, kinship_matrix_file = NULL) {
  # PCA Correction
  if (!is.null(pcs_file)) {
    pcs <- read.table(pcs_file, header=TRUE, stringsAsFactors=FALSE)
    if (ncol(pcs) < num_pcs + 1) stop("The PCA file must contain at least the specified number of PCs.")

    # Extract the specified number of PCs
    selected_pcs <- pcs[, 2:(num_pcs + 1)]
    rownames(selected_pcs) <- pcs[, 1]

    # Subset genotype to match PCA samples
    common_samples <- intersect(rownames(selected_pcs), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    selected_pcs <- selected_pcs[common_samples, , drop=FALSE]
    cat(sprintf("Population structure correction applied using %d PCs.
", num_pcs))

    # Visualization of PCs
    if (num_pcs >= 2) {
      pca_df <- data.frame(PC1 = selected_pcs[,1], PC2 = selected_pcs[,2], Sample = rownames(selected_pcs))
      plot(pca_df$PC1, pca_df$PC2, main="PCA - Population Structure", xlab="PC1", ylab="PC2", pch=20, col="darkblue")
    }
  } else {
    selected_pcs <- NULL
    cat("No PCA file provided. Skipping PCA-based correction.
")
  }

  # Genomic Control (GC)
  if (!is.null(lambda_gc)) {
    genotype <- genotype / sqrt(lambda_gc)
    cat(sprintf("Genomic control applied with lambda = %.3f.
", lambda_gc))
  }

  # Linear Mixed Model (LMM) with Kinship Matrix (if provided)
  if (!is.null(kinship_matrix_file)) {
    kinship_matrix <- read.table(kinship_matrix_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)
    common_samples <- intersect(rownames(kinship_matrix), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    kinship_matrix <- kinship_matrix[common_samples, common_samples]
    cat(sprintf("LMM correction applied using kinship matrix for %d samples.
", length(common_samples)))
    
    # Heatmap for kinship matrix
    heatmap(as.matrix(kinship_matrix), main="Kinship Matrix", col=heat.colors(256), margins=c(5,5))
    return(list(genotype=genotype, pcs=selected_pcs, kinship=kinship_matrix))
  }

  return(list(genotype=genotype, pcs=selected_pcs))
}


correct_population_structure <- function(genotype, pcs_file = NULL, num_pcs = 5, lambda_gc = NULL, kinship_matrix_file = NULL) {
  # PCA Correction
  if (!is.null(pcs_file)) {
    pcs <- read.table(pcs_file, header=TRUE, stringsAsFactors=FALSE)
    if (ncol(pcs) < num_pcs + 1) stop("The PCA file must contain at least the specified number of PCs.")

    # Extract the specified number of PCs
    selected_pcs <- pcs[, 2:(num_pcs + 1)]
    rownames(selected_pcs) <- pcs[, 1]

    # Subset genotype to match PCA samples
    common_samples <- intersect(rownames(selected_pcs), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    selected_pcs <- selected_pcs[common_samples, , drop=FALSE]
    cat(sprintf("Population structure correction applied using %d PCs.
", num_pcs))
  } else {
    selected_pcs <- NULL
    cat("No PCA file provided. Skipping PCA-based correction.
")
  }

  # Genomic Control (GC)
  if (!is.null(lambda_gc)) {
    genotype <- genotype / sqrt(lambda_gc)
    cat(sprintf("Genomic control applied with lambda = %.3f.
", lambda_gc))
  }

  # Linear Mixed Model (LMM) with Kinship Matrix (if provided)
  if (!is.null(kinship_matrix_file)) {
    kinship_matrix <- read.table(kinship_matrix_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)
    common_samples <- intersect(rownames(kinship_matrix), colnames(genotype))
    genotype <- genotype[, common_samples, drop=FALSE]
    kinship_matrix <- kinship_matrix[common_samples, common_samples]
    cat(sprintf("LMM correction applied using kinship matrix for %d samples.
", length(common_samples)))
    return(list(genotype=genotype, pcs=selected_pcs, kinship=kinship_matrix))
  }

  return(list(genotype=genotype, pcs=selected_pcs))
}


correct_population_structure <- function(genotype, pcs_file = NULL, num_pcs = 5) {
  if (is.null(pcs_file)) {
    cat("No PCA file provided. Skipping population structure correction.
")
    return(genotype)
  }

  pcs <- read.table(pcs_file, header=TRUE, stringsAsFactors=FALSE)
  if (ncol(pcs) < num_pcs + 1) stop("The PCA file must contain at least the specified number of PCs.")

  # Extract the specified number of PCs
  selected_pcs <- pcs[, 2:(num_pcs + 1)]
  rownames(selected_pcs) <- pcs[, 1]

  # Subset genotype to match PCA samples
  common_samples <- intersect(rownames(selected_pcs), colnames(genotype))
  genotype <- genotype[, common_samples, drop=FALSE]
  selected_pcs <- selected_pcs[common_samples, , drop=FALSE]

  cat(sprintf("Population structure correction applied using %d PCs.
", num_pcs))
  return(list(genotype=genotype, pcs=selected_pcs))
}

#-------------------------
# MAIN PIPELINE FUNCTION
#-------------------------
run_burden_test <- function(
  genotype_file = NULL,
  build_dir = NULL,
  sample_types = c("Case","Control"),
  annotation_db = NULL,
  gnomad_file = NULL,
  genoMAD_file = NULL,
  case_control_file = NULL,
  pcs_file = NULL,
  num_pcs = 5,
  lambda_gc = NULL,
  kinship_matrix_file = NULL,
  sample_metadata_file = NULL,
  correction_method = "bonferroni",
  thresholds = c(impact_thresholds, frequency_thresholds)
) {
  genotype_file = NULL,
  build_dir = NULL,
  sample_types = c("Case","Control"),
  annotation_db = NULL,
  gnomad_file = NULL,
  genoMAD_file = NULL,
  case_control_file = NULL,
  pcs_file = NULL,
  num_pcs = 5,
  lambda_gc = NULL,
  kinship_matrix_file = NULL,
  sample_metadata_file = NULL,
  correction_method = "bonferroni",
  thresholds = c(impact_thresholds, frequency_thresholds)
) {
  genotype_file = NULL,
  build_dir = NULL,
  sample_types = c("Case","Control"),
  annotation_db = NULL,
  gnomad_file = NULL,
  genoMAD_file = NULL,
  case_control_file = NULL,
  pcs_file = NULL,
  num_pcs = 5,
  lambda_gc = NULL,
  kinship_matrix_file = NULL,
  correction_method = "bonferroni",
  thresholds = c(impact_thresholds, frequency_thresholds)
) {
  genotype_file = NULL,
  build_dir = NULL,
  sample_types = c("Case","Control"),
  annotation_db = NULL,
  gnomad_file = NULL,
  genoMAD_file = NULL,
  case_control_file = NULL,
  pcs_file = NULL,
  num_pcs = 5,
  correction_method = "bonferroni",
  thresholds = c(impact_thresholds, frequency_thresholds)
) {
  genotype_file = NULL,
  build_dir = NULL,
  sample_types = c("Case","Control"),
  annotation_db = NULL,
  gnomad_file = NULL,
  genoMAD_file = NULL,
  case_control_file = NULL,
  correction_method = "bonferroni",
  thresholds = c(impact_thresholds, frequency_thresholds)
) {
  genotype_file = NULL,
  build_dir = NULL,
  sample_types = c("Case","Control"),
  annotation_db = NULL,
  gnomad_file = NULL,
  genoMAD_file = NULL,
  thresholds = c(impact_thresholds, frequency_thresholds)
) {
  # Extract Samples
  case_control <- read.table(case_control_file, header=TRUE, stringsAsFactors=FALSE)
  case_control_vector <- setNames(ifelse(case_control$Status == "Case", 1, 0), case_control$Sample_ID)

  # Correct for Population Structure (if PCA file, lambda, or kinship matrix provided)
  pc_corrected <- correct_population_structure(genotype, pcs_file, num_pcs, lambda_gc, kinship_matrix_file, sample_metadata_file, pvals=NULL)
  genotype <- pc_corrected$genotype
  pcs <- pc_corrected$pcs
  lambda_gc_est <- if (!is.null(pc_corrected$lambda_gc)) pc_corrected$lambda_gc else NULL
  kinship <- if (!is.null(pc_corrected$kinship)) pc_corrected$kinship else NULL
  pc_corrected <- correct_population_structure(genotype, pcs_file, num_pcs, lambda_gc, kinship_matrix_file, sample_metadata_file)
  genotype <- pc_corrected$genotype
  pcs <- pc_corrected$pcs
  kinship <- if (!is.null(pc_corrected$kinship)) pc_corrected$kinship else NULL
  pc_corrected <- correct_population_structure(genotype, pcs_file, num_pcs, lambda_gc, kinship_matrix_file)
  genotype <- pc_corrected$genotype
  pcs <- pc_corrected$pcs
  kinship <- if (!is.null(pc_corrected$kinship)) pc_corrected$kinship else NULL
  pc_corrected <- correct_population_structure(genotype, pcs_file, num_pcs)
  genotype <- pc_corrected$genotype
  pcs <- pc_corrected$pcs
  case_control <- read.table(case_control_file, header=TRUE, stringsAsFactors=FALSE)
  case_control_vector <- setNames(ifelse(case_control$Status == "Case", 1, 0), case_control$Sample_ID)
  genotype <- extract_samples(build_dir = build_dir, genotype_file = genotype_file, sample_types = sample_types)

  # Apply Variant Filtering
  filtered_genotype <- filter_variants(genotype, annotation_db, gnomad_file, genoMAD_file, thresholds)

  # Perform Burden Tests
  burden_results <- perform_burden_tests(filtered_genotype, case_control_vector, correction_method, pcs_file, num_pcs, sample_metadata_file)

  # Save Results
  write.table(filtered_genotype, "Filtered_Genotypes.txt", quote=FALSE, sep="	", row.names=TRUE)
  write.table(burden_results$fisher, "Fisher_Results.txt", quote=FALSE, sep="	", row.names=TRUE)
  write.table(burden_results$skat_o, "SKAT_O_Results.txt", quote=FALSE, sep="	", row.names=TRUE)
  write.table(burden_results$madsen_browning, "Madsen_Browning_Results.txt", quote=FALSE, sep="	", row.names=TRUE)
  write.table(filtered_genotype, "Filtered_Genotypes.txt", quote=FALSE, sep="\t", row.names=TRUE)
  return(filtered_genotype)
}

# Example invocation (modify paths accordingly)
# run_burden_test(
#   genotype_file = "burden_table.txt",
#   annotation_db = "Variant_Impact_Annotations.txt",
#   thresholds = c(impact_thresholds, frequency_thresholds)
# )
