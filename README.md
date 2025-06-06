
# Burden Testing Pipeline

---

### 📉 Overview

This pipeline supports **gene-based burden testing** of rare variants using SKAT and related methods. It is designed for **flexibility** in variant filtering, population structure correction, and real-case sample setups. It includes options for **covariates**, **QC visualization**, and **multiple filtering combinations** to adapt to clinical or family-based studies. 

All scripts are richly commented to guide users through the logic, assumptions, and interpretation of each step.

---

### 🧬 Variant Filtering: `variant_filtering_main.R`
Filtering is a critical step to ensure that the burden tests focus on variants most likely to affect gene function and contribute to disease. This script enables:

- Support for multiple file formats (e.g., subject-based tables or binary presence/absence matrices)
- Filtering by functional consequence (e.g., missense, nonsense)
- Frequency pruning using sources like gnomAD
- Prioritization by functional scores (e.g., CADD)
- Segregation filtering (e.g., ≥ N affected individuals in a family, absent in controls/outgroup)
- Optional filtering by clinical annotations, inheritance model, or panel membership
- Option to aggregate and test variants by custom genomic regions (provided as BED or TSV format), instead of by gene

---

### 🧪 Burden Testing: `burden_testing_pipeline.R`
This script performs **gene-based rare variant association testing** with several statistical methods:

- SKAT and SKAT-O: variance-component tests
- CMC, CAST, and ACAT: burden-style tests
- Mixed Models: correct for relatedness, ancestry, batch using covariates
- Outputs: results are saved in TSV format and ready for downstream analysis

---

### 📊 Additional Modules and Considerations
- **QC Plots**: PCA plot, QQ plot, variant count histograms, volcano plot
- **Sample Size Guidance**: Power estimation functions available in SKAT

---

### 🧾 Usage Example
**Filtering step:**
```bash
Rscript variant_filtering_main.R \
  --input data/input.tsv \
  --output data/filtered.tsv \
  --format subjects \
  --min_carriers 3 \
  --af_max 0.01 \
  --cadd_min 15
```

**Burden test step:**
```bash
Rscript burden_testing_pipeline.R \
  --geno data/geno.tsv \
  --pheno data/pheno.tsv \
  --covar data/covar.tsv \
  --out results/burden_results.tsv \
  --plots results/plots \
  --summary results/summary_stats.tsv
```

---

### 📚 Full Documentation
For complete explanation of methods and advanced usage:  
→ `docs/README.pdf`

---

### 📬 Contact & License
**Author:** Sally Yepes (sallyepes233@gmail.com)  
**License:** MIT

