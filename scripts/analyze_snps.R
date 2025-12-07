#!/usr/bin/env Rscript

# --------------------------------------------------------------
# analyze_snps.R
#
# This script reads the VCF produced by snp-sites:
#   03_snp/TS_like_snps.vcf
# and summarizes SNP calls across Atropa TS isoforms and the
# Datura TS candidate CDS.
#
# Output (in 05_results/):
#   TS_like_snp_counts.tsv
#   TS_like_snp_counts.png
# --------------------------------------------------------------

library(tidyverse)

# Create 05_results directory if it does not exist
if (!dir.exists("05_results")) {
  dir.create("05_results")
}

# 1. Read VCF (skip '##' meta lines)
vcf_raw <- read_tsv(
  "03_snp/TS_like_snps.vcf",
  comment = "##",
  col_types = cols(.default = "c")
)

# Identify sample columns
fixed_cols <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
sample_names <- setdiff(names(vcf_raw), fixed_cols)

# 2. If there are no SNPs → produce zero table
if (nrow(vcf_raw) == 0) {
  message("No SNPs detected by snp-sites for TS-like CDS.")
  snp_summary <- tibble(
    sample = sample_names,
    non_missing_calls = 0
  )
  
} else {
  
  # Long format: one row per (site, sample)
  vcf_long <- vcf_raw %>%
    pivot_longer(
      cols = all_of(sample_names),
      names_to = "sample",
      values_to = "genotype"
    )
  
  # 3. Count non-missing genotype calls per sample
  snp_summary <- vcf_long %>%
    dplyr::filter(genotype != "." & !is.na(genotype)) %>%
    dplyr::count(sample, name = "non_missing_calls")
}

# 4. Save summary table
write_tsv(snp_summary, "05_results/TS_like_snp_counts.tsv")

# 5. Plot SNP counts
p_snps <- snp_summary %>%
  ggplot(aes(x = sample, y = non_missing_calls)) +
  geom_col() +
  labs(
    x = "Sequence (Atropa TS isoforms and Datura TS candidate)",
    y = "Number of non-missing SNP calls",
    title = "SNP counts in TS and TS-like CDS"
  ) +
  theme_minimal()

ggsave("05_results/TS_like_snp_counts.png", p_snps, width = 7, height = 5)
