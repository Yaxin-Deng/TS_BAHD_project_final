#!/usr/bin/env Rscript

# --------------------------------------------------------------
# analyze_blast.R
#
# This script reads the BLASTP output (outfmt 6) from:
#   02_blast/blast_clade3_vs_Ds_proteome.tab
# It extracts the top hit for each query (highest bitscore),
# saves a summary table, and creates a barplot of percent
# identity for the top hits.
#
# Input:
#   02_blast/blast_clade3_vs_Ds_proteome.tab
#
# Output (saved in 05_results/):
#   blast_top_hits_summary.tsv
#   blast_top_hits.png
# --------------------------------------------------------------

# Load tidyverse (as in the course)
library(tidyverse)

# Make sure the results directory exists
if (!dir.exists("05_results")) {
  dir.create("05_results")
}

# 1. Read BLAST outfmt 6 table
blast <- read_tsv(
  "02_blast/blast_clade3_vs_Ds_proteome.tab",
  col_names = c(
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore"
  )
)

# 2. For each query, keep the top hit (highest bitscore)
blast_top <- blast |>
  group_by(qseqid) |>
  slice_max(order_by = bitscore, n = 1) |>
  ungroup()

# 3. Save summary table
write_tsv(blast_top, "05_results/blast_top_hits_summary.tsv")

# 4. Create barplot of percent identity for top hits
plot_top_hits <- blast_top |>
  ggplot(aes(x = pident,
             y = fct_reorder(qseqid, pident))) +
  geom_col() +
  labs(
    x = "Percent identity (%)",
    y = "Query (clade 3 BAHD proteins)",
    title = "Top BLAST hits of clade 3 BAHD proteins in the Datura proteome"
  ) +
  theme_minimal()

# 5. Save plot to 05_results
ggsave(
  filename = "05_results/blast_top_hits.png",
  plot = plot_top_hits,
  width = 8,
  height = 6
)
