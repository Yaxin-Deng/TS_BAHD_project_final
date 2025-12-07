# scripts/plot_bahd_tree.R
# This script reads a multiple-sequence alignment in FASTA format,
# converts it into an alignment object, computes a distance matrix,
# constructs a Neighbor-Joining (NJ) tree, and saves the tree
# as a PDF and a Newick file in 05_results/.

library(seqinr)
library(ape)

# 1. Path to the aligned FASTA file
align_file <- "03_alignments/BAHD_clade3_ref_plus_DsTS_aligned.fa"

# 2. Basic validation: check that the file exists and is not empty
if (!file.exists(align_file)) {
  stop("Alignment file not found: ", align_file)
}

if (file.info(align_file)$size == 0) {
  stop("Alignment file is empty: ", align_file)
}

# 3. Read FASTA file
# Use seqinr::read.fasta to avoid conflicts with other packages
aln <- seqinr::read.fasta(
  file = align_file,
  seqtype = "AA",
  as.string = TRUE,
  set.attributes = FALSE
)

# Convert list of sequences into an alignment object
ali <- seqinr::as.alignment(
  nb  = length(aln),
  nam = names(aln),
  seq = unlist(aln)
)

# 4. Compute distance matrix using identity model
d <- seqinr::dist.alignment(ali, "identity")

# 5. Construct Neighbor-Joining tree
tree <- ape::nj(d)

# 6. Plot the tree to a PDF (not on RStudio graphics device)
out_pdf <- "05_results/BAHD_tree_nj.pdf"
pdf(out_pdf)
ape::plot.phylo(tree, cex = 0.5)
title("BAHD reference proteins - Neighbor-Joining tree")
dev.off()

# 7. Save the tree in Newick format
out_tree <- "05_results/BAHD_tree_nj.nwk"
ape::write.tree(tree, file = out_tree)

cat("Tree figure saved to: ", out_pdf, "\n")
cat("Newick tree saved to: ", out_tree, "\n")
