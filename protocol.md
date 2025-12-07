# Project protocol: BAHD acyltransferases and TS-like enzymes in *Datura stramonium*

## 1. Project overview

This project focuses on BAHD acyltransferases that are closely related to the *Atropa belladonna* tigloyltransferase (TS; aba_locus_5896) and on TS-like enzymes in *Datura stramonium*. Using protein sequences from the published clade 3 BAHD phylogeny (Fig. 2E) and the *D. stramonium* proteome, I will:

1. Identify TS-like BAHD acyltransferases in *D. stramonium* using BLAST.
2. Place *Datura* candidates in a BAHD clade 3 phylogenetic tree.
3. Compare coding sequences of TS and TS-like acyltransferases using alignment-based SNP analysis.

## 2. Directory structure

- `00_raw/` — raw input files (FASTA sequences).
- `02_blast/` — BLAST databases and tab-delimited BLAST results.
- `03_alignments/` — multiple sequence alignments (protein and CDS).
- `03_snp/` — SNP-related files (aligned CDS, VCF).
- `04_phylogeny/` — tree files (Newick) and intermediate outputs.
- `05_results/` — final figures and summary tables used in the report.
- `scripts/` — all shell and R scripts used in the analysis.
- `protocol.md` — this protocol.
- `notebook.md` — detailed log/diary of commands and notes.

## 3. Software used

All tools are from the course materials and OSC environment:

- BLAST+ (makeblastdb, blastp)
- MAFFT (for multiple sequence alignment)
- snp-sites (for SNP detection)
- R (tidyverse, ape / ggtree for plotting)

## 4. Analysis steps

### 4.1 BLAST analysis: clade 3 BAHD vs *Datura* proteome

**Goal.**  
Find *D. stramonium* proteins that are most similar to known clade 3 BAHD acyltransferases, including the *A. belladonna* TS and TS-like candidates.

**Input files (from `00_raw/`):**

- `BAHD_clade3_ref_plus_DsTS.fa` — clade 3 BAHD proteins from Fig. 2E plus TS-like sequences from *Datura*.
- `Datura_genome_protein.fa` — *D. stramonium* proteome.

**Script:** `scripts/run_blast.sh`

This script:

1. Loads the BLAST+ module.
2. Builds a protein BLAST database from `Datura_genome_protein.fa`.
3. Runs `blastp` with `BAHD_clade3_ref_plus_DsTS.fa` as queries against the database.
4. Saves the tab-delimited results in `02_blast/blast_clade3_vs_Ds_proteome.tab`.

**How to run:**

```bash
cd /fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final

chmod +x scripts/run_blast.sh

sbatch scripts/run_blast.sh
```
The output was:
```
Submitted batch job 42368936
```

The output file was:
```
02_blast/BAHD_clade3_vs_Ds_proteome.tab
```

Outputs:

Intermediate:

02_blast/Datura_protein_db.* — BLAST protein database.

02_blast/blast_clade3_vs_Ds_proteome.tab — BLASTP output in outfmt 6.

Final (after R analysis, see Section 4.4):

A summary table of top BLAST hits and at least one plot saved in 05_results/ (e.g. 05_results/blast_top_hits_barplot.png).


### 4.2 MAFFT + phylogeny
Multiple sequence alignment and phylogeny (protein)

Planned steps:

- Combine clade 3 BAHD reference proteins and Datura TS-like sequences into one FASTA.

- Use MAFFT to align the protein sequences and save the alignment in 03_alignments/.

- Build a phylogenetic tree (using a method covered in the course or derived from the alignment).

- Visualize the tree in R and save annotated tree figures in 05_results/.

**Goal.**  
Align clade 3 BAHD protein sequences (including the *A. belladonna* TS and the
putative *D. stramonium* TS-like proteins) and construct a simple phylogenetic
tree to visualize their relationships.

**Input files:**

- `00_raw/BAHD_clade3_ref_plus_DsTS.fa` — clade 3 BAHD proteins from Fig. 2E plus TS-like sequences.

#### 4.2.1 MAFFT alignment (OSC bash)

**Script:** `scripts/mafft_bahd.slurm`

This Slurm script:

1. Loads the MAFFT module.
2. Aligns all sequences in `00_raw/BAHD_clade3_ref_plus_DsTS.fa` using `mafft --auto`.
3. Saves the alignment to `03_alignments/BAHD_clade3_plus_DsTS_aligned.fa`.

#### Note: 
```
OSC does not provide a `mafft` module (attempting `module load mafft`
results in an "unknown module" error). Therefore, following the same approach
used in the course for BLAST+, MAFFT is installed in the existing `blast_env`
Conda environment:

- load `miniconda3/24.1.2-py310`
- activate `blast_env`
- install MAFFT from the bioconda channel (only once)

The Slurm script `scripts/mafft_bahd.slurm` then loads `miniconda3`,
activates `blast_env`, and runs `mafft --auto`.
```

**How to run:**

```bash
cd /fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final
sbatch scripts/mafft_bahd.slurm

Output (intermediate):

03_alignments/BAHD_clade3_plus_DsTS_aligned.fa

03_alignments/mafft_bahd_<jobid>.out (Slurm log)
```

#### 4.2.2 Tree construction and visualization in RStudio

Script: scripts/plot_bahd_tree.R

This R script (run in the OSC RStudio environment):

1. Uses seqinr::read.alignment() to read the MAFFT alignment.

2. Uses ape::dist.alignment() to compute a distance matrix (1 − identity).

3. Builds a neighbor-joining tree with ape::nj().

4. Saves the tree in Newick format to 04_phylogeny/BAHD_tree.nwk.

5. Plots the tree using base R plot() and saves it to 05_results/BAHD_tree_basic.png.

How to run (in RStudio):
```r
setwd("/fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final")
library(seqinr)
library(ape)
source("scripts/plot_bahd_tree.R")
```
### Outputs:

- Intermediate:

    04_phylogeny/BAHD_tree.nwk

- Final (used in the report):

    05_results/BAHD_tree_basic.png

### 4.3 SNP analysis

Planned steps:

- Create a FASTA file with TS and TS-like CDS from A. belladonna and D. stramonium.

- Use MAFFT to align the CDS and save the alignment in 03_snp/.

- Run snp-sites on the aligned CDS to call SNPs and generate a VCF file in 03_snp/.

- Use R to summarize SNP counts per gene and visualize the results, saving plots/tables in 05_results/.

### 4.3 SNP analysis of TS and TS-like CDS

**Goal.**  
Compare coding sequences of the *Atropa belladonna* TS isoforms and the
putative *Datura stramonium* TS candidate, and summarize SNP variation
across these CDS.

**Input files (00_raw/):**

- `Atropa_TS.cds.fa` — four isoforms of the *A. belladonna* TS CDS.
- `Datura_TS_candidate.cds.fa` — CDS of the *D. stramonium* TS-like candidate.
- `TS_like_cds.fa` — combined FASTA (Atropa TS isoforms + Datura TS candidate).

---

#### 4.3.1 MAFFT alignment of TS-like CDS (OSC bash)

**Script:** `scripts/mafft_ts_cds.slurm`

This Slurm script:

1. Loads the Conda `blast_env` environment (which contains MAFFT).
2. Aligns all sequences in `00_raw/TS_like_cds.fa` with `mafft --auto`.
3. Saves the aligned CDS to `03_snp/TS_like_cds_aligned.fa`.

**How to run:**

```bash
cd /fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final
sbatch scripts/mafft_ts_cds.slurm
```

The intermediate output was:

- 03_snp/TS_like_cds_aligned.fa

#### 4.3.2 SNP calling with snp-sites (OSC bash)

Script: scripts/snp_ts_cds.slurm

This script:

1. Activates the same blast_env environment where snp-sites is installed.

2. Runs snp-sites -v on 03_snp/TS_like_cds_aligned.fa.

3. Writes SNP calls in VCF format to 03_snp/TS_like_snps.vcf.

How to run:
```bash
cd /fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final
sbatch scripts/snp_ts_cds.slurm
```
The intermediate output was:
- 03_snp/TS_like_snps.vcf

#### 4.3.3 R-based SNP summary and visualization (RStudio)

Script: scripts/analyze_snps.R

This script (run in the OSC RStudio environment):

1. Reads 03_snp/TS_like_snps.vcf using readr::read_tsv, skipping meta-lines
starting with ##.

2. Identifies sample columns corresponding to the Atropa TS isoforms and the
Datura TS candidate.

3. Reshapes the genotype matrix into a long format using tidyverse.

4. Counts non-missing genotype calls per sequence.

5. Saves a summary table:

    - 05_results/TS_like_snp_counts.tsv

6. Creates a barplot of SNP counts per sequence and saves it to:

    - 05_results/TS_like_snp_counts.png

How to run in RStudio:
```r
setwd("/fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final")
library(tidyverse)
source("scripts/analyze_snps.R")
```

### 4.4 R-based analysis of BLAST results

Planned R scripts (stored in scripts/):

- analyze_blast.R — read the BLAST output, find the top hit for each query, and visualize percent identity and bit scores (e.g. barplots) and save results in 05_results/.

- plot_tree.R — read the BAHD phylogenetic tree and produce annotated tree plots saved in 05_results/.

- analyze_snps.R — read SNP calls from snp-sites, summarize variation across TS and TS-like genes, and save figures/tables in 05_results/.

**Goal.**  
Summarize BLASTP results and visualize the similarity between clade 3 BAHD proteins and the *Datura stramonium* proteome.

**Input:**

- `02_blast/blast_clade3_vs_Ds_proteome.tab` — BLASTP output in outfmt 6 format.

**Script:**

- `scripts/analyze_blast.R`

This script (run in the OSC RStudio environment used in the course):

1. Loads the `tidyverse` package.
2. Reads the BLAST table with column names (`qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`).
3. For each query sequence, selects the top hit based on the highest bitscore.
4. Saves a summary table to `05_results/blast_top_hits_summary.tsv`.
5. Creates a barplot of percent identity for the top hits and saves it to `05_results/blast_top_hits.png`.

**How to run (in RStudio):**

```r
setwd("/fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final")
library(tidyverse)
source("scripts/analyze_blast.R")
```

## To-do list:

 1. Organize raw files in 00_raw/.

 2.  Implement and run BLAST analysis (scripts/run_blast.sh).

 3. Create multiple sequence alignments for BAHD clade 3 and TS-like CDS.

 4. Build a phylogenetic tree and visualize it in R.

 5. Run snp-sites on aligned TS-like CDS and summarize SNPs in R.

 6. Save final plots and tables to 05_results/.

 7. Update this protocol and notebook.md during the analysis.

 8. Write the final report using results from 05_results/.