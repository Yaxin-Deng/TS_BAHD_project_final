# Project notebook
## 2025-12-01

- Created new project directory `TS_BAHD_project_final` on OSC.
- Set up subdirectories: `00_raw`, `02_blast`, `03_alignments`, `03_snp`, `04_phylogeny`, `05_results`, `scripts`.
- Copied raw FASTA files from the older project into `00_raw/` and renamed them:
  - `ANN02033_cds_JP.fna.txt` → `Datura_TS_candidate.cds.fa`
  - `DS_TS_ANN02033.prot.fa` → `Datura_TS_candidate.protein.fa`
  - `Clade3_Protein.fa` → `BAHD_clade3_ref_plus_DsTS.fa`
  - `DS_protein.faa` → `Datura_genome_protein.fa`
- Created `protocol.md` and added sections on project overview, directory structure, and BLAST analysis.

## 2025-12-02

- Wrote BLAST script `scripts/run_blast.sh` (builds Datura protein database and runs BLASTP with clade 3 BAHD queries).
- Submitted the BLAST job on OSC:

```bash
cd /fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final
sbatch scripts/run_blast.sh
```
Next step: check that 02_blast/blast_clade3_vs_Ds_proteome.tab was created and start writing scripts/analyze_blast.R for R-based summaries.

## 2025-12-04

- Opened OSC RStudio and set the working directory:

```r
setwd("/fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final")
```

- Created an R script scripts/analyze_blast.R to process the BLASTP output
(02_blast/blast_clade3_vs_Ds_proteome.tab). The script:

- reads the BLAST outfmt 6 table,

- selects the top hit (highest bitscore) for each query BAHD protein,

- saves a summary table,

- creates a barplot of percent identity.

In RStudio, I run:

```r
library(tidyverse)
source("scripts/analyze_blast.R")
```
New result files were generated in 05_results/:

blast_top_hits_summary.tsv

blast_top_hits.png

## 2025-12-05 (MAFFT and phylogeny)

- Wrote a Slurm script `scripts/mafft_bahd.slurm` to align clade 3 BAHD protein
  sequences (including Atropa TS and Datura TS-like proteins) with MAFFT.
- Submitted the job on OSC:

```bash
cd /fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final
sbatch scripts/mafft_bahd.slurm
```
- The MAFFT alignment was saved to:

  - 03_alignments/BAHD_clade3_plus_DsTS_aligned.fa

- In OSC RStudio, set the working directory:

```
setwd("/fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final")
```
- Created an R script scripts/plot_bahd_tree.R using the seqinr and ape
packages to:

  - read the MAFFT alignment,

  - compute a distance matrix (dist.alignment),

  - build a neighbor-joining tree (nj),

  - save the tree as Newick, and

  - plot the tree to a PNG file.

- Ran the script in RStudio:
```r
library(seqinr)
library(ape)
source("scripts/plot_bahd_tree.R")
```
New files generated:

  - 04_phylogeny/BAHD_tree.nwk

  - 05_results/BAHD_tree_basic.png

## 2025-12-05 (fixing MAFFT alignment)

- Submitted `scripts/mafft_bahd.slurm` and noticed that the output file
  `03_alignments/BAHD_clade3_plus_DsTS_aligned.fa` was empty.
- I Checked the Slurm log file `03_alignments/mafft_bahd_42370166.out`, which showed:

  - `Lmod has detected the following error: The following module(s) are unknown: "mafft"`
  - `/var/spool/slurmd/...: line 22: mafft: command not found`

  This indicated that there is no `mafft` module on OSC.

- Solution:
  - Loaded miniconda and activated the existing `blast_env` environment.
  - Installed MAFFT into that environment:

```bash
module load miniconda3/24.1.2-py310
source activate blast_env
conda install -y -c bioconda mafft
```

- Updated scripts/mafft_bahd.slurm to:

  - load miniconda3/24.1.2-py310

  - activate blast_env

  - then run mafft --auto ...

- Re-submitted the job:

```bash
cd /fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final
sbatch scripts/mafft_bahd.slurm
```

The output was:
```
Submitted batch job 42370364
```
After the job finished, the file 03_alignments/BAHD_clade3_plus_DsTS_aligned.fa
contained the expected multiple sequence alignment.

### 2025-12-07 (Fixing wrong alignment filename)

- When running `scripts/plot_bahd_tree.R`, the script returned an error because
  the alignment filename was still set to an old name from my previous,
  disorganized project (`03_alignments/BAHD_ref_protein_aligned.fa`).

- After checking the new project folder structure, I confirmed the correct file is:

  `03_alignments/BAHD_clade3_ref_plus_DsTS_aligned.fa`

- Updated the script to:

```r
align_file <- "03_alignments/BAHD_clade3_ref_plus_DsTS_aligned.fa"
```

- Re-ran the script (via Rscript due to RStudioGD error) and successfully generated:

  - 04_phylogeny/BAHD_tree.nwk

  - 05_results/BAHD_tree_basic.png

### 2025-12-07 SNP analysis of TS and TS-like CDS

#### Preparing CDS input
  - Combined Atropa belladonna TS CDS (4 isoforms) and the Datura stramonium TS candidate CDS into a single FASTA file:
```bash
cat 00_raw/Atropa_TS.cds.fa 00_raw/Datura_TS_candidate.cds.fa > 00_raw/TS_like_cds.fa
```
#### MAFFT alignment of TS-like CDS
Submitted MAFFT alignment job using the Slurm script

scripts/mafft_ts_cds.slurm:
```bash
sbatch scripts/mafft_ts_cds.slurm
```
The output was:
```
Submitted batch job 42370747

03_snp/TS_like_cds_aligned.fa
```
SNP calling using snp-sites

Installed snp-sites into the course Conda environment (blast_env):

```bash
module load miniconda3/24.1.2-py310
source activate blast_env
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels defaults
conda install -y snp-sites
```
Submitted SNP calling job:
```
sbatch scripts/snp_ts_cds.slurm
```

The output was:
```
Submitted batch job 42370753

03_snp/TS_like_snps.vcf
```

## R analysis and visualization

I run R script scripts/analyze_snps.R in OSC RStudio:
```r
setwd("/fs/ess/PAS2880/users/dengyaxin1156/TS_BAHD_project_final")
library(tidyverse)
source("scripts/analyze_snps.R")
```
New result files were:

- 05_results/TS_like_snp_counts.tsv

- 05_results/TS_like_snp_counts.png