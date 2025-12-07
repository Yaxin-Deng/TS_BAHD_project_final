#!/bin/bash
#SBATCH --job-name=mafft_ts_cds
#SBATCH --account=PAS2880
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --output=03_snp/mafft_ts_cds_%j.out

# mafft_ts_cds.slurm
# Multiple sequence alignment of Atropa TS isoforms and the Datura TS candidate CDS.
#
# Input:
#   00_raw/TS_like_cds.fa   (Atropa_TS.cds.fa + Datura_TS_candidate.cds.fa)
#
# Output (intermediate):
#   03_snp/TS_like_cds_aligned.fa

module load miniconda3/24.1.2-py310
source activate blast_env   # 环境里已经安装了 mafft

mkdir -p 03_snp

mafft --auto 00_raw/TS_like_cds.fa \
  > 03_snp/TS_like_cds_aligned.fa
