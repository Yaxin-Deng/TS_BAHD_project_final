#!/bin/bash
#SBATCH --job-name=snp_ts_cds
#SBATCH --account=PAS2880
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=03_snp/snp_ts_cds_%j.out

# snp_ts_cds.slurm
# Call SNPs from aligned TS and TS-like CDS using snp-sites.
#
# Input:
#   03_snp/TS_like_cds_aligned.fa
#
# Output:
#   03_snp/TS_like_snps.vcf

module load miniconda3/24.1.2-py310
source activate blast_env   # 环境里有 snp-sites

mkdir -p 03_snp

snp-sites -v \
  -o 03_snp/TS_like_snps.vcf \
  03_snp/TS_like_cds_aligned.fa
