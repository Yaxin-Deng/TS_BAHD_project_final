#!/bin/bash
#SBATCH --job-name=mafft_bahd
#SBATCH --account=PAS2880
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=03_alignments/mafft_bahd_%j.out

# mafft_bahd.slurm
# Multiple sequence alignment of clade 3 BAHD protein sequences
# including the Atropa TS and Datura TS-like sequences.
#
# Input:
#   00_raw/BAHD_clade3_ref_plus_DsTS.fa
#
# Output (intermediate):
#   03_alignments/BAHD_clade3_plus_DsTS_aligned.fa

# 加载 miniconda 并激活 blast_env 环境（课堂上使用过的方式）
module load miniconda3/24.1.2-py310
source activate blast_env   # 或 conda activate blast_env，如果你那边是这种写法

mkdir -p 03_alignments

mafft --auto 00_raw/BAHD_clade3_ref_plus_DsTS.fa \
  > 03_alignments/BAHD_clade3_plus_DsTS_aligned.fa
