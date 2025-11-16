# Final Project Proposal

**Student:** Yaxin Deng  
**Course:**  PLNTPTH 5006 - Comput Omics Data 
**Project Title:** BLAST Analyses and SNP Identification in Available *Atropa* BAHD Gene Sequences  
**Date:** November 15, 2025

---

## 1. Background and Motivation

My PhD research focuses on the biosynthesis of tropane alkaloids in Solanaceae species. A key enzymatic step involves BAHD acyltransferases, which show notable sequence diversity among *Atropa* accessions and isoforms.  

My advisor specifically asked me to complete:

1. **BLAST analyses**  
2. **SNP identification in the available *Atropa* sequences**

This Final Project will build a fully reproducible, OSC-based workflow to accomplish these two tasks using the tools taught in this course (shell scripting, Slurm, Git, reproducibility principles).  
The project also directly supports my ongoing PhD research on tigloyltransferase-like BAHD enzymes.

---

## 2. Data and Reference Sequences

### **Input data**
- *Atropa belladonna* BAHD candidate CDS/protein sequences from public genome assemblies
- Lab-provided *Datura stramonium* tigloyltransferase ortholog (ANN22971)
- Published BAHD sequences, including all sequences in **Fig. 2E** of  
  *de la Cruz et al. (2024)*  
  (*Discovering a mitochondrion-localized BAHD acyltransferase involved in calystegine biosynthesis and engineering the production of 3β-tigloyloxytropane*)

### **Expected outputs**
- BLAST similarity tables (`.tab` + headered `.tsv`)  
- MAFFT alignment files (`.fa`)  
- SNP summary table (`.tsv`), including synonymous vs nonsynonymous classification  
- Final Markdown report and reproducible Slurm scripts

---

## 3. Use of BAHD Reference Sequences (Fig. 2E)

An essential component of this project is the curated reference set of BAHD acyltransferase proteins taken from **Fig. 2E of de la Cruz et al. (2024)**. These enzyme sequences represent validated BAHD family members across multiple clades and contain canonical motifs such as **HXXXD** and **DFGWG**.

### **Why these sequences are included**
These sequences form a **high-quality functional reference catalog** required for BLAST analyses. They allow accurate classification of *Atropa* and *Datura* candidate sequences by:

- Providing known functional BAHD family representatives  
- Anchoring BLAST searches in well-characterized enzymatic clades  
- Supporting functional inference about tigloyltransferase-like activity  

### **Where these sequences are used**
1. **BLAST Database Construction**
   - Combined into `05_bahd_catalog/BAHD_ref.fa`
   - Used to build a local protein BLAST database:

     ```bash
     makeblastdb -in BAHD_ref.fa -dbtype prot
     ```

2. **Homology Search**
   All *Atropa* and *Datura* BAHD candidates are BLASTed against this database to determine:
   - closest homologs  
   - % identity  
   - alignment quality  
   - predicted functional similarity  

### **Where they are NOT used**
The Fig. 2E sequences are **not used in SNP analysis**.  
SNP calling is performed **only among *Atropa* isoforms** of the same gene because SNP analysis must compare allelic variants within a species.

---

## 4. Planned Directory Structure

```plaintext
TS_BAHD_project_final/
├── 00_raw/              # Raw FASTA files (reference + Atropa sequences)
├── 01_qc/               # QC and translations
├── 02_blast/            # BLAST database + BLAST output
├── 03_msa/              # MAFFT alignments
├── 04_snp/              # snp-sites TSV outputs
├── 05_reports/          # Plots, tables, README documentation
└── scripts/
    ├── run_blast.sh
    ├── run_mafft.sh
    ├── run_snp.sh
    └── slurm/
        ├── blast_job.sbatch
        ├── mafft_job.sbatch
        └── snp_job.sbatch
