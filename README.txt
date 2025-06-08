# README.txt

Author: Rebecca Key
Adapted from Bryndan Durham

University of Florida
Department of Biology
Durham Lab

## Overview
This folder contains a series of bash scripts designed to automate phylogenetic placement of short-read metatranscripts using pplacer. The pipeline consists of multiple steps, beginning with building reference packages, aligning reads, running pplacer, and summarizing outputs. This is catered to perform using parallel computation, where we provide a list of genes to reun simultaneously, rather than linearly. Plus, we are operating in the HiperGator Cluster. 

## Purpose
The goal is to identify transcripts in both NS (non-selected) and PA (particle-associated) seawater samples that align to a gene of interest (e.g., an enzyme-coding gene). We then pull annotation information to determine which microbial organisms contain the query gene. This generates a gene dendrogram, and the number of occurrences across sample datasets is tallied using the dendrogram as a reference framework.

---

## Scripts

1. **make_newGeneFolders.sh**  
   - Creates directory structures for storing gene-specific outputs and results.

2. **1_pplacer_A.sh**  
   - Prepares input files: runs sequence info gathering, taxonomy updates, builds reference packages, and constructs HMM profiles.

3. **1_pplacer_B.sh**  
   - Submits pplacer jobs to the cluster for each gene using dataset-specific scripts.

4. **1_pplacer_C_*.sh**  
   - Executes HMM searches, alignments, and runs pplacer for individual samples across datasets (g1ns, g1pa, g2ns, g2pa).

5. **2_pplacer.sh**  
   - Merges pplacer outputs and generates XML and CSV files summarizing placements.

6. **3_pplacer.sh**  
   - Further processes placement results, applying normalization factors and creating final count tables.

---

## Usage
- Submit scripts in order:
  1. `make_newGeneFolders.sh`
  2. `1_pplacer_A.sh`
  3. `1_pplacer_B.sh` (which triggers appropriate 1_pplacer_C_*.sh scripts)
  4. `2_pplacer.sh`
  5. `3_pplacer.sh`

- Each script contains SLURM batch directives for high-performance cluster scheduling.
- Adjust file paths, dataset names, and environment modules as needed for your cluster.

---

## Notes
- Dependencies include: python, seqmagick, taxtastic, hmmer, pplacer.
- Some scripts require dataset-specific input files (e.g., genesAbbr.txt, genes.txt, sample lists).
- Use caution with submodules or nested repositories (e.g., tools/guppy).
- Refer to individual scripts for detailed usage instructions.

---

## Contact
For questions or assistance, please contact:  
Rebecca Key, Ph.D.  
rkeyMicrobe@proton.me

---
