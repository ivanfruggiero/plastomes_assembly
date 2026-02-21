# plastomes_assembly

Modular and reproducible workflow for plastome assembly from paired-end reads, integrating quality control, adapter trimming, de novo organelle assembly (GetOrganelle), and anchor-based sequence normalization.

---

## Overview

This repository provides a minimal, production-style workflow designed to be portable across systems.

The focus is on clean scripting practices:
- explicit inputs and outputs
- structured logging
- parameterization via config files
- robust sample discovery
- portability across environments (local, HPC, conda, modules)

---

## Workflow Steps

1. QC (pre-trimming) – FastQC on raw reads  
2. Trimming – Trimmomatic paired-end trimming (parameters should be adapted based on QC results)  
3. Assembly – GetOrganelle plastid assembly from trimmed reads  
4. Normalization – Anchor-based rotation to match a reference start (single-record or quadripartite reference)

---

## Directory Layout

plastomes_assembly/
├── config/
│   ├── config.example.env
│   └── config.env              # local only (not tracked)
├── data/
│   └── raw/                    # raw FASTQ (not tracked)
├── scripts/
│   ├── 01_qc.sh
│   ├── 02_trimming.sh
│   ├── 03_getorganelle.sh
│   └── 04_normalize_plastomes.py
├── results/
└── logs/

---

## Requirements

- FastQC  
- Trimmomatic (jar) + adapter FASTA  
- GetOrganelle  
- Python 3  

Tools can be installed system-wide, via conda/mamba, or via module systems on HPC.  
The workflow assumes executables are available in PATH unless configured otherwise.

---

## Setup

Copy the example config and edit it locally:

cp config/config.example.env config/config.env
nano config/config.env

Place raw paired-end reads in data/raw/ using the naming convention:

SAMPLE_R1_001.fastq.gz  
SAMPLE_R2_001.fastq.gz  

---

## Running the Workflow

1) QC (pre-trimming)

bash scripts/01_qc.sh --stage pre

2) Trimming

bash scripts/02_trimming.sh

3) QC (post-trimming)

bash scripts/01_qc.sh --stage post

4) GetOrganelle assembly

bash scripts/03_getorganelle.sh

If GetOrganelle is installed in a conda environment:

bash scripts/03_getorganelle.sh --conda --conda-env getorganelle

5) Normalization (anchor-based)

Provide a reference FASTA and an input folder containing single-record FASTA files to normalize.

The reference can be:
- a single-record FASTA, or  
- a quadripartite FASTA with 4 records (e.g., LSC, IRb, SSC, IRa), which will be concatenated in file order.

Example:

python3 scripts/04_normalize_plastomes.py \
  --ref data/reference/reference_quadripartition.fasta \
  --in_dir results/03_assembly/getorganelle_fastas \
  --out_dir results/04_normalized \
  --pattern "*.fasta"

A normalization_report.tsv file is produced in the output directory.

---

## Notes on Trimming Parameters

Trimming settings are dataset-dependent.

Use FastQC (and optionally MultiQC) to determine whether additional steps such as:

- HEADCROP  
- LEADING  
- TRAILING  
- stricter SLIDINGWINDOW  

are required.

This workflow provides a generic starting configuration and should be adapted to the characteristics of the sequencing dataset.
