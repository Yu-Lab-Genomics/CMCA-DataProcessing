# CMCA-DataProcessing
This repository offers standardized data processing pipelines for the CMCA database. We hope this could assist users in comprehending the uniform processing applied to the data. 

`01_Data_Processing_SRA.md`
This script provides examples on processing SRA data and preparing it for downstream single-cell multiomics analysis using Seurat. The workflow encompasses downloading SRA data, converting it to FASTQ format, compressing the FASTQ files, generating input files for Cell Ranger, and executing Cell Ranger to produce such as gene expression and chromatin accessibility data.

`02_Data_Processing_ADT.md`
This script provides examples for comprehensive preprocessing and analysis techniques for single-cell RNA sequencing (scRNA-seq) and Antibody-Derived Tag (ADT) data using the Seurat package in R.

`03_Data_Processing_ATAC.md`
This script provides examples for comprehensive preprocessing and analysis of single-cell RNA sequencing (scRNA-seq) and single-cell Assay for Transposase-Accessible Chromatin with high-throughput sequencing (scATAC-seq) data using the Seurat package in R.

The software we used includes:
- cellranger-arc (version 2.0.2)
- SAMtools (version 1.19.2)
- BCFtools (version 1.19)
- HTSlib (version 1.19.1)
- fastq-dump (version 3.1.0)
- R (version 4.3.3)
- MACS2 (version 2.2.9.1)
- Seurat (version 5.0.3)
- Signac (version 1.13.0)


