
# Project Description
This repository contains scripts for processing SRA data and preparing it for downstream single-cell multiomics analysis using Seurat. The workflow includes **downloading SRA data**, **converting it to FASTQ format**, **compressing the FASTQ files**, **generating input files for Cell Ranger**, and **running Cell Ranger** to generate gene expression and chromatin accessibility data.

# Workflow Steps

## 1. Download SRA Data

The script downloads SRA data specified in the `SRR_Acc_List.txt` file and saves it in the 01_sra directory.
```
# Download SRA data to the 01_sra directory
prefetch --option-file SRR_Acc_List.txt -c --max-size 200G -O ./01_sra
```
## 2. Convert SRA to FASTQ Format
The script converts SRA files to FASTQ format and saves them in corresponding folders within the `02_fastq` directory.
```
# Define SRA ID array
SRA_ids=("SRR18907481" "SRR18907482" "SRR18907480")

# Convert SRA to FASTQ and compress the files
for SRA_id in "${SRA_ids[@]}"
do
    fasterq-dump "./01_sra/$SRA_id" --split-files --include-technical --threads 128 -O "./02_fastq/$SRA_id"
    echo "$SRA_id sra2fastq completed at: $(date +"%Y-%m-%d %H:%M:%S")"
done

for SRA_id in "${SRA_ids[@]}"
do
    pigz -p 128 "./02_fastq/$SRA_id"/*
    echo "$SRA_id fastq2gz completed at: $(date +"%Y-%m-%d %H:%M:%S")"
done
```
## 3. Generate CSV File for Cell Ranger Input

The script creates a CSV file required by Cell Ranger to specify the paths to the FASTQ files and their corresponding sample and library types.
```
fastq_path="/your_path/02_fastq"

# Define SRA ID arrays
SRA_GEX=("SRR18907480")
SRA_ATAC=("SRR18907481" "SRR18907482")

# Output header to CSV file
echo "fastqs,sample,library_type" > ./03_cellranger/libraries.csv

# Output GEX and ATAC data to CSV file
for sample in "${SRA_GEX[@]}"; do
    echo "$fastq_path/$sample,$sample,Gene Expression" >> ./03_cellranger/libraries.csv
done

for sample in "${SRA_ATAC[@]}"; do
    echo "$fastq_path/$sample,$sample,Chromatin Accessibility" >> ./03_cellranger/libraries.csv
done
```

## 4. Run Cell Ranger
The script runs Cell Ranger with the specified reference genome and libraries CSV file to generate the multiomics data.
```
cd ./03_cellranger

sample_name="your_sample"
ref_genome="/your_software_path/refdata-cellranger-arc-mm10-2020-A-2.0.0" # mm10
# ref_genome="/your_software_path/refdata-cellranger-arc-GRCh38-2020-A-2.0.0" # hg38

cellranger-arc count --id=${sample_name} \
                     --reference=${ref_genome} \
                     --libraries=./libraries.csv \
                     --localcores=128 \
                     --localmem=128
```
## Usage
#### 1. Update the `your_path` and `your_sample` in the script as necessary.

#### 2. Execute the script step by step to process the SRA data and prepare it for Seurat analysis.
