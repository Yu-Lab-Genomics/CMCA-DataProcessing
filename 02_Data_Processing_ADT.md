# Single-Cell RNA and ADT Analysis Pipeline
This script performs comprehensive preprocessing and analysis of single-cell RNA sequencing (scRNA-seq) and Antibody-Derived Tag (ADT) data using the Seurat package in R. 
## 1. Data Input
Set the working directory and output path
```
rm(list = ls())
set.seed(2024)

library(Seurat)
library(hdf5r)
library(tidyverse)
library(data.table)

# Set working directory
setwd("/file_path/cmcaID/sampleID")
getwd()

# Set variable
output_path.prefix <- '/file_path/output'

filename <- "sampleID_filtered_feature_bc_matrix.h5"
```
Load the RNA and ADT data
```
counts <- Read10X_h5(filename)
```

## 2. Seurat Object Creation

Create Seurat object containing the RNA data
```
mt <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  min.cells = 3)
```
Add ADT assay to Seurat object
```
adt <- CreateAssayObject(counts$`Antibody Capture`)
mt[["ADT"]] <- adt
```
Calculate the mitochondrial gene percentage for quality control.
```
DefaultAssay(mt) <- "RNA"
mt[["percent.mt"]] <- PercentageFeatureSet(mt, pattern = "^mt-")
```
## 3. Quality Control
Filter out low-quality cells based on specific thresholds for RNA counts, feature counts, and mitochondrial gene percentage.
```
cell.number1 <- dim(mt)
mt_filted <- subset( 
  x = mt,
  subset = nCount_RNA > 500 &
    nCount_RNA < 25000 &
    nFeature_RNA > 100 &
    percent.mt < 40)
cell.number2 <- dim(mt_filted) # sample cell numbers
mt <- mt_filted
```
## 4. RNA Data Processing
```
DefaultAssay(mt) <- "RNA"
# Normalize, find variable features, and scale RNA data
nPC = 1:20
mt <- NormalizeData(mt)
mt <- FindVariableFeatures(mt, nfeatures = 2000)
mt <- ScaleData(mt)
mt <- RunPCA(mt, dims = nPC)

# Find neighbors and clusters, and run UMAP for RNA data
mt <- FindNeighbors(mt, dims = nPC)
mt <- FindClusters(mt, resolution = 0.9, verbose = FALSE)
mt <- RunUMAP(mt, dims = nPC)
```

## 5. ADT Data Processing
```
# Normalize ADT data
DefaultAssay(mt) <- "ADT"
mt <- NormalizeData(mt, normalization.method = "CLR", margin = 2)
```
## 6. Data Output
Save Seurat object to RDS
```
saveRDS(mt, file = paste0(output_path.prefix,"sample.rds"))
```
Export data to CSV
```
mt.rna <- as.data.frame(mt@assays$RNA@layers$data) # RNA
rownames(mt.rna) <- rownames(mt@assays$RNA)
colnames(mt.rna) <- colnames(mt@assays$RNA)

mt.adt <- as.data.frame(mt@assays$ADT@data) # ADT

umap.rna <- as.data.frame(mt@reductions$umap@cell.embeddings)

mt.ct <- Idents(mt)

df.umap_rna <- cbind(umap.rna, mt.ct)
colnames(df.umap_rna) <- c("x","y","celltype")

df.trans_rna <- as.data.frame(t(mt.rna))
df.trans_adt <- as.data.frame(t(mt.adt))

df.merged_umap.rna <- merge(df.umap_rna, df.trans_rna, by = "row.names", all.x = TRUE) %>%
  dplyr::rename(name = "Row.names")
fwrite(df.merged_umap.rna, paste0(output_path.prefix,"cmcaID-sampleID-RNA.csv"), row.names = FALSE)

df.merged_umap.adt <- merge(df.umap_rna, df.trans_adt, by = "row.names", all.x = TRUE) %>%
  dplyr::rename(name = "Row.names")
fwrite(df.merged_umap.adt, paste0(output_path.prefix,"cmcaID-sampleID-ADT.csv"), row.names = FALSE)
```

## Usage
- Set the working directory and output path: Update the `setwd` and `output_path.prefix` variables with your file paths.
