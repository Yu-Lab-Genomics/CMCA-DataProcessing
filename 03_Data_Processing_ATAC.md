# Single-Cell Multi-Omics Integration Pipeline

This repository hosts a robust analysis pipeline designed to integrate single-cell RNA sequencing (scRNA-seq) and Assay for Transposase-Accessible Chromatin with high-throughput sequencing (ATAC-seq) data using Seurat, a popular R package for single-cell genomics analysis.
rm(list = ls())
set.seed(2024)

## 1. Data Input
```
# Set working directory
setwd("/file_path/sampel_ID/")
getwd()

# Set variable
GENOME <- "hg38"
counts_path <- "filtered_feature_bc_matrix/"
frag_path <- "sample_atac_fragments.tsv.gz"
mac2_path <- "/your_path/macs2"
output_path.prefix <- "/your_path/output/"

# Load necessary packages
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
#library(EnsDb.Mmusculus.v79) #mm10
#library(BSgenome.Mmusculus.UCSC.mm10) #mm10
library(EnsDb.Hsapiens.v86) # hg38
library(BSgenome.Hsapiens.UCSC.hg38) #hg38

# Set genome-specific variables
if (GENOME == "hg38") {
  ref_genome <- EnsDb.Hsapiens.v86
  mt_prefix <- "^MT-"
  peak.blacklist <- blacklist_hg38_unified
  print("genome=hg38")
} else if (GENOME == "mm10") {
  ref_genome <- EnsDb.Mmusculus.v79
  mt_prefix <- "^mt-"
  peak.blacklist <- blacklist_mm10
  print("genome=mm10")
} else {
  stop("genome not defined")
}
```
## 2. Data Process
```
# Load the RNA and ATAC data
counts <- Read10X(counts_path)
fragpath <- frag_path

# Get gene annotations
annotation <- GetGRangesFromEnsDb(ensdb = ref_genome)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# Create Seurat object containing the RNA data
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  min.cells = 3)

# Add ATAC assay to Seurat object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation)

# Perform quality control
DefaultAssay(pbmc) <- "RNA"
mt.genes <- rownames(pbmc)[grep(mt_prefix,rownames(pbmc))]
QC <- GetAssayData(object = pbmc, layer = "counts")
percent.mito <- Matrix::colSums(QC[mt.genes,])/Matrix::colSums(QC) * 100
pbmc <- AddMetaData(pbmc, percent.mito, col.name = "percent.mt")

DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

# Filter low-quality cells
cell.number1 <- dim(pbmc)
#[1] 124914   3361
pbmc_filted <- subset( 
  x = pbmc,
  subset = nCount_RNA > 500 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    nFeature_RNA > 100 &
    percent.mt < 40)
cell.number2 <- dim(pbmc_filted)
# [1] 111857  11546

# Call peaks using MACS2
pbmc <- pbmc_filted
DefaultAssay(pbmc) <- "ATAC"
peaks <- CallPeaks(pbmc, macs2.path = mac2_path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = peak.blacklist, invert = TRUE)

# Quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# Create and add peak assay to Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# Calculate gene activity scores
DefaultAssay(pbmc) <- 'peaks'
# Create gene activity matrix
gene.activity <- GeneActivity(pbmc)
# Add the gene activity matrix to the Seurat object
pbmc[['Activity']] <- CreateAssayObject(counts = gene.activity)
# Normalize data
pbmc <- NormalizeData(
object = pbmc,
assay = 'Activity',
normalization.method = 'LogNormalize',
scale.factor = median(pbmc$nCount_Activity)
)

# Perform RNA UMAP
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
nPC=1:20
pbmc <- FindNeighbors(pbmc, dims = nPC) %>%
  RunUMAP(dims = nPC)
pbmc = FindClusters(pbmc, resolution = 0.1)

# Perform ATAC UMAP
DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
nPC = 2:30
pbmc <- RunUMAP(pbmc,reduction = 'lsi', dims = nPC, reduction.name = "lsiumap") %>%
  FindNeighbors(reduction = 'lsi', dims = nPC)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

```
## 3. Data Output
```
# Save Seurat object to RDS
saveRDS(pbmc, file = paste0(output_path.prefix,"sample.rds"))

# Export RNA data to CSV
pbmc.rna <- as.data.frame(pbmc@assays$SCT@data) # RNA
pbmc.atac <- as.data.frame(pbmc@assays$peaks@data) # ATAC
pbmc.activity <- as.data.frame(pbmc@assays$Activity@data) # activity score

umap.rna <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
umap.atac <- as.data.frame(pbmc@reductions$lsiumap@cell.embeddings)

pbmc.ct <- Idents(pbmc)

df.umap_rna <- cbind(umap.rna, pbmc.ct)
colnames(df.umap_rna) <- c("x","y","celltype")
df.umap_atac <- cbind(umap.atac, pbmc.ct)
colnames(df.umap_atac) <- c("x","y","celltype")

df.trans_rna <- as.data.frame(t(pbmc.rna))
df.trans_atac <- as.data.frame(t(pbmc.atac))
df.trans_activity <- as.data.frame(t(pbmc.activity))

df.merged_umap.rna <- merge(df.umap_rna, df.trans_rna, by = "row.names", all.x = TRUE) %>%
  dplyr::rename(name = "Row.names")
fwrite(df.merged_umap.rna, paste0(output_path.prefix,"sample-RNA.csv"), row.names = FALSE)

df.merged_umap.atac <- merge(df.umap_atac, df.trans_atac, by = "row.names", all.x = TRUE) %>%
  dplyr::rename(name = "Row.names")
fwrite(df.merged_umap.atac, paste0(output_path.prefix,"sample-ATAC.csv"), row.names = FALSE)

df.merged_umap.activity <- merge(df.umap_rna, df.trans_activity, by = "row.names", all.x = TRUE) %>%
  dplyr::rename(name = "Row.names")
fwrite(df.merged_umap.activity, paste0(output_path.prefix,"sample-Activity.csv"), row.names = FALSE)
```
## Usage
Set the working directory and sample in the script.
