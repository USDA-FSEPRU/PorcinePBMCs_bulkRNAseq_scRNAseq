---
Title:  "Porcine immune cell atlas single cell RNA sequencing analysis pipeline"
Author: "Sathesh K. Sivasankaran"
Date:   "Feb 15, 2021"
---

# Porcine Immune Single Cell Atlas

## Doublet cells filtering.
progressed cells and genes from each sample were analyzed using the Python package, Scrublet (v0.2; https://github.com/AllonKleinLab/scrublet), to remove highly probable neotypic doublets. To estimate the doublet rate of each cell Scrublet was accessed via JupyterHub. For each sample, synthetic doublets were predicted using the  doublet rate score various from 0.060 to 0.080 with 0.005 intervals. Jupyter notebook scripts to estimate doublet rate score for each sample was located in:\
![A](Notebook/Doublet/A.ipynb)\
![B](Notebook/Doublet/B.ipynb)\
![C](Notebook/Doublet/C.ipynb)\
![CT21NOV](Notebook/Doublet/CT21NOV.ipynb)\
![CT230OCT](Notebook/Doublet/CT230OCT.ipynb)\
![PBMC1](Notebook/Doublet/PBMC1.ipynb)\
![PBMC2](Notebook/Doublet/PBMC2.ipynb)

A Scrublet object was created using doublet rate of 0.070. Doublet scores were calculated, and a doublet score threshold of 0.25 was set based on a bimodal distribution observed between embedded and neotypic doublets.\
![DoubletEstimated](Notebook/Doublet/DoubletEstimated.ipynb)

Finally, a lists of doublet probability scores were generated and imported into R, where cells with corresponding doublet probability scores grater than 0.25 were removed from the dataset. Filtered count matrix for each sample was outputted into CellRanger format using DropletUtils.

Load required packages into R.
```
pacman::p_load(dplyr, DropletUtils, Seurat)
setwd ("C:/Your/Filtered/count/matrix/location/")
```

Load data into R and filter poor quality cells from the datset.
```
All <- readRDS("Filtered/PBMCQC.rds")
All$counts <- All$counts[,All$phenoData$PassAll]
All$phenoData <- All$phenoData[All$phenoData$PassAll,]
```

Import the estimated Scrublet scores into R and merge scores into cell feature (phenoData or cDat). Finally, cells with scrublet score <0.25 removed form the dataset.
```
A <- read.csv("Scrublet/A_ScrubScore.csv")
B <- read.csv("Scrublet/B_ScrubScore.csv")
C <- read.csv("Scrublet/C_ScrubScore.csv")
CT21NOV <- read.csv("Scrublet/CT21NOV_ScrubScore.csv")
CT230OCT <- read.csv("Scrublet/CT230OCT_ScrubScore.csv")
PBMC1 <- read.csv("Scrublet/PBMC1_ScrubScore.csv")
PBMC2 <- read.csv("Scrublet/PBMC2_ScrubScore.csv")

scrub_combined <- rbind(A, B, C, CT21NOV, CT230OCT, PBMC1, PBMC2)
All$phenoData$Scrublet <- scrub_combined$X0
All$phenoData <- mutate(All$phenoData, PassScrub = Scrublet < 0.25)
rownames(All$phenoData) <- All$phenoData$Loupe

All$counts <- All$counts[,All$phenoData$PassScrub]
```

Again, gene filtration step was performed here after all cell filtering steps to remove genes that showed no expression in the dataset.
```
rownames(All$featureData) <- rownames(All$counts)
keep <- rowSums(All$counts) > 0
All$counts <- All$counts[keep,]

All$phenoData <- All$phenoData[All$phenoData$PassScrub,]
All$featureData <- All$featureData[rownames(All$counts),]
```
Finally, dataset was save in RDS format and also count matrix (*.mtx and *.tsv) for downstream analyses including integration, clustering and differential gene expression.
```
stopifnot(identical(as.character(rownames(All$phenoData)),colnames(All$counts)))
stopifnot(identical(as.character(All$featureData$Name),rownames(All$counts)))
out <- list()
out[["counts"]] <- All$counts
out[["phenoData"]] <- All$phenoData
out[["featureData"]] <- All$featureData
saveRDS(out,file=file.path("Scrublet/", "PBMCScrubbedDFs.rds"))

seurat <- CreateSeuratObject(counts = All$counts, meta.data = All$phenoData)
write10xCounts(x = seurat@assays$RNA@counts, path = "Scrublet/PBMCScrubbedCount", version = "3")
saveRDS(seurat,file=file.path("Scrublet/", "PBMCScrubbedSeurat.rds"))
```

### All filtering steps were ended here.
## ___END___
