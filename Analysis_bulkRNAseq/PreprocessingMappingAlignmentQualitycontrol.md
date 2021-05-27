Preprocessing, mapping, alignment, quality control
================
### Background
Data processing was performed as previously reported (Herrera-Uribe et al., 2020) using Sscrofa 11.1 genome and annotation v11.1.97 were used. Counts per gene of each sample in the two count tables were added together to get the final count table. Given that different types of immune cells have different transcriptome profiles (Hicks and Irizarry, 2015), YARN (Paulson et al., 2017), a tissue type-aware RNA-seq data normalization tool, was used to filter and normalize the count table. Genes with extremely low expression levels (<4 counts in at least one cell type) were filtered out using filterLowGenes(). The final count table contained 12,261 genes across 15 samples, which was then normalized using normalizeTissueAware(), which leverages the smooth quantile normalization method (Hicks et al., 2018).
Data quality control was performed using DESeq2 (v1.24.0) (Love et al., 2014) within  RStudio s (v1.2.1335). Regularized log-transformation was applied to the normalized count table with the rld function. Then principal component analysis (PCA) and sample similarity analyses were carried out and visualized using plotPCA() and distancePlot(), respectively. Heatmaps to display enriched genes were created using pheatmap (v1.0.12) within RStudio.

### Preprocessing
``` r
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # 1 processor core(s) per node
#SBATCH --mem=4G   # maximum memory per node
#SBATCH --job-name="cmb.fastq"
#SBATCH --output="logs/out.%j.cmb.fastq" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.%j.cmb.fastq" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

sh combine.R1.fastq.sh
sh combine.R2.fastq.sh
```
``` r
mkdir -p combined.fastq; cat 28Oct19/LIB105109_S1_L001_R1_001.fastq.gz	29Oct19/LIB105109_S1_L001_R1_001.fastq.gz	28Oct19/LIB105109_S1_L002_R1_001.fastq.gz	29Oct19/LIB105109_S1_L002_R1_001.fastq.gz	28Oct19/LIB105109_S1_L003_R1_001.fastq.gz	29Oct19/LIB105109_S1_L003_R1_001.fastq.gz	28Oct19/LIB105109_S1_L004_R1_001.fastq.gz	29Oct19/LIB105109_S1_L004_R1_001.fastq.gz > combined.fastq/LIB105109_S1_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105110_S2_L001_R1_001.fastq.gz	29Oct19/LIB105110_S2_L001_R1_001.fastq.gz	28Oct19/LIB105110_S2_L002_R1_001.fastq.gz	29Oct19/LIB105110_S2_L002_R1_001.fastq.gz	28Oct19/LIB105110_S2_L003_R1_001.fastq.gz	29Oct19/LIB105110_S2_L003_R1_001.fastq.gz	28Oct19/LIB105110_S2_L004_R1_001.fastq.gz	29Oct19/LIB105110_S2_L004_R1_001.fastq.gz > combined.fastq/LIB105110_S2_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105111_S3_L001_R1_001.fastq.gz	29Oct19/LIB105111_S3_L001_R1_001.fastq.gz	28Oct19/LIB105111_S3_L002_R1_001.fastq.gz	29Oct19/LIB105111_S3_L002_R1_001.fastq.gz	28Oct19/LIB105111_S3_L003_R1_001.fastq.gz	29Oct19/LIB105111_S3_L003_R1_001.fastq.gz	28Oct19/LIB105111_S3_L004_R1_001.fastq.gz	29Oct19/LIB105111_S3_L004_R1_001.fastq.gz > combined.fastq/LIB105111_S3_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105112_S4_L001_R1_001.fastq.gz	29Oct19/LIB105112_S4_L001_R1_001.fastq.gz	28Oct19/LIB105112_S4_L002_R1_001.fastq.gz	29Oct19/LIB105112_S4_L002_R1_001.fastq.gz	28Oct19/LIB105112_S4_L003_R1_001.fastq.gz	29Oct19/LIB105112_S4_L003_R1_001.fastq.gz	28Oct19/LIB105112_S4_L004_R1_001.fastq.gz	29Oct19/LIB105112_S4_L004_R1_001.fastq.gz > combined.fastq/LIB105112_S4_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105113_S5_L001_R1_001.fastq.gz	29Oct19/LIB105113_S5_L001_R1_001.fastq.gz	28Oct19/LIB105113_S5_L002_R1_001.fastq.gz	29Oct19/LIB105113_S5_L002_R1_001.fastq.gz	28Oct19/LIB105113_S5_L003_R1_001.fastq.gz	29Oct19/LIB105113_S5_L003_R1_001.fastq.gz	28Oct19/LIB105113_S5_L004_R1_001.fastq.gz	29Oct19/LIB105113_S5_L004_R1_001.fastq.gz > combined.fastq/LIB105113_S5_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105147_S6_L001_R1_001.fastq.gz	29Oct19/LIB105147_S6_L001_R1_001.fastq.gz	28Oct19/LIB105147_S6_L002_R1_001.fastq.gz	29Oct19/LIB105147_S6_L002_R1_001.fastq.gz	28Oct19/LIB105147_S6_L003_R1_001.fastq.gz	29Oct19/LIB105147_S6_L003_R1_001.fastq.gz	28Oct19/LIB105147_S6_L004_R1_001.fastq.gz	29Oct19/LIB105147_S6_L004_R1_001.fastq.gz > combined.fastq/LIB105147_S6_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105148_S7_L001_R1_001.fastq.gz	29Oct19/LIB105148_S7_L001_R1_001.fastq.gz	28Oct19/LIB105148_S7_L002_R1_001.fastq.gz	29Oct19/LIB105148_S7_L002_R1_001.fastq.gz	28Oct19/LIB105148_S7_L003_R1_001.fastq.gz	29Oct19/LIB105148_S7_L003_R1_001.fastq.gz	28Oct19/LIB105148_S7_L004_R1_001.fastq.gz	29Oct19/LIB105148_S7_L004_R1_001.fastq.gz > combined.fastq/LIB105148_S7_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105149_S8_L001_R1_001.fastq.gz	29Oct19/LIB105149_S8_L001_R1_001.fastq.gz	28Oct19/LIB105149_S8_L002_R1_001.fastq.gz	29Oct19/LIB105149_S8_L002_R1_001.fastq.gz	28Oct19/LIB105149_S8_L003_R1_001.fastq.gz	29Oct19/LIB105149_S8_L003_R1_001.fastq.gz	28Oct19/LIB105149_S8_L004_R1_001.fastq.gz	29Oct19/LIB105149_S8_L004_R1_001.fastq.gz > combined.fastq/LIB105149_S8_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105150_S9_L001_R1_001.fastq.gz	29Oct19/LIB105150_S9_L001_R1_001.fastq.gz	28Oct19/LIB105150_S9_L002_R1_001.fastq.gz	29Oct19/LIB105150_S9_L002_R1_001.fastq.gz	28Oct19/LIB105150_S9_L003_R1_001.fastq.gz	29Oct19/LIB105150_S9_L003_R1_001.fastq.gz	28Oct19/LIB105150_S9_L004_R1_001.fastq.gz	29Oct19/LIB105150_S9_L004_R1_001.fastq.gz > combined.fastq/LIB105150_S9_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105151_S10_L001_R1_001.fastq.gz	29Oct19/LIB105151_S10_L001_R1_001.fastq.gz	28Oct19/LIB105151_S10_L002_R1_001.fastq.gz	29Oct19/LIB105151_S10_L002_R1_001.fastq.gz	28Oct19/LIB105151_S10_L003_R1_001.fastq.gz	29Oct19/LIB105151_S10_L003_R1_001.fastq.gz	28Oct19/LIB105151_S10_L004_R1_001.fastq.gz	29Oct19/LIB105151_S10_L004_R1_001.fastq.gz > combined.fastq/LIB105151_S10_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105152_S11_L001_R1_001.fastq.gz	29Oct19/LIB105152_S11_L001_R1_001.fastq.gz	28Oct19/LIB105152_S11_L002_R1_001.fastq.gz	29Oct19/LIB105152_S11_L002_R1_001.fastq.gz	28Oct19/LIB105152_S11_L003_R1_001.fastq.gz	29Oct19/LIB105152_S11_L003_R1_001.fastq.gz	28Oct19/LIB105152_S11_L004_R1_001.fastq.gz	29Oct19/LIB105152_S11_L004_R1_001.fastq.gz > combined.fastq/LIB105152_S11_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105153_S12_L001_R1_001.fastq.gz	29Oct19/LIB105153_S12_L001_R1_001.fastq.gz	28Oct19/LIB105153_S12_L002_R1_001.fastq.gz	29Oct19/LIB105153_S12_L002_R1_001.fastq.gz	28Oct19/LIB105153_S12_L003_R1_001.fastq.gz	29Oct19/LIB105153_S12_L003_R1_001.fastq.gz	28Oct19/LIB105153_S12_L004_R1_001.fastq.gz	29Oct19/LIB105153_S12_L004_R1_001.fastq.gz > combined.fastq/LIB105153_S12_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105154_S13_L001_R1_001.fastq.gz	29Oct19/LIB105154_S13_L001_R1_001.fastq.gz	28Oct19/LIB105154_S13_L002_R1_001.fastq.gz	29Oct19/LIB105154_S13_L002_R1_001.fastq.gz	28Oct19/LIB105154_S13_L003_R1_001.fastq.gz	29Oct19/LIB105154_S13_L003_R1_001.fastq.gz	28Oct19/LIB105154_S13_L004_R1_001.fastq.gz	29Oct19/LIB105154_S13_L004_R1_001.fastq.gz > combined.fastq/LIB105154_S13_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105155_S14_L001_R1_001.fastq.gz	29Oct19/LIB105155_S14_L001_R1_001.fastq.gz	28Oct19/LIB105155_S14_L002_R1_001.fastq.gz	29Oct19/LIB105155_S14_L002_R1_001.fastq.gz	28Oct19/LIB105155_S14_L003_R1_001.fastq.gz	29Oct19/LIB105155_S14_L003_R1_001.fastq.gz	28Oct19/LIB105155_S14_L004_R1_001.fastq.gz	29Oct19/LIB105155_S14_L004_R1_001.fastq.gz > combined.fastq/LIB105155_S14_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105156_S15_L001_R1_001.fastq.gz	29Oct19/LIB105156_S15_L001_R1_001.fastq.gz	28Oct19/LIB105156_S15_L002_R1_001.fastq.gz	29Oct19/LIB105156_S15_L002_R1_001.fastq.gz	28Oct19/LIB105156_S15_L003_R1_001.fastq.gz	29Oct19/LIB105156_S15_L003_R1_001.fastq.gz	28Oct19/LIB105156_S15_L004_R1_001.fastq.gz	29Oct19/LIB105156_S15_L004_R1_001.fastq.gz > combined.fastq/LIB105156_S15_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105157_S16_L001_R1_001.fastq.gz	29Oct19/LIB105157_S16_L001_R1_001.fastq.gz	28Oct19/LIB105157_S16_L002_R1_001.fastq.gz	29Oct19/LIB105157_S16_L002_R1_001.fastq.gz	28Oct19/LIB105157_S16_L003_R1_001.fastq.gz	29Oct19/LIB105157_S16_L003_R1_001.fastq.gz	28Oct19/LIB105157_S16_L004_R1_001.fastq.gz	29Oct19/LIB105157_S16_L004_R1_001.fastq.gz > combined.fastq/LIB105157_S16_L001_R1_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105158_S17_L001_R1_001.fastq.gz	29Oct19/LIB105158_S17_L001_R1_001.fastq.gz	28Oct19/LIB105158_S17_L002_R1_001.fastq.gz	29Oct19/LIB105158_S17_L002_R1_001.fastq.gz	28Oct19/LIB105158_S17_L003_R1_001.fastq.gz	29Oct19/LIB105158_S17_L003_R1_001.fastq.gz	28Oct19/LIB105158_S17_L004_R1_001.fastq.gz	29Oct19/LIB105158_S17_L004_R1_001.fastq.gz > combined.fastq/LIB105158_S17_L001_R1_001.fastq.gz
```
``` r
mkdir -p combined.fastq; cat 28Oct19/LIB105109_S1_L001_R2_001.fastq.gz	29Oct19/LIB105109_S1_L001_R2_001.fastq.gz	28Oct19/LIB105109_S1_L002_R2_001.fastq.gz	29Oct19/LIB105109_S1_L002_R2_001.fastq.gz	28Oct19/LIB105109_S1_L003_R2_001.fastq.gz	29Oct19/LIB105109_S1_L003_R2_001.fastq.gz	28Oct19/LIB105109_S1_L004_R2_001.fastq.gz	29Oct19/LIB105109_S1_L004_R2_001.fastq.gz > combined.fastq/LIB105109_S1_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105110_S2_L001_R2_001.fastq.gz	29Oct19/LIB105110_S2_L001_R2_001.fastq.gz	28Oct19/LIB105110_S2_L002_R2_001.fastq.gz	29Oct19/LIB105110_S2_L002_R2_001.fastq.gz	28Oct19/LIB105110_S2_L003_R2_001.fastq.gz	29Oct19/LIB105110_S2_L003_R2_001.fastq.gz	28Oct19/LIB105110_S2_L004_R2_001.fastq.gz	29Oct19/LIB105110_S2_L004_R2_001.fastq.gz > combined.fastq/LIB105110_S2_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105111_S3_L001_R2_001.fastq.gz	29Oct19/LIB105111_S3_L001_R2_001.fastq.gz	28Oct19/LIB105111_S3_L002_R2_001.fastq.gz	29Oct19/LIB105111_S3_L002_R2_001.fastq.gz	28Oct19/LIB105111_S3_L003_R2_001.fastq.gz	29Oct19/LIB105111_S3_L003_R2_001.fastq.gz	28Oct19/LIB105111_S3_L004_R2_001.fastq.gz	29Oct19/LIB105111_S3_L004_R2_001.fastq.gz > combined.fastq/LIB105111_S3_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105112_S4_L001_R2_001.fastq.gz	29Oct19/LIB105112_S4_L001_R2_001.fastq.gz	28Oct19/LIB105112_S4_L002_R2_001.fastq.gz	29Oct19/LIB105112_S4_L002_R2_001.fastq.gz	28Oct19/LIB105112_S4_L003_R2_001.fastq.gz	29Oct19/LIB105112_S4_L003_R2_001.fastq.gz	28Oct19/LIB105112_S4_L004_R2_001.fastq.gz	29Oct19/LIB105112_S4_L004_R2_001.fastq.gz > combined.fastq/LIB105112_S4_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105113_S5_L001_R2_001.fastq.gz	29Oct19/LIB105113_S5_L001_R2_001.fastq.gz	28Oct19/LIB105113_S5_L002_R2_001.fastq.gz	29Oct19/LIB105113_S5_L002_R2_001.fastq.gz	28Oct19/LIB105113_S5_L003_R2_001.fastq.gz	29Oct19/LIB105113_S5_L003_R2_001.fastq.gz	28Oct19/LIB105113_S5_L004_R2_001.fastq.gz	29Oct19/LIB105113_S5_L004_R2_001.fastq.gz > combined.fastq/LIB105113_S5_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105147_S6_L001_R2_001.fastq.gz	29Oct19/LIB105147_S6_L001_R2_001.fastq.gz	28Oct19/LIB105147_S6_L002_R2_001.fastq.gz	29Oct19/LIB105147_S6_L002_R2_001.fastq.gz	28Oct19/LIB105147_S6_L003_R2_001.fastq.gz	29Oct19/LIB105147_S6_L003_R2_001.fastq.gz	28Oct19/LIB105147_S6_L004_R2_001.fastq.gz	29Oct19/LIB105147_S6_L004_R2_001.fastq.gz > combined.fastq/LIB105147_S6_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105148_S7_L001_R2_001.fastq.gz	29Oct19/LIB105148_S7_L001_R2_001.fastq.gz	28Oct19/LIB105148_S7_L002_R2_001.fastq.gz	29Oct19/LIB105148_S7_L002_R2_001.fastq.gz	28Oct19/LIB105148_S7_L003_R2_001.fastq.gz	29Oct19/LIB105148_S7_L003_R2_001.fastq.gz	28Oct19/LIB105148_S7_L004_R2_001.fastq.gz	29Oct19/LIB105148_S7_L004_R2_001.fastq.gz > combined.fastq/LIB105148_S7_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105149_S8_L001_R2_001.fastq.gz	29Oct19/LIB105149_S8_L001_R2_001.fastq.gz	28Oct19/LIB105149_S8_L002_R2_001.fastq.gz	29Oct19/LIB105149_S8_L002_R2_001.fastq.gz	28Oct19/LIB105149_S8_L003_R2_001.fastq.gz	29Oct19/LIB105149_S8_L003_R2_001.fastq.gz	28Oct19/LIB105149_S8_L004_R2_001.fastq.gz	29Oct19/LIB105149_S8_L004_R2_001.fastq.gz > combined.fastq/LIB105149_S8_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105150_S9_L001_R2_001.fastq.gz	29Oct19/LIB105150_S9_L001_R2_001.fastq.gz	28Oct19/LIB105150_S9_L002_R2_001.fastq.gz	29Oct19/LIB105150_S9_L002_R2_001.fastq.gz	28Oct19/LIB105150_S9_L003_R2_001.fastq.gz	29Oct19/LIB105150_S9_L003_R2_001.fastq.gz	28Oct19/LIB105150_S9_L004_R2_001.fastq.gz	29Oct19/LIB105150_S9_L004_R2_001.fastq.gz > combined.fastq/LIB105150_S9_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105151_S10_L001_R2_001.fastq.gz	29Oct19/LIB105151_S10_L001_R2_001.fastq.gz	28Oct19/LIB105151_S10_L002_R2_001.fastq.gz	29Oct19/LIB105151_S10_L002_R2_001.fastq.gz	28Oct19/LIB105151_S10_L003_R2_001.fastq.gz	29Oct19/LIB105151_S10_L003_R2_001.fastq.gz	28Oct19/LIB105151_S10_L004_R2_001.fastq.gz	29Oct19/LIB105151_S10_L004_R2_001.fastq.gz > combined.fastq/LIB105151_S10_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105152_S11_L001_R2_001.fastq.gz	29Oct19/LIB105152_S11_L001_R2_001.fastq.gz	28Oct19/LIB105152_S11_L002_R2_001.fastq.gz	29Oct19/LIB105152_S11_L002_R2_001.fastq.gz	28Oct19/LIB105152_S11_L003_R2_001.fastq.gz	29Oct19/LIB105152_S11_L003_R2_001.fastq.gz	28Oct19/LIB105152_S11_L004_R2_001.fastq.gz	29Oct19/LIB105152_S11_L004_R2_001.fastq.gz > combined.fastq/LIB105152_S11_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105153_S12_L001_R2_001.fastq.gz	29Oct19/LIB105153_S12_L001_R2_001.fastq.gz	28Oct19/LIB105153_S12_L002_R2_001.fastq.gz	29Oct19/LIB105153_S12_L002_R2_001.fastq.gz	28Oct19/LIB105153_S12_L003_R2_001.fastq.gz	29Oct19/LIB105153_S12_L003_R2_001.fastq.gz	28Oct19/LIB105153_S12_L004_R2_001.fastq.gz	29Oct19/LIB105153_S12_L004_R2_001.fastq.gz > combined.fastq/LIB105153_S12_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105154_S13_L001_R2_001.fastq.gz	29Oct19/LIB105154_S13_L001_R2_001.fastq.gz	28Oct19/LIB105154_S13_L002_R2_001.fastq.gz	29Oct19/LIB105154_S13_L002_R2_001.fastq.gz	28Oct19/LIB105154_S13_L003_R2_001.fastq.gz	29Oct19/LIB105154_S13_L003_R2_001.fastq.gz	28Oct19/LIB105154_S13_L004_R2_001.fastq.gz	29Oct19/LIB105154_S13_L004_R2_001.fastq.gz > combined.fastq/LIB105154_S13_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105155_S14_L001_R2_001.fastq.gz	29Oct19/LIB105155_S14_L001_R2_001.fastq.gz	28Oct19/LIB105155_S14_L002_R2_001.fastq.gz	29Oct19/LIB105155_S14_L002_R2_001.fastq.gz	28Oct19/LIB105155_S14_L003_R2_001.fastq.gz	29Oct19/LIB105155_S14_L003_R2_001.fastq.gz	28Oct19/LIB105155_S14_L004_R2_001.fastq.gz	29Oct19/LIB105155_S14_L004_R2_001.fastq.gz > combined.fastq/LIB105155_S14_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105156_S15_L001_R2_001.fastq.gz	29Oct19/LIB105156_S15_L001_R2_001.fastq.gz	28Oct19/LIB105156_S15_L002_R2_001.fastq.gz	29Oct19/LIB105156_S15_L002_R2_001.fastq.gz	28Oct19/LIB105156_S15_L003_R2_001.fastq.gz	29Oct19/LIB105156_S15_L003_R2_001.fastq.gz	28Oct19/LIB105156_S15_L004_R2_001.fastq.gz	29Oct19/LIB105156_S15_L004_R2_001.fastq.gz > combined.fastq/LIB105156_S15_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105157_S16_L001_R2_001.fastq.gz	29Oct19/LIB105157_S16_L001_R2_001.fastq.gz	28Oct19/LIB105157_S16_L002_R2_001.fastq.gz	29Oct19/LIB105157_S16_L002_R2_001.fastq.gz	28Oct19/LIB105157_S16_L003_R2_001.fastq.gz	29Oct19/LIB105157_S16_L003_R2_001.fastq.gz	28Oct19/LIB105157_S16_L004_R2_001.fastq.gz	29Oct19/LIB105157_S16_L004_R2_001.fastq.gz > combined.fastq/LIB105157_S16_L001_R2_001.fastq.gz
mkdir -p combined.fastq; cat 28Oct19/LIB105158_S17_L001_R2_001.fastq.gz	29Oct19/LIB105158_S17_L001_R2_001.fastq.gz	28Oct19/LIB105158_S17_L002_R2_001.fastq.gz	29Oct19/LIB105158_S17_L002_R2_001.fastq.gz	28Oct19/LIB105158_S17_L003_R2_001.fastq.gz	29Oct19/LIB105158_S17_L003_R2_001.fastq.gz	28Oct19/LIB105158_S17_L004_R2_001.fastq.gz	29Oct19/LIB105158_S17_L004_R2_001.fastq.gz > combined.fastq/LIB105158_S17_L001_R2_001.fastq.gz
```

### 1.0 fastqc
``` r
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=1:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=4   # 36 processor core(s) per node
#SBATCH --mem=16G   # maximum memory per node
#SBATCH --job-name="fastqc"
#SBATCH --array=1-34
#SBATCH --output="logs/out.fastqc.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.fastqc.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

i=$(($SLURM_ARRAY_TASK_ID - 1))

module load fastqc/0.11.7-d5mgqc7

fastq=(`ls data/combined.fastq/*.gz`)

out=results/1.0.fastqc.out/
mkdir -p $out

fastqc -o $out -t 4 --extract ${fastq[$i]}
```

### 1.1 trimmomatic
``` r
#!/bin/bash


# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 36 processor core(s) per node
#SBATCH --mem=32G   # maximum memory per node
#SBATCH --job-name="trimming"
#SBATCH --array=1-17
#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


i=$(($SLURM_ARRAY_TASK_ID - 1))

module load java/1.8.0_171
module load trimmomatic/0.36-lkktrba

in_dir=data/combined.fastq

R1=(`ls ${in_dir}/*_R1_001.fastq.gz`)
R2=(`ls ${in_dir}/*_R2_001.fastq.gz`)

name=(`ls ${in_dir}/*_R1_001.fastq.gz | perl -p -e   's{.+/(.+?)_R1_001.fastq.gz}{$1}g' `)

trimmomatic  PE -threads 8 \
          -phred33 ${R1[$i]} ${R2[$i]} \
          ${name[$i]}.R1Ptrim.fastq.gz  ${name[$i]}.R1Utrim.fastq.gz \
          ${name[$i]}.R2Ptrim.fastq.gz  ${name[$i]}.R2Utrim.fastq.gz  \
          ILLUMINACLIP:adapter.fa:2:30:10:1:true LEADING:0 TRAILING:3 \
          SLIDINGWINDOW:4:20 LEADING:0 TRAILING:3 MINLEN:25
```

### 1.2 trimPolyATG
``` r
#!/bin/bash


# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=8:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # 36 processor core(s) per node
#SBATCH --mem=2G   # maximum memory per node
#SBATCH --job-name="trimming"
#SBATCH --array=1-68
#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


i=$(($SLURM_ARRAY_TASK_ID - 1))


in_dir=data/trimmed.fastq
out_dir=data/trimPolyATG.data
mkdir -p ${out_dir}


R1=(`ls ${in_dir}/*.fastq.gz`)

name=(`ls ${in_dir}/*.fastq.gz | perl -p -e   's{.+/(.+?)trim.fastq.gz}{$1}g' `)

zcat ${R1[$i]} | perl scripts/trimPolyATG.pl | gzip   > ${out_dir}/${name[$i]}.fastq.gz
```

### trimPolyATG.pl
``` r
#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $trim_polyG = 1;
my $trim_polyA = 1;
my $trim_polyT = 1;

my ($windowSize, $stepSize, $min_content, $min_length) = (10, 1, 90, 25);
sub get_base_content;
sub get_kmer_base_content;
sub get_trimming_site;
sub trim_read;

while (my $header =  <STDIN>)
{
   chomp $header;
   chomp (my $seq = <STDIN>);
   chomp (my $conj = <STDIN>);
   chomp (my $qual = <STDIN>);

   ## trimming polyG from the 3' end, for Nextseq data
   if ($seq =~/G{10,}/ && $trim_polyG && length($seq) >= 25)
    {
        ($seq, $qual) = @{trim_read($seq, $qual, "G")};
    }

    ### trimming polyT at the 5' end only if the remaining length is >=25
    if (length($seq) >= 25)
    {
        if ($seq =~/T{10,}/ && $trim_polyT)
        {
            ## reverse sequence and quality
            $seq = reverse $seq;
            $qual = reverse $qual;

           ($seq, $qual) = @{trim_read($seq, $qual, "T")};

            ## reverse sequence and quality back to the read directions
            $seq = reverse $seq;
            $qual = reverse $qual;
        }
    }

    ### trimming polyA at the 3' end only if the remaining length is >=25
    ### even after trimming polyG
    if (length($seq) >= 25)
    {
        if ($seq =~/A{10,}/ && $trim_polyA)
        {
            ($seq, $qual) = @{trim_read($seq, $qual, "A")};
        }
    }

    ### output reads if the current length is >=25
    if (length($seq) >= 25)
    {
        my @fastq = ($header, $seq, $conj, $qual);
        print join "\n", @fastq;
        print "\n";
    }
}

sub get_base_content
{
    my ($kmer, $bases) = @_;

    no strict;
    no warnings;
    my @nt_count = eval "\$kmer =~ tr/$bases//";
    my @nt_pct = map {$_ / length($kmer)} @nt_count;
    my %freq = ("$bases" => $nt_pct[0]);
    return \%freq;
}


sub get_kmer_base_contents
{
    my ($sequence, $windowSize, $stepSize, $base) = @_;
    my @kmers=();
    for( my $windowStart = 0 ; $windowStart <= (length($sequence) - $windowSize);
        $windowStart += $stepSize)
    {
        my $kmer = substr($sequence, $windowStart, $windowSize);
        my $kemr_nt_content = get_base_content($kmer, $base); ## reference to a hash
        push @kmers, {$kmer => $kemr_nt_content};  ## array of anonymous hash of anonymous hash
    }
    return \@kmers;
}



sub get_trimming_site
{
    ## $min_content: min G contents,considered as polyG stretch
    ## $base to trim
    my ($kmer_base_content_ref, $base, $min_content, $read_length, $stepSize) =  @_;
    my $highest_content = $min_content;
    my @base_pct = ();
    my @kmers = ();


    my @kmer_base_content_tmp = @$kmer_base_content_ref; ## dereference the first level to get
                                                         ## array of anonymous hash of anonymous hash


    ## add kmer sequence and content of base to be trimmed to two parallel arrays
    for my $kmer (@kmer_base_content_tmp)
    {
       push @kmers, keys(%{$kmer});

       # print Dumper(%{(values (%{$kmer}))[0]});  ### values(%{$kmer}) return an array of a
                                                   ### single element, which is an anonymous hash
       push @base_pct, ${(values (%{$kmer}))[0]}{$base};
    }


    my $trimming_start = $read_length;

    for (my $i=0; $i < scalar(@base_pct); $i++)
    {
        if ($kmers[$i] =~/^$base/ && $base_pct[$i]*100 >= $min_content)
        {
            ## scan the following kemrs to make sure the minumal baseq% is satisfied
            my $all_high = 1;
            for (my $j = $i +1; $j < scalar(@base_pct); $j++)
            {
                if ($base_pct[$j]*100 < $min_content)
                {
                    $all_high = 0;
                    last;
                }
            }
            if ($all_high)
            {
                $trimming_start = $i * $stepSize; ## trimming site
                ## print $trimming_start."\n";
                return $trimming_start;
            }
        }
    }
    return $trimming_start;
}

sub trim_read
{
    my ($seq, $qual, $trim_base) = @_;
    my $read_length = length($seq);
    my $kmer_base_content_ref = get_kmer_base_contents($seq, $windowSize, $stepSize, $trim_base);
    my $triming_start = get_trimming_site($kmer_base_content_ref, $trim_base,
                                          $min_content, $read_length, $stepSize);
    $seq  = substr($seq, 0, $triming_start);
    $qual = substr($qual, 0, $triming_start);

    ## return anomymous array containing $seq and quality after trimming
    return [$seq, $qual];
}
```

### 1.3 repair
``` r
#!/bin/sh

# Copy/paste this job script into a text file and submit with the command:
# sbatch thefilename
#SBATCH --time=12:00:00 # walltime limit (HH:MM:SS)
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks-per-node=1 # 36 processor core(s) per node
#SBATCH --mem=32G # maximum memory per node
#SBATCH --job-name="repair"
#SBATCH --array=7,14
#SBATCH --output="logs/out.rcorrect.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.rcorrect.%A_%a.txt" # job standard error file (%j replaced by job id)
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load jdk/8u141-b15-xf726oe
module load bbmap


i=$(($SLURM_ARRAY_TASK_ID - 1))

in_dir=data/trimPolyATG.data
out=data/repaired.fastq

R1=(`ls  $in_dir/*R1P.fastq.gz`)
R2=(`ls  $in_dir/*R2P.fastq.gz`)
name=(`ls $in_dir/*R1P.fastq.gz | perl -p -e 's{.+/(.+?).R1P.fastq.gz}{$1}'`)

mkdir -p ${out}

## using BBmap to re-pair the trimmed reads
repair.sh in1=${R1[$i]} in2=${R2[$i]} \
          out1=${out}/${name[$i]}.PE.R1.fastq.gz   out2=${out}/${name[$i]}.PE.R2.fastq.gz \
          outs=${out}/${name[$i]}.SE.fastq.gz  \
          repair overwrite  ziplevel=9  -Xmx16G
```

### 1.4 reversecomplement
``` r
#!/bin/sh

# Copy/paste this job script into a text file and submit with the command:
# sbatch thefilename
#SBATCH --time=6:00:00 # walltime limit (HH:MM:SS)
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks-per-node=1 # 36 processor core(s) per node
#SBATCH --mem=2G # maximum memory per node
#SBATCH --job-name="repair"
#SBATCH --array=1-17
#SBATCH --output="logs/out.rcorrect.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.rcorrect.%A_%a.txt" # job standard error file (%j replaced by job id)
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load jdk/8u141-b15-xf726oe
module load bbmap


i=$(($SLURM_ARRAY_TASK_ID - 1))

uread1=(`ls data/trimPolyATG.data/*R1U.fastq.gz`)
uread2=(`ls data/trimPolyATG.data/*R2U.fastq.gz`)

name=(`ls data/trimPolyATG.data/*R2U.fastq.gz | perl -p -e 's{.+/(.+?).R2U.fastq.gz}{$1}'`)
out=data/SE.fastq

mkdir -p $out

cp ${uread1[$i]} $out/${name[$i]}.SE.fastq.gz

seqtk seq -r ${uread2[$i]} | perl -p -e 's/2(:N:0:[ATCGN]+)/1$1/ if ($. % 4 ==1)' | gzip >> $out/${name[$i]}.SE.fastq.gz

in=data/repaired.fastq

seread=(`ls $in/*.SE.fastq.gz `)

zcat ${seread[$i]} | paste - - - - | awk 'BEGIN{FS=OFS="\t"} /2:N:0:[ATCGN]+/ {print $1"\n"$2"\n"$3"\n"$4}' | seqtk seq -r - | \
                    perl -p -e 's/2(:N:0:[ATCGN]+)/1$1/ if ($. % 4 ==1)'| gzip  >> $out/${name[$i]}.SE.fastq.gz

zcat ${seread[$i]} | paste - - - - | awk 'BEGIN{FS=OFS="\t"} /1:N:0:[ATCGN]+/ {print $1"\n"$2"\n"$3"\n"$4}' | gzip  >> $out/${name[$i]}.SE.fastq.gz
```

### 1.5 PE.Star
``` r
#!/bin/bash

##copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 36 processor core(s) per node
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --job-name="STAR"
#SBATCH --array=1-17
#SBATCH --output="logs/out.star.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.star.%A_%a.txt" # job standard error file (%j replaced by job id)
##SBATCH --dependency=afterok:2655991

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load star/2.5.3a-rdzqclb

i=$(($SLURM_ARRAY_TASK_ID - 1))

in=data/repaired.fastq

R1=(`ls $in/*PE.R1.fastq.gz`)
R2=(`ls $in/*PE.R2.fastq.gz`)
names=(`ls $in/*PE.R2.fastq.gz | perl -p -e s'{.+/(.+)\.R2.fastq.gz}{$1}g' `)
out=results/STAR.out

mkdir -p $out

genomeDir=~/juber/pig.genome/SSC11.1.overhang.76
gtf=~/juber/pig.genome/Sus_scrofa.Sscrofa11.1.97.gtf

     STAR  --runThreadN 16  \
      --readFilesCommand  zcat  \
      --outFileNamePrefix  ${out}/${names[${i}]}. \
      --genomeDir  $genomeDir  \
      --readFilesIn  ${R1[${i}]} ${R2[${i}]} \
      --sjdbGTFfile  $gtf \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20  \
      --alignSJoverhangMin 8  \
      --outFilterIntronMotifs RemoveNoncanonical \
      --alignSJDBoverhangMin 1  \
      --outFilterMismatchNmax 30  \
      --seedSearchStartLmax 50 \
      --seedSearchStartLmaxOverLread 0.5 \
      --alignIntronMin 20  \
      --alignIntronMax 1000000  \
      --alignMatesGapMax 1000000  \
      --outSAMstrandField intronMotif \
      --outSAMtype BAM Unsorted  \
      --outReadsUnmapped Fastx
```



### 1.6 SE.Star
``` r
#!/bin/bash

##copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 36 processor core(s) per node
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --job-name="STAR"
#SBATCH --array=1-17
#SBATCH --output="logs/out.star.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.star.%A_%a.txt" # job standard error file (%j replaced by job id)
##SBATCH --dependency=afterok:2655991

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load star/2.5.3a-rdzqclb

i=$(($SLURM_ARRAY_TASK_ID - 1))

in=data/SE.fastq

R1=(`ls $in/*fastq.gz`)
names=(`ls $in/*fastq.gz | perl -p -e s'{.+/(.+)\.fastq.gz}{$1}g' `)
out=results/STAR.out

mkdir -p $out

genomeDir=~/juber/pig.genome/SSC11.1.overhang.76
gtf=~/juber/pig.genome/Sus_scrofa.Sscrofa11.1.97.gtf

     STAR  --runThreadN 16  \
      --readFilesCommand  zcat  \
      --outFileNamePrefix  ${out}/${names[${i}]}. \
      --genomeDir  $genomeDir  \
      --readFilesIn  ${R1[${i}]} \
      --sjdbGTFfile  $gtf \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20  \
      --alignSJoverhangMin 8  \
      --outFilterIntronMotifs RemoveNoncanonical \
      --alignSJDBoverhangMin 1  \
      --outFilterMismatchNmax 30  \
      --seedSearchStartLmax 50 \
      --seedSearchStartLmaxOverLread 0.5 \
      --alignIntronMin 20  \
      --alignIntronMax 1000000  \
      --alignMatesGapMax 1000000  \
      --outSAMstrandField intronMotif \
      --outSAMtype BAM Unsorted  \
      --outReadsUnmapped Fastx
```

### 1.7 featurecount
``` r
#!/bin/bash

##copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 36 processor core(s) per node
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --job-name="STAR"
#SBATCH --array=1
#SBATCH --output="logs/out.star.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.star.%A_%a.txt" # job standard error file (%j replaced by job id)
#SBATCH --dependency=afterok:2656183,2656200

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load subread/1.6.0-ak6vxhs

i=$(($SLURM_ARRAY_TASK_ID - 1))

export GTF=~/juber/pig.genome/Sus_scrofa.Sscrofa11.1.97.gtf
export genome=~/juber/pig.genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

PE_bam=(`ls results/STAR.out/*PE*.bam `)
SE_bam=(`ls results/STAR.out/*SE*.bam `)

PE_table=FAANG.sortedcell.PE.count.table
SE_table=FAANG.sortedcell.SE.count.table

#### paired end mapping
featureCounts  -p -a $GTF  -d 25 -s 2  -F GTF -g gene_id -G $genome   -o ${PE_table} -t exon -T 8 ${PE_bam[@]}

cut -f 2-6 --complement ${PE_table}  | perl -p -e 's/^#.+\n$//' > 1.0.formatted.${PE_table}


#### single end mapping
featureCounts -a $GTF  -d 25  -s 2  -F GTF -g gene_id -G $genome   -o ${SE_table} -t exon -T 8  ${SE_bam[@]}

cut -f 2-6 --complement ${SE_table}  | perl -p -e 's/^#.+\n$//' > 1.0.formatted.${SE_table}
```




### References
Herrera-Uribe, J., Liu, H., Byrne, K.A., Bond, Z.F., Loving, C.L., and Tuggle, C.K. (2020). Changes in H3K27ac at Gene Regulatory Regions in Porcine Alveolar Macrophages Following LPS or PolyIC Exposure. Front Genet 11, 817. doi: 10.3389/fgene.2020.00817.

Hicks, S.C., Okrah, K., Paulson, J.N., Quackenbush, J., Irizarry, R.A., and Bravo, H.C. (2018). Smooth quantile normalization. Biostatistics (Oxford, England) 19(2), 185-198. doi: 10.1093/biostatistics/kxx028.

Paulson, J.N., Chen, C.Y., Lopes-Ramos, C.M., Kuijjer, M.L., Platig, J., Sonawane, A.R., et al. (2017). Tissue-aware RNA-Seq processing and normalization for heterogeneous and sparse data. BMC Bioinformatics 18(1), 437. doi: 10.1186/s12859-017-1847-x.

Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15(12), 550. doi: 10.1186/s13059-014-0550-8.
