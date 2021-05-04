---
Title:  "Porcine immune cell atlas single cell RNA sequencing analysis pipeline"
Author: "Sathesh K. Sivasankaran"
Date:   "Feb 15, 2021"
---

# Porcine Immune Single Cell Atlas

## Pre-sequence processing
### sequence quality
The quality of the fastq reads were analyzed using FastQC (version v0.11.8) [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/]. Jobs were submitted in SciNet using below mentioned bash script which uses the parallel function.

```
module load fastqc
module load parallel
parallel -j 6 "fastqc {} -t 6 -o FQC/" ::: *.fastq.gz
```

Results were stored in one location and combined using MultiQC (https://multiqc.info/) using the below mentioned command.
```
module load python_2
multiqc .
```
MultiQC results could be found.\
![A](Notebook/FastQC/A.html)\
![B](Notebook/FastQC/B.html)\
![C](Notebook/FastQC/C.html)\
![CT21NOV](Notebook/FastQC/CT21NOV.html)\
![CT230OCT](Notebook/FastQC/CT230OCT.html)\
![PBMC1](Notebook/FastQC/PBMC1.html)\
![PBMC2](Notebook/FastQC/PBMC2.html)

### R-correction
Reads 2 (R2) were corrected for errors using Rcorrector (https://github.com/mourisl/Rcorrector).
```
module load java/1.8.0_121
module load zlib/1.2.7
libs = CT2-1NOV
R2=(`ls  -1 fastq/${libs[$i]}*_R2_001.fastq.gz |perl -p -e 's/\n/,/' |perl -p -e  's/,$//'` )
perl ~/tools/rcorrector/run_rcorrector.pl -s  ${R2[0]}  -k 25  -od ./fastq/Rcorrected  -t 72
```

### 3' AT trimming
3’ polyA tails longer than 10 bases were trimmed using custom Perl scripts. Perl script is located in Notebook/Scripts//
```
gzip -d CT2-1NOV-1_S5_L001_R2_001.cor.fq.gz
perl -p  -i -e  's/([^ ]+s+[^ ]+).+/$1/ if ($. % 4 == 1)' CT2-1NOV-1_S5_L001_R2_001.cor.fq
gzip CT2-1NOV-1_S5_L001_R2_001.cor.fq
zcat CT2-1NOV-1_S5_L001_R2_001.cor.fq.gz | perl trimming.AT.pl RNA-seq | gzip > CT2-1NOV-1_S5_L001_R2_001.fastq.gz
```

### Read correction
R-correction and AT trimming can drop some R2 reads due to its poor quality. To match R1 and R2 reads, R2 longer than 25 bases were re-paired and matched with R1 using BBMap (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/).
```
module load java/1.8.0_121
~/tools/bbmap/./repair.sh in1=CT2-1NOV-1_S5_L001_R1_001.fastq.gz in2=CT2-1NOV-1_S5_L001_R2_001.fastq.gz out1=Repaired/CT2-1NOV-1_S5_L001_R1_001.fastq.gz out2=Repaired/CT2-1NOV-1_S5_L001_R2_001.fastq.gz outs=Repaired/CT2-1NOV-1_S5_L001.SE.fastq.gz repair -da ziplevel=2 -Xmx20g
```

### Genome and annotation
Sus scrofa genome Sscrofa 11.1 (V97; http://ftp.ensembl.org/pub/release-97/fasta/sus_scrofa/dna/) and annotation GTF (v11.1.97; http://ftp.ensembl.org/pub/release-97/gtf/sus_scrofa/) from Ensemble were used to build the reference genome index (Yates et al., 2020).
The annotation GTF file was modified to include both gene symbol (if available) and Ensembl ID as a gene reference (e.g. GZMA_ENSSSCG00000016903) using custom Perl scripts.
```
perl -p -e 'if (/gene_name/) {s{(gene_id\s+"([^"]+).+?gene_name\s+")([^"]+)}{$1$3_$2}} elsif (!/^#/ && /gene_id/) {s/(gene_id\s+"([^"]+)";\s+)/$1gene_name "$2"; /}'  Sus_scrofa.Sscrofa11.1.97.gtf > Sus_scrofa.Sscrofa11.1.97.M.gtf

```
### Genome indexing and alignment
Pig genome index was created using CellRanger V4.0 (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count).
```
cellranger mkref --genome=ssc97 \
	--fasta=Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	--genes=Sus_scrofa.Sscrofa11.1.97.M.gtf \
	--nthreads=38 \
	--ref-version=11.1.97
  ```
Reads were mapped to reference genome and read counts were quantified (count matrix) using cellRanger V4.0.
```
cellranger count --id=CT21NOV \
                   --transcriptome=ssc97 \
                   --fastqs=ForRep/Repaired \
                   --sample=CT2-1NOV-1,CT2-1NOV-2,CT2-1NOV-3,CT2-1NOV-4 \
                   --localcores=38
```

Table: Mapping statistics.
| Samples                        | A           | B           | C          | CT2-01NOV   | CT2-30OCT   | PBMC1       | PBMC2    |
|--------------------------------|-------------|-------------|------------|-------------|-------------|-------------|----------|
| Estimated No. of Cells         | 4,581       | 4,204       | 5,064      | 4,040       | 4,638       | 8,283       | 7,544    |
| Mean per Cell                  | 26,105      | 28,414      | 16,686     | 40,038      | 32,025      | 19,730      | 26,748   |
| Median Genes per Cell          | 618         | 615         | 595        | 765         | 841         | 957         | 1,138    |
| No. of Reads                   | 119,590,111 | 119,453,412 | 84,500,871 | 161,757,334 | 148,532,329 | 163,426,564 | ######## |
| Valid Barcodes                 | 98.30%      | 98.20%      | 98.10%     | 98.30%      | 98.30%      | 98.20%      | 98.10%   |
| Sequencing Saturation          | 80.20%      | 82.80%      | 71.10%     | 80.30%      | 78.40%      | 61.70%      | 65.80%   |
| Q30 Bases in Barcode           | 97.60%      | 97.50%      | 97.50%     | 97.10%      | 97.20%      | 99.10%      | 99.10%   |
| Q30 Bases in RNA Read          | 80.50%      | 79.40%      | 80.60%     | 78.40%      | 84.50%      | 76.00%      | 76.60%   |
| Q30 Bases in UMI               | 97.80%      | 97.80%      | 97.70%     | 97.40%      | 97.50%      | 99.20%      | 99.20%   |
| ֎ to Genome                    | 94.00%      | 93.40%      | 93.50%     | 91.00%      | 93.30%      | 95.40%      | 95.20%   |
| ֎ Confidently to Genome        | 90.50%      | 90.30%      | 90.20%     | 87.80%      | 90.00%      | 92.60%      | 92.60%   |
| ֎ Confidently to Intergenic    | 8.30%       | 8.80%       | 9.00%      | 8.50%       | 8.30%       | 7.10%       | 7.60%    |
| ֎ Confidently to Intronic      | 28.70%      | 30.00%      | 33.10%     | 27.50%      | 26.90%      | 22.00%      | 24.10%   |
| ֎ Confidently to Exonic        | 53.40%      | 51.60%      | 48.10%     | 51.90%      | 54.80%      | 63.50%      | 60.90%   |
| ֎ Confidently to Transcriptome | 49.90%      | 47.80%      | 44.30%     | 48.00%      | 50.20%      | 59.90%      | 57.00%   |
| ֎ Antisense to Gene            | 1.20%       | 1.30%       | 1.40%      | 1.20%       | 1.30%       | 1.00%       | 1.20%    |
| Fraction in Cells              | 91.20%      | 85.60%      | 86.20%     | 79.30%      | 91.90%      | 92.30%      | 90.60%   |
| Total Genes Detected           | 13,478      | 13,424      | 13,468     | 13,607      | 13,707      | 14,217      | 14,425   |
| Median UMI Counts per Cell     | 1,555       | 1,374       | 1,263      | 1,953       | 2,098       | 3,172       | 3,418    |

֎ = Reads Mapped

### Alignment and read counting was successfully completed.
### ___END___
