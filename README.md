# eutopsQC

Preprocessing pipeline for raw epigenetic data from raw IDAT (.idat) to a preprocessed beta matrix, including provision of a quality control (QC) report.

input: input folder direction to raw data

output: folder, where should data be stored?

report: folder, where should QC report be placed input?

array: EPIC, mouse or other

pheno: optional, point to .csv, .Rdata, .txt file. Important - needs to contain one column "basename" corresponding to idat names; Each column should not be uniform and contain only the same variable

path_to_bad_sample_list: optional, point to .csv list of flagged samples to be removed


## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("chiaraherzog/eutopsQC")
```