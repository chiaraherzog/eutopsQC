# eutopsQC

Preprocessing pipeline for raw epigenetic data from raw IDAT (.idat) to a preprocessed beta matrix, including provision of a quality control (QC) report.

* input: input folder direction to raw data
* output: folder, where should data be stored?
* report: folder, where should QC report be placed input?
* array: EPIC, mouse or other
* pheno: optional, point to .csv, .Rdata, .txt file. REQUIRES column basename! each column should have multiple values (non-uniform), for QC
* path_to_bad_sample_list: optional, point to .csv list of flagged samples to be removed
* cores: number of core to use for ChAMP normalisation. 4 by default.
* by.dir: process input by directory (instead of all at once) - FALSE by default, but recommended for larger projects (>10 beadchips; depending on available RAM)
* save.rs: extract SNP (rs) probe values. FALSE by default.
* find.files: only selects those files specified in basename from a given folder. REQUIRES a pheno file with basename column, run.name, and sets by.dir to F (not applicable)

## Installation

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("chiaraherzog/eutopsQC")
```
