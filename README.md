# eutopsQC

Preprocessing pipeline for raw epigenetic data from raw IDAT (.idat) or signal intensity files (e.g., obtained from GEO) to a preprocessed beta matrix, including a quality control (QC) report.

* `input`: input folder direction to raw data (idat) or .txt.gz file (signal intensity file)
* `output`: folder where beta matrix should be stored
* `report`: folder where the QC report should be generated
* `array`: EPIC, mouse or other
* `pheno`:
  - optional for .idat input; can point to .csv, .Rdata, .txt file. Important: REQUIRES column basename! each column should have multiple values (non-uniform) for QC
  - mandatory for signal intensity files
* `path_to_bad_sample_list`: optional, point to .csv list of flagged samples to be removed
* `cores`: number of core to use for ChAMP normalisation. 4 by default.
* `overwrite`: if a beta matrix is present in the output directory already, overwrite = F prevents accidental overwriting.
* `run.name`: provide run name
* `beta.subset.compatible`: for EPIC version 2 only: should a version1/version2 compatible matrix be generated using the same rownames and probes?

IDAT processing only:

* `by.dir`: process input by directory (instead of all at once) - FALSE by default, but recommended for larger projects (>10 beadchips; depending on available RAM)
* `save.rs`: extract SNP (rs) probe values. FALSE by default.
* find.files: only selects those files specified in basename from a given folder. REQUIRES a pheno file with basename column, run.name, and sets by.dir to F (not applicable)

Signal matrix processing only:

* `Uname`: column grep name for Unmethylated values
* `Mname`: column grep name for Methylated values
* `detPname`: column grep name for detection P values
* `sep`: separator for signal intensity file (\t by default)

## Installation

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("chiaraherzog/eutopsQC")
```
