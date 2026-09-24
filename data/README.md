# Input data

The analysis uses the public STRAP bulk RNA-seq cohort from accession [E-MTAB-13733](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13733).

The sample metadata file is included as `E-MTAB-13733.sdrf.txt`. Before running the pipeline, place the expression-count object at:

```text
data/strap_counts.RData
```

The RData file is expected to contain the raw gene-by-sample count object accepted by `load_counts()` in `pipeline/functions/01_preprocessing.R`. Raw counts are deliberately excluded from Git because of file size and data-distribution considerations.
