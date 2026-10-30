# Public-source worked data

These small derived data files support examples 14 and 15. They contain previously
curated expression for the study's frozen 60-gene panels, public sample identifiers,
histology labels and overall-survival endpoints. The source data are public TCGA/GDC and
GEO records, and the original investigators are cited when they are reused.

| File | Development cohorts | Development n | Worked external cohort | External n |
|---|---|---:|---|---:|
| lung_worked_data.rds | TCGA-LUAD (503), TCGA-LUSC (493) | 996 | GSE3141 | 110 |
| glioma_worked_data.rds | TCGA-LGG (511), TCGA-GBM (153) | 664 | GSE43378 | 50 |

Each RDS contains `train`, `test`, `train_cohorts`, `test_cohort`, and `provenance`.
Each split contains `x` (samples x 60 genes), `domain` (histology, not assay batch),
and `survival` (`time` in months; `event` 1=death, 0=censored). The training and test
gene order is identical; sample IDs are matrix and endpoint-table row names.

## Sources and derivation

- [TCGA program](https://www.cancer.gov/ccg/research/genome-sequencing/tcga) and
  [Genomic Data Commons](https://portal.gdc.cancer.gov/): the previously curated
  eligible primary-tumor patients and public clinical endpoints.
- [UCSC Xena GDC hub](https://gdc.xenahubs.net/): TCGA-LUAD, TCGA-LUSC, TCGA-LGG and
  TCGA-GBM `star_tpm` expression, on the supplied `log2(TPM+1)` scale. Public source
  downloads follow `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz`
  (replace the cohort name for the other three). Ensembl versions were removed,
  protein-coding symbols mapped using GENCODE v36, and duplicate entries/symbols
  combined by median. One eligible primary expression sample per patient was selected
  in the original curation.
- [GSE3141](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3141): GPL570
  primary lung tumors, processed MAS5 signals. The existing workflow applied log2
  to positive deposited signals, mapped probes consistently and aggregated by gene.
  Of 111 arrays, GSM70230 lacked survival time, leaving 110 with usable OS.
- [GSE43378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43378): processed
  public glioma microarray expression and public survival annotations, using the
  existing gene mapping and eligibility curation. The export retains all 50 eligible
  samples with the curated expression as-is.

The frozen clinical files express OS in months. Original day-based endpoints were
converted with 365.25/12 days per month in the existing curation. Both data files are
subsets of the existing processed matrices, not raw assay releases. No extra gene-wise
z-scoring, pair encoding or domain adjustment has been applied to these `x` matrices.
The functions in the examples perform those steps after loading the data.

`source_manifest.csv` records source-file checksums and exported data-file checksums.
`feature_panels.csv` lists the frozen genes in order. These trace this package export.

## Intended use

The real examples use fixed lambda 0.4 and ridge penalty 10, and walk through the complete
workflow on public cohorts.

Cohort-level eligibility and finite positive OS follow the frozen analysis. The original
studies and accession records document the complete study design, treatment context and
endpoint definitions.
