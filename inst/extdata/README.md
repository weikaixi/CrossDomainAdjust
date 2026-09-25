# Public-source worked data

These small derived data files support examples 14 and 15. They contain previously
curated expression for the study's frozen 60-gene panels, public sample identifiers,
histology labels and overall-survival endpoints. No new patient recruitment or private
medical records are involved. The MIT software license does not replace upstream data
conditions or attribution obligations.

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
  in the original curation; this package export does not repeat sample selection.
- [GSE3141](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3141): GPL570
  primary lung tumors, processed MAS5 signals. The existing workflow applied log2
  to positive deposited signals, mapped probes consistently and aggregated by gene.
  Of 111 arrays, GSM70230 lacked survival time, leaving 110 with usable OS.
- [GSE43378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43378): processed
  public glioma microarray expression and public survival annotations, using the
  existing gene mapping and eligibility curation. The export retains all 50 eligible
  samples; it does not renormalize the curated expression.

The frozen clinical files express OS in months. Original day-based endpoints were
converted with 365.25/12 days per month in the existing curation. Both data files are
subsets of the existing processed matrices, not raw assay releases. No extra gene-wise
z-scoring, pair encoding or domain adjustment has been applied to these `x` matrices.
The functions in the examples perform those steps after loading the data.

`source_manifest.csv` records source-file checksums and exported data-file checksums.
`feature_panels.csv` lists the frozen genes in order. These trace this package export;
they are not a claim to reproduce the entire upstream raw-assay processing pipeline.

## Intended use and limitations

The real examples use fixed illustrative lambda 0.4 and ridge penalty 10. They do not
reproduce the manuscript's chosen parameters, CV curves, nine-baseline comparison or
frozen performance table. The 60-gene panels and cohorts were established and inspected
in earlier exploratory work. Thus they are examples of using the package on real data,
not newly untouched external validation or proof that a particular method must win.

No demographic attribute absent from the public source is inferred. Cohort-level
eligibility and finite positive OS follow the frozen analysis. Users should consult
the original studies and public accession records for their complete study design,
treatment context, consent/data conditions and endpoint limitations, and cite the
original investigators when reusing these data.
