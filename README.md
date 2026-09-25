# CrossDomainAdjust

**Reusable feature encoding and partial cross-domain adjustment in R — version 2.1.0.**

[中文说明](README_zh.md) · [Detailed user guide](inst/doc/USER_GUIDE.md) ·
[中文完整手册](inst/doc/USER_GUIDE_zh.md) · [15 runnable examples](inst/examples) ·
[Data provenance](inst/extdata/README.md) · [Release validation](RELEASE_VALIDATION.md)

CrossDomainAdjust learns preprocessing on training samples and applies the same frozen
transformation to a new cohort or an individual sample. It provides four convenient
methods for cross-domain modeling, plus a lower-level composable API. It does not assume
that domain differences are entirely technical: cancer type and genetic background can
contribute both to expression differences and to prognosis.

| Method | `method` | Transformation | Feature scaling | Target domain label |
|---|---|---|---|---|
| Within-sample rank | `"rank"` | Fractional ranks over a fixed gene universe | Frozen training z-score | Not needed |
| Binary gene pairs | `"pair"` | `I(gene_i > gene_j)`; retain training prevalence 20–80% | Frozen training z-score | Not needed |
| Partial center translation | `"center"` | Subtract a fraction of a fitted domain offset | None | Known training domain required at nonzero λ |
| Partial domain projection | `"project"` | Remove a fraction of a learned low-rank direction | None | Not needed |

Raw input is assumed to have already received appropriate assay-level processing. These
functions do not turn raw sequencing counts into normalized expression. Centering and
projection use the supplied scale, without adding a gene-wise z-score in the presets.

## Installation

R >= 4.1.0. Core code uses only R's `stats` and `utils`; no compiled code or Rtools is
needed for this package. Survival examples use the recommended package `survival`.

```r
install.packages("remotes")
remotes::install_github("weikaixi/CrossDomainAdjust", upgrade = "never")
library(CrossDomainAdjust)
packageVersion("CrossDomainAdjust")
```

For an archived source tarball:

```r
install.packages("CrossDomainAdjust_2.1.0.tar.gz", repos = NULL, type = "source")
```

## A complete first example

```r
library(CrossDomainAdjust)
dat <- simulate_domain_data(p = 8, seed = 42)
train <- dat$train
test <- dat$test

prep <- fit_domain_method(train$x, train$domain, method = "project", lambda = 0.4)
x_train <- transform_cross_domain(prep, train$x)
x_test <- transform_cross_domain(prep, test$x)

# Example downstream model; preprocessing and predictive modeling are separate.
model <- lm.fit(cbind(1, x_train), train$regression)
prediction <- as.numeric(cbind(1, x_test) %*% model$coefficients)
sqrt(mean((test$regression - prediction)^2))

# A named single sample receives exactly the same frozen transformation.
one_sample <- transform_cross_domain(prep, test$x[1, ])
saveRDS(prep, "preprocessing.rds")
```

Replace `method = "project"` with `"rank"`, `"pair"` or `"center"`. For center translation,
also pass the corresponding labels to `transform_cross_domain(prep, test$x, test$domain)`.
For pair encoding, `frequency = c(0.2, 0.8)` and `max_pairs = 128` control feature growth.

## Choosing λ

`tune_cross_domain()` fits preprocessing independently inside every training fold.
The default grid is **0, 0.2, 0.4, 0.6, 0.8, 1**. Supply model-fitting, prediction and
scoring callbacks; higher scores are better. Vector, matrix, data-frame and `Surv`
outcomes are supported. See the complete [regression](inst/examples/07_regression_lambda_cv.R)
and [survival](inst/examples/09_survival_lambda_cv.R) examples.

```r
geometry <- lambda_diagnostics(prep, train$x, train$domain)
geometry$summary
```

This diagnostic uses a PCA basis fixed at λ = 0, and full-feature between-domain /
total variation. It describes geometry; it does **not** choose a prognostically optimal
λ. At λ = 1, mean-based full-rank adjustment can align the fitting-domain centroids.
It does not generally make distributions identical or preserve all useful biology.
An interior optimum is possible, not imposed. Keep the complete grid and all seeds.

## Examples included with the installed package

| No. | Scenario | Main lesson |
|---|---|---|
| 01 | Four methods | One interface, distinct representations |
| 02 | Rank normalization | Monotone invariance and frozen z-scores |
| 03 | Gene-pair features | Training-only prevalence, explicit candidates, feature caps |
| 04 | Center translation | Known domain labels and rejection of unseen centers |
| 05 | Single new sample | Frozen projection without a target cohort |
| 06 | λ visualization | Fixed PCA coordinates and matching axes |
| 07 | Regression | Fold-local λ selection and held-out RMSE |
| 08 | Binary classification | Logistic model and held-out Brier score |
| 09 | Survival | Ridge Cox callbacks and pooled C-index |
| 10 | Multiple domains | Robust centers and limited-rank projection |
| 11 | Seed sensitivity | All seeds and all λ values retained |
| 12 | Model deployment | RDS persistence and feature-schema validation |
| 13 | Advanced hybrid | Center-then-project composition, with overlapping effects |
| 14 | Real lung data | TCGA LUAD/LUSC → GSE3141 survival workflow |
| 15 | Real glioma data | TCGA LGG/GBM → GSE43378 survival workflow |

```r
source(system.file("examples", "01_four_methods.R", package = "CrossDomainAdjust"))
source(system.file("examples", "14_real_lung_survival.R", package = "CrossDomainAdjust"))
source(system.file("examples", "run_all.R", package = "CrossDomainAdjust"))
```

The first thirteen examples use explicitly synthetic educational data. The two real
examples include compact processed public-source data, a frozen 60-gene feature set
and fixed illustrative tuning; their output is not presented as a reproduction of
the manuscript's selected model estimates. 

## Reproducible use and interpretation

- Use samples in rows and unique, finite numeric gene columns; align outcomes by sample ID.
- Freeze the gene universe before deployment, and fit pair filtering/scaling only on training data.
- Group repeated samples from one patient within a fold. Use an outer holdout or nested validation.
- Rank/pair invariance applies to a strictly increasing transform within each sample, not to
  arbitrary gene-specific platform distortions.
- A frozen center needs a known domain; a frozen projection can accept an unseen label, but
  cannot guarantee removal of the unseen shift.
- A pooled C-index mixes within- and between-subtype comparisons. Report per-domain metrics
  when the intended claim concerns performance within a cancer subtype.
- Core methods do not implement ComBat, limma, quantile normalization or CORAL. Comparisons
  using those baselines belong to separate study code.

The [full guide](inst/doc/USER_GUIDE.md) documents formulas, every important parameter,
fold-local tuning, memory limits, realistic data preparation, error handling and deployment.

## Citation and license

```r
citation("CrossDomainAdjust")
```

MIT software license. Public-source worked data retain their original attribution and
source conditions; see [provenance](inst/extdata/README.md). This repository does not
claim that a peer-reviewed software article has already been published.
