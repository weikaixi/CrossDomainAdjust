# CrossDomainAdjust 2.1.0 — detailed user guide

## Contents

1. Purpose and scope
2. Installation and file layout
3. Preparing an expression matrix
4. Four methods and equations
5. Pair screening and computational limits
6. Frozen scaling and feature schemas
7. Choosing domains, anchors and projection rank
8. Lambda selection without target leakage
9. Regression, classification and survival callbacks
10. Pooled and domain-specific evaluation
11. Fixed-PCA lambda diagnostics
12. Real lung and glioma workflows
13. Repeated splits and sensitivity analysis
14. Saving and deploying a complete prediction pipeline
15. Advanced composition and low-level functions
16. Troubleshooting and compatibility
17. Reproducibility checklist and examples index

## 1. Purpose and scope

CrossDomainAdjust fits reusable transformations for cross-domain prediction. Four presets
represent the workflows used in the associated exploratory lung and glioma analyses.
The package can also process other finite numeric feature matrices; gene-expression
terminology is used because it is the motivating application.

A domain label encodes a grouping supplied by the analyst. It may reflect histology,
platform, laboratory, patient ancestry or a combination of factors. Neither a difference
between domain means nor a fitted projection direction identifies a purely technical
effect. In particular, histology-associated expression can carry prognostic information.

The package separates preprocessing from the downstream model. It neither performs
differential expression analysis nor automatically imputes missing genes, normalizes
RNA-seq counts, maps probes, selects prognostic genes or fits a mandatory prediction
algorithm. Those operations need their own reproducible, training-only rules.

## 2. Installation and file layout

```r
install.packages("remotes")
remotes::install_github("weikaixi/CrossDomainAdjust", upgrade="never")
library(CrossDomainAdjust)
packageVersion("CrossDomainAdjust")
help(package="CrossDomainAdjust")
```

The package requires R >= 4.1.0. Core dependencies are bundled with R. Install `survival`
if your installation does not already include it: `install.packages("survival")`.
There is no compiled code. A local source archive can be installed with
`install.packages("CrossDomainAdjust_2.1.0.tar.gz", repos=NULL, type="source")`.

After installation, `system.file("examples", package="CrossDomainAdjust")` contains
15 numbered scripts and `run_all.R`; `system.file("doc", package="CrossDomainAdjust")`
contains the two guides. `system.file("extdata", package="CrossDomainAdjust")`
contains the two real worked-data files and their provenance.

For long-term reproducibility, record the package version and Git commit or archive
checksum. The default GitHub branch may change after this release.

## 3. Preparing an expression matrix

Input has **samples in rows and genes/features in columns**. Column names must be unique
and nonempty. Values must be numeric and finite. Named vectors are accepted for a single
new sample. A numeric data frame is accepted, but a sample-ID text column must first be
moved to row names. Sparse and delayed matrices are not currently supported.

```r
# Starting from a conventional gene-by-sample table:
# expression <- read.delim("expression.tsv", row.names=1, check.names=FALSE)
# x <- t(as.matrix(expression))
# clinical <- read.csv("clinical.csv")
# stopifnot(!anyDuplicated(clinical$sample_id))
# clinical <- clinical[match(rownames(x), clinical$sample_id), ]
# stopifnot(!anyNA(clinical$sample_id), identical(rownames(x), clinical$sample_id))
# domain <- clinical$histology
```

Do not silently interpret absent genes as zero. Resolve identifiers, duplicate probes
and required genes upstream. If imputation is necessary, estimate its rule using only
training data and apply that same rule to new data. Appropriate assay-level transforms
(for example a predefined log-expression scale) must be handled before the package.

The same frozen gene universe is required at prediction time, even when a pair encoder
ultimately uses fewer genes. Columns may arrive in a different order; the package aligns
them by name. Additional finite numeric columns are ignored after input validation.
Missing fitted columns are errors. Row order is not inferred from an outcome vector:
the caller must align outcomes and domains to the matrix before fitting.

## 4. Four methods and equations

For a row vector x, p fitted features, a domain d, a reference anchor a and a strength
0 <= λ <= 1:

### Rank

`r_j = (rank(x_j) - 1) / (p - 1)`, with average ranks for ties. Ranks are computed within
each sample across the same p genes, not across samples for each gene. At least two genes
are required. The rank preset then uses the training mean and sample SD of each rank
feature to obtain a frozen z-score. New samples do not alter these values.

Strictly increasing transformations applied to all genes within a sample preserve their
ordering. Thus an additive sample offset, a positive sample multiplier or another common
monotone transformation does not change its ranks. Gene-specific platform distortions,
ties, missing genes and a changed gene universe can change ranks. This explains a useful
invariance, not a universal guarantee of improved external prediction.

```r
d <- simulate_domain_data(p=8)
rank_fit <- fit_domain_method(d$train$x, d$train$domain, "rank")
rank_test <- transform_cross_domain(rank_fit, d$test$x)
```

### Pair

`b_ij = I(x_i > x_j)`. A tie returns zero. Each unordered pair is evaluated once in
the orientation determined by the original training-column order. Explicit reversed
pairs are canonicalized to that order; duplicated/reversed redundant pairs are errors.
The proportion of ones is computed in training data, and inclusive bounds retain
`0.2 <= mean(b_ij) <= 0.8` by default. The pair preset then applies frozen training z-scores.
After scaling, the output is continuous and should not be described as binary input.

```r
pair_fit <- fit_domain_method(d$train$x, d$train$domain, "pair",
  frequency=c(.2,.8), max_pairs=128)
pair_test <- transform_cross_domain(pair_fit, d$test$x)
pair_fit$encoder$frequencies
```

### Center translation

Let μ_d be a fitted domain center, δ_d = μ_d - a. The output is
`x* = x - λ δ_d`. The center preset uses the input expression scale directly, with no
additional z-score. With arithmetic-mean centers, fitting-domain means become
`a + (1 - λ) δ_d`. At λ=1 they coincide. Subtracting a constant within a domain leaves
its covariance unchanged. For new samples, a known fitted domain label is required
whenever λ is nonzero; the package does not estimate a new center from the target cohort.

```r
center_fit <- fit_domain_method(d$train$x, d$train$domain, "center", lambda=.4)
center_test <- transform_cross_domain(center_fit, d$test$x, d$test$domain)
```

### Domain projection

Stack the δ_d rows and obtain an orthonormal basis B from their singular value
decomposition. The output is `x* = x - λ (x - a) B B'`. Only the low-rank products are
computed; a dense p-by-p projection matrix is not allocated. Components orthogonal to
B are preserved; components in B are multiplied by 1-λ around the anchor.

Unlike a within-domain translation, projection also changes within-domain variation
along B. Frozen projection can be applied to a single new sample without knowing its
domain label. It removes only components in the learned subspace; a new platform shift
outside that subspace remains. With default mean centers and all fitted directions,
training-domain means coincide at λ=1, but covariances need not.

```r
projection_fit <- fit_domain_method(d$train$x, d$train$domain, "project", lambda=.4)
projection_test <- transform_cross_domain(projection_fit, d$test$x)
```

λ has no role in rank or pair presets, which have no subsequent center/projection step.
It is ignored there. Advanced compositions are explicit separate calls, not the default.

## 5. Pair screening and computational limits

With p genes there are p(p-1)/2 possible unordered pairs. Frequency filtering reduces
the retained output, but does not remove the computational cost of considering candidates.
The implementation evaluates blocks (`chunk_size=1000` by default), keeping memory below
that of a full sample-by-all-candidates matrix.

`max_candidates=1000000` is a guard evaluated before enumeration. It raises an error
if the candidate space is too large. `max_pairs` caps the retained features **after**
frequency screening; it does not make a huge candidate space cheap to enumerate.
When capped, pairs closest to a training prevalence of 0.5 are prioritized; ties follow
candidate order. This is not a supervised association test, independence screen or
redundancy-removal algorithm. Correlated/redundant retained pairs can remain.

For thousands of genes, define a training-derived gene panel or a biologically specified
two-column candidate matrix first. A 60-gene panel contains 1770 candidates; 20,000 genes
would contain almost 200 million. Increasing the guard alone is rarely a good remedy.

```r
candidate_pairs <- rbind(c("gene1","gene2"), c("gene3","gene4"))
small <- fit_domain_method(d$train$x, d$train$domain, "pair",
  pairs=candidate_pairs, frequency=c(0,1))
```

The prevalence is pooled over training samples, not required separately in every domain.
A large domain can dominate it. If no pair passes, fitting fails; reconsider candidate
definitions or bounds using training data instead of silently relaxing rules on test data.

## 6. Frozen scaling and feature schemas

There are two distinct operations in the rank preset: sample-level fractional ranking
and feature-level z-scoring. The latter is optional in the general API and is included
in the study-aligned rank/pair presets. It makes feature scales more comparable for
scale-sensitive predictors; it does not inherently improve every model.

The training feature mean m_j and sample SD s_j give `(z_j - m_j) / s_j` for every future
sample. A constant training feature uses denominator one. Training-constant features
are retained unless removed upstream. Do not call `scale(test_x)` independently, because
that would estimate target-specific parameters and change single-sample behavior.

```r
# Custom rank-only encoding is available, although not one of the four study presets.
custom <- fit_cross_domain(d$train$x, d$train$domain,
  representation="rank", correction="none", standardize="none")
```

Changing the gene universe changes fractional ranks and candidate pairs. Even for shared
genes, encoding a newly chosen subset is not equivalent to subsetting an encoded result.
Freeze the schema during development and distribute it with the model.

## 7. Choosing domains, anchors and projection rank

Choose a grouping that matches the scientific question. The lung worked data use LUAD
and LUSC; glioma uses LGG and GBM. These labels describe biological strata, not GEO study
IDs. A GEO LUAD sample can use the frozen LUAD center learned from TCGA because its label
is known; this does not mean a GEO-specific platform shift has been estimated or removed.

If the task is technical batch correction, supply actual technical batch labels. A
new site then becomes an unseen domain for center translation and is rejected at
nonzero λ. Projection remains computable, but its generalization must be evaluated.
If biology and batch are perfectly confounded, the observed means alone cannot identify
which component should be removed.

`anchor="domain_mean"` is the equal-weight mean of fitted domain centers.
`anchor="sample_mean"` weights those centers by domain sample counts. With default
mean centers, this is the overall sample mean. With median/trimmed centers it is a
weighted average of robust centers, not necessarily the raw arithmetic sample mean.
A custom named vector is also allowed, in the coordinates passed to the adjuster.

`center="mean"` is the default; `"median"` and `"trimmed"` are available, with
`trim=.1` removing that fraction from each tail for trimmed means. Robust centers can
reduce sensitivity to outliers; they do not guarantee alignment of arithmetic means.

For K domains and a built-in anchor, fitted rank is at most K-1. `n_components=NULL`
retains all nonzero directions; an integer between zero and fitted rank retains leading
directions only. A custom anchor outside the centers' affine hull can give rank K.
Truncating the basis may prevent complete centroid alignment at λ=1.

## 8. Lambda selection without target leakage

Use training outcomes and inner folds to select λ, then evaluate a final frozen pipeline
on independent samples. `tune_cross_domain()` accepts fold labels and callbacks:

```r
fit_model <- function(x, y) {
  z <- cbind(1, x)
  solve(crossprod(z) + diag(c(0, rep(1, ncol(x)))), crossprod(z, y))
}
predict_model <- function(model, x) as.numeric(cbind(1, x) %*% model)
score <- function(y, prediction, domain) -mean((y-prediction)^2)
set.seed(7)
folds <- sample(rep(1:3, length.out=nrow(d$train$x)))
tuned <- tune_cross_domain(d$train$x, d$train$domain, d$train$regression,
  folds, fit_model, predict_model, score, correction="project",
  lambda_grid=c(0,.2,.4,.6,.8,1))
tuned$best_lambda
tuned$scores
tuned$fold_scores
```

For each fold, the encoder, pair screen, scaler, centers and basis are fitted only on
the remaining samples. The same fold-local preprocessing is applied across strengths.
The downstream callback is refitted for every λ. Scores are averaged equally across
folds, and the **largest** mean wins. Exact ties choose the smaller λ. Erroring fits or
nonfinite scores stop the procedure rather than silently removing difficult folds.

Return values are `best_lambda`, `scores`, `fold_scores`, `fit` and a reminder `note`.
The returned `fit` is a preprocessing object refitted to all supplied development data;
it is not a final predictive model. Fit that separately on its transformed training set.
The inner selection score is optimistic as a performance estimate: use an outer holdout
or repeat this whole selection inside outer folds.

All observations from a patient must share a fold. Maintain the necessary training
domains in each fold for center translation. Leave-one-domain-out validation can use
projection; nonzero frozen centering cannot infer a center for a held-out domain.
Gene selection performed before this function is not automatically repeated inside it;
supervised feature selection belongs inside the corresponding training-only validation
workflow. The helper does not tune λ and a model penalty jointly; implement an explicit
nested/grid workflow if both are being selected, and record the full search.

Do not exclude λ=0 or λ=1 to obtain an interior optimum. Partial correction may be useful,
but a valid analysis must permit data to favor either endpoint.

## 9. Regression, classification and survival callbacks

Regression examples use negative MSE, so larger is better. For binary classification,
the same principle permits negative Brier score or a properly defined AUC. Thresholds
must be fixed or selected using training data; test-set threshold optimization inflates
apparent performance. Example 08 uses a fixed 0.5 threshold and reports Brier score too.

`y` may be a vector, matrix, data frame or `survival::Surv` object. Matrix-like outcomes
are subset by rows with dimensions retained. Callbacks must handle the supplied type.
Example 09 uses `data.frame(time,event)` and a ridge Cox model:

```r
library(survival)
fit_cox <- function(x, y) {
  z <- x; time <- y$time; event <- y$event
  coef(coxph(Surv(time,event) ~ ridge(z, theta=10, scale=FALSE)))
}
predict_cox <- function(model, x) as.numeric(x %*% model)
cindex <- function(y, prediction, domain) {
  time <- y$time; event <- y$event; risk <- prediction
  as.numeric(concordance(Surv(time,event) ~ risk, reverse=TRUE)$concordance)
}
```

The survival examples use positive time, event 1=death and event 0=censored. Time units
must be consistent. `reverse=TRUE` means a higher Cox linear predictor indicates greater
risk and earlier events. If using `Surv` directly, adapt the callbacks to receive that
object rather than reading `$time` and `$event`. A fold with no comparable survival pairs
cannot yield a meaningful C-index; redesign folds prospectively instead of dropping it.

The examples use a fixed penalty for clarity. A research analysis should justify or
select its penalty within development data. A constant model penalty can interact with
feature scale, so record whether scaling or Cox ridge internal scaling was enabled.

## 10. Pooled and domain-specific evaluation

A pooled metric for a single external cohort uses all eligible samples in that cohort,
including mixed histologies when present. It is distinct from concatenating all external
cohorts into one pooled dataset, and from averaging subtype-specific metrics.

Pooled survival C-index includes both within-subtype and between-subtype comparable pairs.
A model can obtain discrimination from subtype-related differences in prognosis. To
assess within-subtype discrimination, report LUAD/LUSC or LGG/GBM separately as well.
An unweighted mean over subtype metrics is a domain-macro summary; small strata and
non-estimable C-indices must be explicitly reported. Neither quantity equals pooled C-index.



## 11. Fixed-PCA lambda diagnostics

`lambda_diagnostics(fit, x, domain)` first transforms at λ=0, fits mean-centered PCA
without extra gene scaling, then uses those fixed axes for every λ. The returned
`coordinates` table has sample, domain, lambda, PC1 and PC2; `summary` has strength,
domain R-squared, mean centroid distance and mean displacement. The other entries
preserve the baseline PCA center, rotation, variance fractions and interpretation note.

Domain R-squared is `sum_d n_d ||mean_d - overall_mean||^2 / sum_i ||x_i-overall_mean||^2`
in the full transformed feature space. It is not limited to the first two PCs and is not
the average of gene-wise z-scored R-squared values. The latter would define a different
metric. Centroid distance is ordinary Euclidean distance; the older
`diagnose_cross_domain()` function reports distance divided by sqrt(feature count).

```r
diag <- lambda_diagnostics(projection_fit, d$train$x, d$train$domain)
diag$summary
head(diag$coordinates)
source(system.file("examples","06_lambda_fixed_pca.R",package="CrossDomainAdjust"))
```

At full mean-based adjustment of the fitting data, centroid R-squared can become zero
while covariance, multimodality and other distributional differences remain. Fixed
external samples are not guaranteed zero centroid difference even at λ=1. A diagram
showing means approaching each other therefore does not prove removal of all batch effects
or preservation of prognosis. 

## 12. Real lung and glioma workflows

Two compact processed public-source worked datasets are included. Lung has 996 TCGA
training samples (LUAD/LUSC) and 110 GSE3141 samples. Glioma has 664 TCGA training samples
(LGG/GBM) and 50 GSE43378 samples. Each has 60 named genes frozen from the existing study
workflow. Sources, processing provenance and limitations are in `extdata/README.md`.

```r
lung <- readRDS(system.file("extdata","lung_worked_data.rds",package="CrossDomainAdjust"))
names(lung)
dim(lung$train$x)
table(lung$train$domain)
head(lung$train$survival)
source(system.file("examples","14_real_lung_survival.R",package="CrossDomainAdjust"))
source(system.file("examples","15_real_glioma_survival.R",package="CrossDomainAdjust"))
```

Both scripts compare unadjusted input with the four methods, apply the same training-only
preprocessing to the external cohort, fit ridge Cox models and print pooled C-indices.
They use illustrative λ=.4 and θ=10, not the manuscript's selected model configuration.
They do not rerun feature discovery or claim de novo external confirmation. They are
practical worked applications of the package, not an unbiased contest in which a specific
method is guaranteed to outperform every baseline.

To use your own cohorts, replace the data lists while preserving feature names, domain
definitions and endpoint coding. Never infer a target center from validation outcomes.
If a target gene is absent, adapt the feature panel during a new development analysis
and refit the full pipeline; do not silently change a deployed model's feature universe.

## 13. Repeated splits and sensitivity analysis


Repeated observations on the same held-out cohort are not independent external
validations. Repeatedly choosing cohorts based on which method wins is performance-based
selection. Use the complete candidate-cohort log and separate exploratory comparisons
from future independent confirmation.

## 14. Saving and deploying a complete prediction pipeline

Save the preprocessing object, downstream coefficients/model, exact feature schema,
domain definitions, selected λ, model penalty, training preprocessing steps and package
version together. Record endpoint/time units for survival models.

```r
prep <- projection_fit
z <- transform_cross_domain(prep, d$train$x)
model <- fit_model(z, d$train$regression)
bundle <- list(preprocessing=prep, model=model, genes=colnames(d$train$x),
               package_version=as.character(packageVersion("CrossDomainAdjust")))
file <- tempfile(fileext=".rds")
saveRDS(bundle, file)
restored <- readRDS(file)
z_new <- transform_cross_domain(restored$preprocessing, d$test$x[1, ])
prediction <- predict_model(restored$model, z_new)
unlink(file)
```

An explicit `lambda=` in `transform_cross_domain()` temporarily overrides the stored
strength without altering the object. This is useful for diagnostics. Do not change λ
at deployment while reusing coefficients trained under another λ; refit/evaluate the
downstream model for the same transformation.

## 15. Advanced composition and low-level functions

`fit_cross_domain()` exposes `representation` (raw/rank/pair), `standardize` (none/zscore),
and `correction` (none/center/project/hybrid). Processing order is encoding, scaling,
then correction. The default general API uses no feature standardization; the rank/pair
presets explicitly enable it. Preset-fixed arguments cannot be overridden through `...`.

Hybrid applies center translation followed by projection with the same frozen direction.
`lambda=c(center_strength, projection_strength)` specifies separate strengths. A scalar
applies to both. The two operations overlap and are not independent evidence of benefit.
Hybrid is an advanced interface, not a fifth method added to the main four-method study.

Low-level exported functions include `fit_feature_encoder()` / `encode_features()`,
`fit_feature_scaler()` / `scale_features()`, `fit_domain_adjuster()` / `adjust_features()`,
and the retained projection wrappers `fit_domain_projector()`, `project_features()` and
`project_new_sample()`. Use their R help pages for exact argument definitions. The single
sample projector wrapper returns a named vector; pipeline transformation returns a matrix.

## 16. Troubleshooting and compatibility

| Symptom | Meaning and remedy |
|---|---|
| Missing feature names | Set unique gene names on matrix columns; transpose if genes are rows. |
| Missing fitted genes | Supply the full training universe or refit during development; do not zero-fill silently. |
| Nonfinite values | Resolve NA/Inf upstream using documented training-fitted rules. |
| Candidate limit exceeded | Restrict candidate genes/pairs before enumerating; `max_pairs` alone is insufficient. |
| No pairs retained | Inspect training prevalence and candidate definitions. |
| Unknown domain in center | Use a known training domain label; projection can accept unseen labels but offers no guarantee. |
| Rank Z-values outside [0,1] | Expected after feature z-scoring; only fractional ranks lie in [0,1]. |
| A larger λ performs worse | Domain directions may contain useful biology; validation, not visual overlap, selects strength. |
| Callback score is NA | Check sample/fold alignment, events and comparable pairs; failures are not ignored. |
| Different results after shuffling columns | Named columns should align; verify feature names, input values and object version. |
| External cohort means remain separated at λ=1 | Frozen transforms do not estimate new external domain shifts. |

Version 2.1.0 preserves frozen 2.0.1 transformation behavior. New diagnostics and presets
do not change the manuscript's previously frozen models. Legacy 2.0.0 pipelines without
a scaler use identity scaling. Serialized 0.1.x projector objects should be refitted.
Keep historical archives to reproduce older analyses rather than overwriting them.


Run all numbered scripts with:

```r
source(system.file("examples","run_all.R",package="CrossDomainAdjust"))
```

01 four methods; 02 rank/scaling; 03 pair screening; 04 known-domain centering;
05 unseen-sample projection; 06 fixed-PCA λ plot; 07 regression CV; 08 classification;
09 survival CV; 10 robust/multidomain; 11 full seed grid; 12 save/load/schema;
13 hybrid; 14 real lung; 15 real glioma. Synthetic results illustrate execution, not
clinical performance. Public-source worked data illustrate use, not a new confirmatory study.
