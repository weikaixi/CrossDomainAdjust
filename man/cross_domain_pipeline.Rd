\name{fit_cross_domain}
\alias{fit_cross_domain}
\alias{transform_cross_domain}
\alias{diagnose_cross_domain}
\title{Fit a Feature Encoding, Standardization and Domain Adjustment Pipeline}
\description{Freeze a feature representation, training-set pair screening, optional feature standardization and domain adjustment in a single reusable object.}
\usage{
fit_cross_domain(x, domain, representation = c("raw", "rank", "pair"),
  correction = c("none", "center", "project", "hybrid"), lambda = 1,
  frequency = c(0.2, 0.8), pairs = NULL, max_pairs = NULL,
  max_candidates = 1000000, chunk_size = 1000,
  anchor = "domain_mean", center = "mean", trim = 0.1,
  n_components = NULL, standardize = c("none", "zscore"))
transform_cross_domain(object, x, domain = NULL, lambda = NULL)
diagnose_cross_domain(object, x, domain)
}
\arguments{
\item{x}{Numeric samples-by-features matrix, numeric data frame or named sample vector.}
\item{domain}{Domain labels, one per row.}
\item{representation}{Raw values, within-sample fractional ranks or binary gene-pair features.}
\item{correction}{Identity, center translation, low-rank projection or sequential hybrid.}
\item{lambda}{Strength in [0,1]; hybrid accepts one shared or two separate strengths. In transformation, NULL uses the stored value; an explicit override does not modify the fitted object.}
\item{frequency}{Inclusive training prevalence bounds for retained pairs, default 20 to 80 percent ones.}
\item{pairs}{Optional explicit two-column matrix of gene names.}
\item{max_pairs}{Optional cap on retained pairs; see \code{fit_feature_encoder}.}
\item{max_candidates}{Maximum number of candidate pairs allowed before evaluation.}
\item{chunk_size}{Candidate comparisons per block.}
\item{anchor, center, trim, n_components}{Passed to \code{fit_domain_adjuster}.}
\item{standardize}{Optional feature-wise scaling, fitted after encoding and before domain adjustment. \code{"none"} preserves encoded values; \code{"zscore"} uses training means and sample SDs. Constant features use scale 1. See \code{\link{fit_feature_scaler}}.}
\item{object}{A fitted cross_domain_fit.}
}
\details{
The order is encoding, optional feature standardization, then domain correction.
Rank encoding already performs within-sample normalization as
\eqn{(r-1)/(p-1)}, where \eqn{r} uses average ranks for ties and \eqn{p} is
the frozen gene-universe size. This differs from optional z-score scaling,
which estimates each feature's mean and sample SD across training samples.
Standardized rank and pair features need not remain in \eqn{[0,1]}.

New data never change the feature universe, pair prevalence screen, feature
means or scales, or fitted domain centers. A single new sample uses the same
frozen parameters. Version 2.0.0 fitted objects without a scaler retain
identity scaling. When \code{standardize = "zscore"} is passed through
\code{\link{tune_cross_domain}}, each fold fits its own scaler using only
fold-training samples; heldout samples do not contribute to scaling.

Pair selection is pooled over training samples; it does not guarantee balance
within every domain or remove redundant pairs. Corrected pair features are
continuous, not binary. Rank and pair encodings are alternatives, not
sequential steps. For rank encoding the full fitted gene universe must be
supplied. The same requirement is deliberately applied to pair encoders.

Standardization changes the geometry used to estimate domain directions and
is optional rather than inherently beneficial. When enabled, a numeric custom
anchor must be expressed in standardized feature coordinates. Diagnostics
compare encoded-and-standardized values against corrected values in the same
coordinate system; they are descriptive and do not establish biological
validity or prediction performance.
}
\value{\code{fit_cross_domain} returns a cross_domain_fit containing an encoder, scaler and adjuster. \code{transform_cross_domain} returns a numeric matrix. \code{diagnose_cross_domain} returns centroid-distance and within-domain variation summaries and per-sample displacement caused by domain correction, measured after optional scaling.}
\examples{
set.seed(20)
x <- matrix(rnorm(240), 40, 6, dimnames = list(NULL, paste0("g", 1:6)))
d <- rep(c("A", "B"), each = 20)
obj <- fit_cross_domain(x, d, representation = "pair",
  correction = "project", lambda = 0.5, frequency = c(0.2, 0.8))
transform_cross_domain(obj, x[1:2, , drop = FALSE])
diagnose_cross_domain(obj, x, d)$summary
rank_obj <- fit_cross_domain(x, d, representation = "rank",
  standardize = "zscore", correction = "project", lambda = 0.5)
transform_cross_domain(rank_obj, x[1, ])
}
