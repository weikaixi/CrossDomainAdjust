\name{tune_cross_domain}
\alias{tune_cross_domain}
\title{Select Adjustment Strength Using Fold-Local Preprocessing}
\description{Fit the encoder and adjuster independently in each supplied training fold, select lambda using validation scores and refit a final preprocessing object on all supplied data.}
\usage{
tune_cross_domain(x, domain, y, folds, fit_model, predict_model, score,
  lambda_grid = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  representation = "raw", correction = "project", ...)
}
\arguments{
\item{x}{Samples-by-features numeric matrix.}
\item{domain}{One domain label per sample.}
\item{y}{A nonmissing outcome vector, matrix, data frame or survival::Surv object with one row per sample. Structured outcomes are subset by row with dimensions retained.}
\item{folds}{User-supplied validation fold labels. All observations from the same patient should share a fold.}
\item{fit_model}{Function of (x, y) returning a model. Any scaling or model selection in this callback must use only these training inputs.}
\item{predict_model}{Function of (model, x) returning validation predictions.}
\item{score}{Function of (y, predictions, domain) returning one finite scalar; larger is better.}
\item{lambda_grid}{Candidate scalar strengths in [0,1]. Exact ties choose the smaller strength. Hybrid candidates apply the same scalar to both operations.}
\item{representation, correction}{Passed to \code{fit_cross_domain}.}
\item{...}{Additional arguments for \code{fit_cross_domain}.}
}
\details{Fold scores have equal weight in selection. A score callback can compute domain-macro performance. The function does not create an outer test split and selection scores are not unbiased performance estimates. Reserve an untouched outer holdout or nest this entire function in an outer validation loop. Leave-one-domain validation is compatible with frozen projection, but not nonzero center adjustment unless the held-out domain has an independently estimated center; this function never estimates one from validation data. Failures are errors rather than silently discarded folds. No final predictive model is returned; fit it on the returned pipeline's transformed training data.}
\value{A list containing best_lambda, scores, fold_scores, fit (the final preprocessing object) and an evaluation reminder.}
