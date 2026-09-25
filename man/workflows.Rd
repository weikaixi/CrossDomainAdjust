\name{fit_domain_method}
\alias{fit_domain_method}
\alias{lambda_diagnostics}
\alias{simulate_domain_data}
\title{Four-Method Presets, Lambda Diagnostics and Educational Data}
\description{Convenient entry points for frozen rank and pair encodings, partial center translation and domain projection, with fixed-axis diagnostics and reproducible synthetic examples.}
\usage{
fit_domain_method(x, domain, method = c("rank", "pair", "center", "project"),
  lambda = 1, ...)
lambda_diagnostics(object, x, domain,
  lambda_grid = c(0, 0.2, 0.4, 0.6, 0.8, 1))
simulate_domain_data(n_train = 180, n_test = 90, p = 20,
  domains = 3, seed = 42)
}
\arguments{
\item{x}{Finite numeric samples-by-features matrix with unique column names.}
\item{domain}{One domain label per row.}
\item{method}{One of the four presets. Rank and pair use frozen training z-scores after encoding; center and project use raw input without extra standardization.}
\item{lambda}{Strength in [0,1] for center or project. Ignored by rank and pair presets, which have no subsequent domain correction.}
\item{...}{Additional arguments to \code{fit_cross_domain}, for example frequency, max_pairs or anchor. Representation, correction and standardize cannot be overridden through a preset.}
\item{object}{A fitted center or project pipeline.}
\item{lambda_grid}{Diagnostic strengths in [0,1].}
\item{n_train, n_test}{Numbers of simulated training and test samples, each at least domains.}
\item{p}{Number of simulated features, at least four.}
\item{domains}{Number of simulated domains, at least two.}
\item{seed}{A nonnegative integer seed. The caller's random-number state is restored on return.}
}
\details{Diagnostics fit PCA once to lambda-zero values and project every transformed matrix into those fixed axes, without extra feature scaling. Domain R-squared is full-feature between-domain sum of squares divided by total sum of squares, using sample-count weights. It is not an average of gene-wise standardized R-squared values. Centroid distances are Euclidean (not divided by the square root of the number of features). Mean-based, full-rank corrections can align fitting-domain means at lambda one, but need not equalize covariance, remove unseen shifts or preserve prognostic information. Simulated data contain both shared domain offsets and domain-associated biology; they are educational, not a demonstration that any method must win.}
\value{The preset returns a cross_domain_fit. Diagnostics return a list with summary, coordinates, pca_center, pca_rotation, baseline_variance_explained and note. Simulation returns train and test lists, each containing x, domain, regression, binary and survival (time and event), plus a description.}
\examples{
d <- simulate_domain_data(n_train = 60, n_test = 30, p = 6)
fit <- fit_domain_method(d$train$x, d$train$domain, "project", lambda = .4)
z <- transform_cross_domain(fit, d$test$x)
lambda_diagnostics(fit, d$train$x, d$train$domain)$summary
}
