\name{fit_domain_adjuster}
\alias{fit_domain_adjuster}
\alias{adjust_features}
\alias{fit_domain_projector}
\alias{project_features}
\alias{project_new_sample}
\title{Fit and Apply Frozen Domain Adjustments}
\description{Learn training-domain centers and a low-rank domain subspace, then apply partial centering or projection without refitting on new samples.}
\usage{
fit_domain_adjuster(x, domain,
  method = c("project", "center", "hybrid", "none"),
  anchor = "domain_mean", center = c("mean", "median", "trimmed"),
  trim = 0.1, n_components = NULL)
adjust_features(adjuster, x, domain = NULL, lambda = 1)
fit_domain_projector(x, domain, anchor = "domain_mean")
project_features(projector, x, lambda = 1)
project_new_sample(projector, sample, domain = NULL, lambda = 1)
}
\arguments{
\item{x}{Finite numeric matrix with samples in rows and unique named features in columns; numeric data frames and named sample vectors are also accepted. Extra prediction columns are ignored; missing fitted columns are errors.}
\item{domain}{One nonempty domain label per sample. Required for fitting and for nonzero centering. Optional labels supplied for projection are checked for length but need not be seen in training.}
\item{method}{Projection, center translation, sequential center-then-project, or identity.}
\item{anchor}{\code{"domain_mean"} for equal-weight mean of fitted centers; \code{"sample_mean"} for sample-count-weighted mean of fitted centers; or a named custom vector. With robust centers the second option is not the arithmetic mean of all observations.}
\item{center}{Coordinatewise center estimator.}
\item{trim}{Fraction trimmed from each tail for trimmed centers.}
\item{n_components}{Number of leading singular directions to retain, from zero to estimated rank; NULL retains all.}
\item{adjuster}{Object returned by \code{fit_domain_adjuster}.}
\item{projector}{Object returned by \code{fit_domain_projector}.}
\item{sample}{Exactly one named sample.}
\item{lambda}{Strength between zero and one. Hybrid accepts two values, in center/project order; a scalar is applied to both.}
}
\details{
For row vector z, center adjustment is z - lambda * (mu[d,] - a).
Projection is z - lambda * (z - a) B B', where B has orthonormal columns
from the singular value decomposition of center offsets. A dense feature-by-feature
projection matrix is never materialized. A built-in anchor gives rank at most
K-1 for K domains. A custom anchor outside the affine hull can give rank K.
Full projection removes only the estimated subspace, not all batch effects.
Hybrid first translates and then projects using the same frozen basis; these
operations overlap and are not assumed to provide independent benefits.
An unseen-domain individual can receive the frozen projection, but removal of
unseen shifts is not guaranteed. Nonzero centering rejects unseen domains;
no target-domain mean is estimated silently. Missing values must be handled
upstream with training-fitted procedures. The old projector function names are
retained, but old serialized 0.1.x objects should be refitted with this version.
}
\value{Fitters return a domain_adjuster (or domain_projector) list with centers, anchor, basis, rank, feature names and domain labels. Matrix transformation returns a numeric matrix. The single-sample wrapper returns a named numeric vector.}
\examples{
set.seed(10)
x <- matrix(rnorm(120), 30, 4, dimnames = list(NULL, paste0("g", 1:4)))
d <- rep(c("A", "B", "C"), each = 10)
x[d == "B", ] <- x[d == "B", ] + 1
obj <- fit_domain_projector(x, d)
z <- project_features(obj, x, lambda = 0.5)
project_new_sample(obj, x[1, ], domain = "new_site", lambda = 0.5)
}
