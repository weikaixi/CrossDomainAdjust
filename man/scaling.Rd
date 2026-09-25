\name{fit_feature_scaler}
\alias{fit_feature_scaler}
\alias{scale_features}
\title{Fit and Apply Frozen Feature-Wise Standardization}
\description{
Learns feature-wise means and standard deviations from training samples and
reuses them without updating on new data. Identity scaling is the default.
}
\usage{
fit_feature_scaler(x, method = c("none", "zscore"))
scale_features(scaler, x)
}
\arguments{
  \item{x}{Finite numeric matrix or numeric data frame with samples in rows
  and uniquely named features in columns, or a named numeric vector for a
  single sample. Z-score fitting requires at least two training samples;
  applying a fitted scaler accepts a single new sample. All fitted feature
  names must occur in new data; additional columns are ignored after matching.}
  \item{method}{\code{"none"} for identity scaling, or \code{"zscore"} for
  feature-wise standardization with training means and sample standard deviations.}
  \item{scaler}{An object returned by \code{fit_feature_scaler}.}
}
\details{
Z-scores are computed as \eqn{(x_j-\mu_j)/s_j}, where \eqn{\mu_j} and
\eqn{s_j} are estimated across the training samples for feature \eqn{j}.
The standard deviation uses denominator \eqn{n-1}. For a constant training
feature the stored scale is 1, so training values become zero. A future value
different from that training constant remains a nonzero deviation; it is not
forced to zero or used to update the fitted center or scale.

This operation is distinct from \code{representation = "rank"} in
\code{\link{fit_feature_encoder}}, which ranks genes within each sample and
maps ranks to \eqn{[0,1]} as \eqn{(r-1)/(p-1)}, using average ranks for ties.
Feature-wise z-score scaling can be applied after raw, rank, or binary pair
encoding. Scaled ranks and pairs need not remain within \eqn{[0,1]}.

No statistics are estimated from new or validation samples. In a cross-domain
pipeline this optional scaling is fitted before domain correction; it changes
the geometry in which domain directions are estimated. Standardization does
not guarantee improved prediction or better biological preservation and can
amplify features with very small training variability.
}
\value{
\code{fit_feature_scaler} returns a \code{feature_scaler} containing
\code{method}, \code{feature_names}, named \code{center} and \code{scale}
vectors, \code{training_samples}, and \code{zero_variance_features}.
The latter lists constant features detected during z-score fitting. For
identity scaling, center and scale are zero and one and this list is empty.
\code{scale_features} returns a samples-by-features numeric matrix with
fitted feature order and preserved sample row names.
}
\examples{
train <- rbind(s1 = c(A = 1, B = 4), s2 = c(A = 3, B = 4),
               s3 = c(A = 5, B = 4))
scaler <- fit_feature_scaler(train, "zscore")
scale_features(scaler, train)
scaler$zero_variance_features
scale_features(scaler, c(B = 5, A = 7))
}
\seealso{\code{\link{fit_cross_domain}}, \code{\link{fit_feature_encoder}}}
