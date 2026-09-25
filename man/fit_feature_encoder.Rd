\name{fit_feature_encoder}
\alias{fit_feature_encoder}
\title{Fit a Frozen Raw, Rank, or Gene-Pair Feature Encoder}
\description{
Fits a feature representation using training samples only. Gene-pair encoders
retain comparisons whose observed frequency of one lies within specified
inclusive bounds, optionally capping the retained feature count.
}
\usage{
fit_feature_encoder(x, representation = c("raw", "rank", "pair"),
  frequency = c(0.2, 0.8), pairs = NULL, max_pairs = NULL,
  max_candidates = 1000000, chunk_size = 1000)
}
\arguments{
  \item{x}{A finite numeric matrix or numeric data frame with samples in rows
  and genes in columns. Columns must have unique, nonempty names. A named
  numeric vector is accepted as a single sample.}
  \item{representation}{Raw values, normalized within-sample ranks, or binary
  within-sample gene-pair comparisons. Rank and pair representations require
  at least two input genes.}
  \item{frequency}{Two ordered numbers in \eqn{[0,1]}. For pairs, the lower and
  upper inclusive training-sample frequencies of one. The default retains
  frequencies from 20 to 80 percent. Unused by raw and rank encoders.}
  \item{pairs}{Optional two-column character matrix or data frame of gene
  names. If omitted, all unordered pairs are candidates. Explicit pairs are
  oriented according to the fitted input-column order; self pairs and duplicate
  or reversed redundant pairs are rejected. Used only for pair representation.}
  \item{max_pairs}{Optional positive integer cap applied after frequency
  filtering. Pairs closest to frequency 0.5 are preferred; ties are resolved
  by candidate order. This is an unsupervised balance criterion, not a measure
  of biological relevance or predictive importance.}
  \item{max_candidates}{Positive integer computational guard. Fitting stops
  before enumerating all pairs if their count exceeds this limit. Preselect
  genes or supply explicit pairs to reduce the candidate search.}
  \item{chunk_size}{Maximum number of comparisons processed in one block.}
}
\details{
For genes \eqn{A} and \eqn{B}, a pair is \eqn{I(x_A > x_B)}; ties produce zero.
The earlier training column is always \eqn{A}. With all pairs, candidates are
enumerated by the first column, then the second column. Explicit candidates
follow the supplied row order after orientation. Output names such as
\code{pair_000001} record candidate indices and may have gaps after filtering.
The pair-to-gene mapping is stored in the encoder.

Frequencies are estimated over all supplied training samples. Frequency
filtering controls near-constant features but does not eliminate redundant
pairs, guarantee a particular retained count, or select outcome-associated
pairs. Use \code{max_pairs} for an additional count cap. Filtering is fitted
anew within each training fold during cross-validation, never on a held-out
test cohort. Since a small cohort gives discrete frequencies, bounds may
retain no pairs; in that case fitting stops with an informative error.

Candidate frequencies are computed in blocks, without materializing a
samples-by-all-candidates matrix or a list of all possible pairs. The search
still requires time proportional to the number of samples times the number
of candidate pairs. Generating the retained matrix with
\code{\link{encode_features}} requires memory proportional to samples times
retained features.

For ranks, each sample is ranked within the entire fitted gene universe using
average ranks for ties, then transformed as \eqn{(r-1)/(p-1)}, where \eqn{p}
is the number of fitted genes. No cross-sample normalization is performed.
Raw encoding preserves the fitted column order.
}
\value{
An object of class \code{feature_encoder} containing \code{representation},
the input \code{feature_names}, \code{output_names}, and the stored
\code{frequency} bounds. For pairs, \code{pairs} is a two-column character
matrix with columns \code{first} and \code{second} and output feature names
as row names; \code{frequencies} is the corresponding named numeric vector.
\code{candidate_count}, \code{retained_count}, and \code{filtered_count}
report the total candidate, retained, and excluded counts. The excluded count
includes both frequency filtering and the optional cap;
\code{frequency_filtered_count} and \code{capped_count} distinguish these.
\code{frequency_pass_count} is the count before the cap. For raw and rank
encoders, counts refer to input genes and \code{pairs} and \code{frequencies}
are \code{NULL}.
}
\examples{
x <- cbind(A = c(1, 0, 0, 0, 0),
           C = c(1, 1, 1, 1, 0), B = rep(0, 5))
encoder <- fit_feature_encoder(x, "pair", frequency = c(0.2, 0.8))
encoder$pairs
encoder$frequencies
encode_features(encoder, x)
# The fitted comparison set is unchanged for a new individual.
encode_features(encoder, c(A = 2, C = 0, B = 1))
}
\seealso{\code{\link{encode_features}}}
