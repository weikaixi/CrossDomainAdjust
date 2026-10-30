\name{encode_features}
\alias{encode_features}
\title{Apply a Frozen Feature Encoder}
\description{
Encodes training or new samples with an encoder fitted by
\code{\link{fit_feature_encoder}}, preserving its fitted feature universe,
pair orientations, frequency selection, and output order.
}
\usage{
encode_features(encoder, x)
}
\arguments{
  \item{encoder}{An object returned by \code{\link{fit_feature_encoder}}.}
  \item{x}{Finite numeric matrix or numeric data frame, with samples in rows
  and uniquely named genes in columns, or a named numeric vector for one
  sample. All fitted input genes must be present, including fitted genes not
  used in retained pairs. Extra columns are ignored and fitted genes are
  matched by name, so input-column order may differ.}
}
\details{
No frequency threshold is re-estimated from these samples. A fitted pair can
therefore be constant in new data and remains in the output. Pairs use strict
greater-than comparison, with ties encoded as zero. Rank representations use
average within-sample ranks over the complete fitted gene universe and map
them to \eqn{[0,1]} by \eqn{(r-1)/(p-1)}. The returned matrix is dense: a
cap in \code{fit_feature_encoder} can be useful when retaining many pairs.
}
\value{
A numeric matrix with one row per sample and columns in
\code{encoder$output_names} order. Input sample row names are preserved.
Pair matrices contain integer zeros and ones.
}
\examples{
train <- rbind(s1 = c(A = 3, B = 1, C = 2),
               s2 = c(A = 1, B = 3, C = 2))
encoder <- fit_feature_encoder(train, "pair")
encode_features(encoder, c(C = 4, B = 2, A = 1))
rank_encoder <- fit_feature_encoder(train, "rank")
encode_features(rank_encoder, train)
}
\seealso{\code{\link{fit_feature_encoder}}}
