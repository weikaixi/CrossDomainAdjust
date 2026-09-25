library(CrossDomainAdjust)

expect_error <- function(expr, pattern) {
  caught <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(caught, "error"),
            grepl(pattern, conditionMessage(caught), fixed = TRUE))
}

# Chunking and the streaming cap must agree with a direct exhaustive reference.
set.seed(7401)
x <- matrix(sample(-2:4, 35 * 9, replace = TRUE), 35, 9,
            dimnames = list(paste0("s", 1:35), paste0("g", 1:9)))
all_pairs <- combn(ncol(x), 2)
brute <- vapply(seq_len(ncol(all_pairs)), function(j) {
  as.integer(x[, all_pairs[1, j]] > x[, all_pairs[2, j]])
}, integer(nrow(x)))
f <- colMeans(brute)
qualify <- which(f >= 0.2 & f <= 0.8)
enc1 <- fit_feature_encoder(x, "pair", chunk_size = 1)
enc7 <- fit_feature_encoder(x, "pair", chunk_size = 7)
stopifnot(identical(enc1$pairs, enc7$pairs),
          identical(enc1$frequencies, enc7$frequencies),
          identical(encode_features(enc1, x), encode_features(enc7, x)),
          identical(unname(encode_features(enc1, x)),
                    unname(brute[, qualify, drop = FALSE])),
          identical(unname(enc1$frequencies), unname(f[qualify])),
          enc1$candidate_count == choose(ncol(x), 2))
best <- qualify[order(abs(f[qualify] - 0.5), qualify)][1:5]
capped <- fit_feature_encoder(x, "pair", chunk_size = 2, max_pairs = 5)
stopifnot(identical(as.integer(capped$candidate_indices), sort(best)),
          capped$retained_count == 5,
          capped$filtered_count == capped$candidate_count - 5,
          capped$capped_count == length(qualify) - 5,
          identical(capped$pairs,
                    fit_feature_encoder(x, "pair", chunk_size = 99,
                                        max_pairs = 5)$pairs))

# Exact inclusive endpoints, ties as zero, and frozen selection on new samples.
boundary <- cbind(A = c(1, 0, 0, 0, 0),
                  C = c(1, 1, 1, 1, 0), B = rep(0, 5))
bound_enc <- fit_feature_encoder(boundary, "pair")
stopifnot(identical(unname(bound_enc$frequencies), c(0.2, 0.8)),
          bound_enc$retained_count == 2,
          bound_enc$frequency_filtered_count == 1)
before <- bound_enc
new <- rbind(c(A = 1, C = 1, B = 1), c(A = 2, C = 2, B = 0))
encoded_new <- encode_features(bound_enc, new)
stopifnot(identical(bound_enc, before),
          all(encoded_new[1, ] == 0), all(encoded_new[2, ] == 1),
          identical(colnames(encoded_new), bound_enc$output_names))
constant <- matrix(1, 5, 3, dimnames = list(NULL, c("A", "C", "B")))
expect_error(fit_feature_encoder(constant, "pair"), "No pair features")
stopifnot(all(encode_features(
  fit_feature_encoder(constant, "pair", frequency = c(0, 1)), constant) == 0))

# Pair orientations, explicit mappings, name alignment and one-sample inputs.
explicit <- rbind(c("B", "A"), c("C", "B"))
ex <- fit_feature_encoder(boundary, "pair", pairs = explicit)
stopifnot(identical(unname(ex$pairs), rbind(c("A", "B"), c("C", "B"))),
          identical(encode_features(ex, boundary),
                    encode_features(ex, cbind(extra = 100,
                                              boundary[, c("B", "A", "C")]))),
          identical(unname(encode_features(ex, new[1, ])),
                    unname(encoded_new[1, , drop = FALSE])))
expect_error(fit_feature_encoder(boundary, "pair",
                                pairs = rbind(c("A", "B"), c("B", "A"))),
             "Duplicate or reversed")
expect_error(fit_feature_encoder(boundary, "pair", pairs = rbind(c("A", "A"))),
             "Self pairs")
expect_error(fit_feature_encoder(boundary, "pair", pairs = rbind(c("A", "Z"))),
             "must occur")
expect_error(encode_features(ex, boundary[, c("A", "B")]), "Missing fitted input")

# Raw/rank preserve genes, frozen universe, average ties and sample names.
rank_x <- rbind(s1 = c(A = 5, B = 5, C = 1),
                s2 = c(A = 2, B = 2, C = 2))
rank_enc <- fit_feature_encoder(rank_x, "rank")
rank_out <- encode_features(rank_enc, rank_x)
stopifnot(identical(unname(rank_out), rbind(c(0.75, 0.75, 0), rep(0.5, 3))),
          identical(dimnames(rank_out), dimnames(rank_x)),
          identical(rank_out,
                    encode_features(rank_enc, cbind(Z = 999, rank_x[, 3:1]))),
          identical(encode_features(fit_feature_encoder(x, "raw"), x[, 9:1]), x))

# Validation and the quadratic guard fail early and informatively.
wide <- matrix(0, 1, 2000, dimnames = list(NULL, paste0("g", 1:2000)))
expect_error(fit_feature_encoder(wide, "pair"), "exceeds max_candidates")
expect_error(fit_feature_encoder(boundary, "pair", max_pairs = 0), "max_pairs")
expect_error(fit_feature_encoder(boundary, "pair", chunk_size = 1.5), "chunk_size")
expect_error(fit_feature_encoder(boundary, "pair", frequency = c(0.8, 0.2)),
             "frequency")
bad <- boundary; colnames(bad)[2] <- "A"
expect_error(fit_feature_encoder(bad), "unique, nonempty")
bad <- boundary; bad[1, 1] <- NA_real_
expect_error(fit_feature_encoder(bad), "finite values")
bad[1, 1] <- Inf
expect_error(encode_features(bound_enc, bad), "finite values")
expect_error(fit_feature_encoder(unname(boundary)), "unique, nonempty")
expect_error(fit_feature_encoder(matrix(1, 2, 1,
                                      dimnames = list(NULL, "A")), "rank"),
             "at least two genes")
