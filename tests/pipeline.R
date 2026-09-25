library(CrossDomainAdjust)

near <- function(a, b, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance)))
}
must_error <- function(expr, pattern = NULL) {
  got <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(got, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(got)))
}

set.seed(731)
x <- cbind(g1 = c(rep(10, 8L), rep(-10, 16L)), g2 = 0,
           g3 = rep(c(-1, -1, 1, 1), 6L), g4 = rnorm(24L))
rownames(x) <- paste0("sample", seq_len(nrow(x)))
domain <- rep(c("A", "B"), 12L)
folds <- rep(1:3, each = 8L)
y <- seq_len(nrow(x)) %% 2L
grid <- c(0, 0.5, 1)

# Independent reference implementation: enumerate pairs on each training fold,
# apply inclusive prevalence bounds, then project onto its two-domain contrast.
ij <- utils::combn(seq_len(ncol(x)), 2L)
all_pairs <- 1L * (x[, ij[1L, ], drop = FALSE] > x[, ij[2L, ], drop = FALSE])
manual_fold <- function(fold, lambda) {
  valid <- folds == fold
  prevalence <- colMeans(all_pairs[!valid, , drop = FALSE])
  keep <- prevalence >= 0.2 & prevalence <= 0.8
  z <- all_pairs[, keep, drop = FALSE]
  ma <- colMeans(z[!valid & domain == "A", , drop = FALSE])
  mb <- colMeans(z[!valid & domain == "B", , drop = FALSE])
  a <- (ma + mb) / 2
  b <- ma - mb
  if (sum(b^2) > 0) {
    b <- b / sqrt(sum(b^2))
    z <- z - lambda * (sweep(z, 2L, a) %*% b) %*% t(b)
  }
  list(train = z[!valid, , drop = FALSE], valid = z[valid, , drop = FALSE],
       selected = which(keep))
}
feature_counts <- vapply(1:3, function(f) ncol(manual_fold(f, 0)$train), integer(1))
stopifnot(length(unique(feature_counts)) > 1L)
calls <- 0L
score_calls <- 0L
fit_callback <- function(train, outcome) {
  calls <<- calls + 1L
  fold <- (calls - 1L) %/% length(grid) + 1L
  lambda <- grid[(calls - 1L) %% length(grid) + 1L]
  expected <- manual_fold(fold, lambda)
  near(train, expected$train)
  stopifnot(identical(rownames(train), rownames(expected$train)))
  near(outcome, y[folds != fold])
  list(fold = fold, lambda = lambda, expected = expected$valid)
}
predict_callback <- function(model, valid) {
  near(valid, model$expected)
  stopifnot(identical(rownames(valid), rownames(model$expected)))
  rep(model$lambda, nrow(valid))
}
score_callback <- function(outcome, prediction, domains) {
  score_calls <<- score_calls + 1L
  fold <- (score_calls - 1L) %/% length(grid) + 1L
  near(outcome, y[folds == fold])
  stopifnot(identical(domains, domain[folds == fold]))
  # Deliberately maximize at 0.5 to test the direction of grid selection.
  1 - abs(prediction[1L] - 0.5)
}
tuned <- tune_cross_domain(x, domain, y, folds, fit_callback, predict_callback,
  score_callback, lambda_grid = c(1, 0.5, 0, 0.5), representation = "pair",
  correction = "project", frequency = c(0.2, 0.8), chunk_size = 2)
stopifnot(calls == 9L, score_calls == 9L, tuned$best_lambda == 0.5,
          nrow(tuned$fold_scores) == 9L,
          identical(tuned$scores$lambda, grid),
          all(tuned$fold_scores$features == rep(feature_counts, each = length(grid))),
          grepl("untouched outer test", tuned$note))
full <- fit_cross_domain(x, domain, "pair", "project", lambda = 0.5, chunk_size = 2)
near(transform_cross_domain(tuned$fit, x), transform_cross_domain(full, x))
stopifnot(tuned$fit$training_samples == nrow(x),
          tuned$fit$encoder$training_samples == nrow(x))
# A held-out sample can be transformed using the fitted training schema.
outer <- x[1:2, , drop = FALSE] + 0.37
names_before <- tuned$fit$encoder$output_names
invisible(transform_cross_domain(tuned$fit, outer, c("external", "external")))
stopifnot(identical(names_before, tuned$fit$encoder$output_names))

for (repr in c("raw", "rank", "pair")) {
  for (correction in c("none", "center", "project", "hybrid")) {
    z <- fit_cross_domain(x, domain, repr, correction, lambda = 0)
    near(transform_cross_domain(z, x), encode_features(z$encoder, x))
    stopifnot(identical(rownames(transform_cross_domain(z, x)), rownames(x)))
  }
}
hybrid <- fit_cross_domain(x, domain, "raw", "hybrid", lambda = c(0.2, 0.7))
near(transform_cross_domain(hybrid, x, domain),
     adjust_features(hybrid$adjuster, x, domain, c(0.2, 0.7)))
diagnostic <- diagnose_cross_domain(full, x, domain)
stopifnot(nrow(diagnostic$summary) == 2L,
          length(diagnostic$sample_displacement) == nrow(x),
          identical(names(diagnostic$sample_displacement), rownames(x)),
          grepl("not a measure", diagnostic$note))
one <- matrix(c(1, 3, 7, 9), ncol = 1L, dimnames = list(NULL, "only"))
one_domain <- c("A", "A", "B", "B")
one_fit <- fit_cross_domain(one, one_domain, "raw", "project")
one_diagnostic <- diagnose_cross_domain(one_fit, one, one_domain)
near(one_diagnostic$summary$mean_centroid_distance_per_sqrt_feature, c(6, 0))
single_pair <- fit_cross_domain(x, domain, "pair", "project",
  pairs = matrix(c("g2", "g3"), nrow = 1L))
stopifnot(ncol(transform_cross_domain(single_pair, x)) == 1L)
invisible(diagnose_cross_domain(single_pair, x, domain))

# Leave-one-domain validation is valid for projection with two remaining domains.
# Frozen centering must reject the unknown validation center rather than infer it.
three_domain <- rep(c("A", "B", "C"), each = 8L)
train_known <- fit_cross_domain(x[1:16, , drop = FALSE], three_domain[1:16],
                               correction = "project")
stopifnot(nrow(transform_cross_domain(train_known, x[17:24, , drop = FALSE],
                                      three_domain[17:24])) == 8L)
train_center <- fit_cross_domain(x[1:16, , drop = FALSE], three_domain[1:16],
                                correction = "center")
must_error(transform_cross_domain(train_center, x[17:24, , drop = FALSE],
                                 three_domain[17:24]), "Unknown domain")

# Tie resolution chooses the lowest strength; scores use larger-is-better.
simple_fit <- function(z, outcome) mean(outcome)
simple_predict <- function(m, z) rep(m, nrow(z))
simple_score <- function(outcome, prediction, domains) 1
ties <- tune_cross_domain(x, domain, y, folds, simple_fit, simple_predict,
                         simple_score, lambda_grid = c(1, 0.25, 0), correction = "none")
stopifnot(ties$best_lambda == 0)
must_error(tune_cross_domain(x, domain, y, rep(1, 24), simple_fit, simple_predict,
                            simple_score), "at least two")
must_error(tune_cross_domain(x, domain, y, folds, simple_fit, simple_predict,
                            simple_score, lambda_grid = c(0, 2)), "lambda_grid")
must_error(tune_cross_domain(x, domain, y, folds, simple_fit, simple_predict,
                            function(a, b, c) NA_real_), "finite number")
must_error(transform_cross_domain(NULL, x), "fitted cross_domain_fit")
