library(CrossDomainAdjust)
near <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-10)))
set.seed(111)
rng <- .Random.seed
d <- simulate_domain_data(n_train = 90, n_test = 30, p = 8)
stopifnot(identical(rng, .Random.seed), identical(d, simulate_domain_data(90, 30, 8)))
x <- d$train$x; g <- d$train$domain
for (method in c("rank", "pair", "center", "project")) {
  a <- fit_domain_method(x, g, method, lambda = .4)
  b <- fit_cross_domain(x, g,
    representation = if (method %in% c("rank", "pair")) method else "raw",
    correction = if (method %in% c("rank", "pair")) "none" else method,
    standardize = if (method %in% c("rank", "pair")) "zscore" else "none", lambda = .4)
  near(transform_cross_domain(a, d$test$x, d$test$domain),
       transform_cross_domain(b, d$test$x, d$test$domain))
}
for (method in c("center", "project")) {
  a <- fit_domain_method(x, g, method, lambda = .4)
  saved <- serialize(a, NULL)
  diag <- lambda_diagnostics(a, x, g)
  stopifnot(identical(saved, serialize(a, NULL)), diag$summary$domain_r2[6] < 1e-20)
  near(diag$summary$mean_centroid_distance / diag$summary$mean_centroid_distance[1],
       1 - diag$summary$lambda)
  near(transform_cross_domain(a, x, g, lambda = 0), x)
  near(unname(as.matrix(diag$coordinates[diag$coordinates$lambda == .4, c("PC1", "PC2")])),
       unname(sweep(transform_cross_domain(a, x, g), 2, diag$pca_center) %*% diag$pca_rotation))
}
# Structured outcomes must retain row identity and class across fold-local fits.
folds <- rep(1:3, each = 30)
for (y in list(data.frame(id = seq_len(nrow(x)), event = d$train$survival$event),
               cbind(id = seq_len(nrow(x)), event = d$train$survival$event))) {
  calls <- 0
  tuned <- tune_cross_domain(x, g, y, folds,
    fit_model = function(z, yy) {
      calls <<- calls + 1
      fold <- (calls - 1) %/% 6 + 1
      stopifnot(identical(yy, y[folds != fold, , drop = FALSE]))
      0
    }, predict_model = function(m, z) rep(m, nrow(z)),
    score = function(yy, pred, domain) {
      fold <- (calls - 1) %/% 6 + 1
      stopifnot(identical(yy, y[folds == fold, , drop = FALSE]))
      1
    })
  stopifnot(calls == 18, tuned$best_lambda == 0)
}
if (requireNamespace("survival", quietly = TRUE)) {
  y <- survival::Surv(d$train$survival$time, d$train$survival$event)
  tune_cross_domain(x, g, y, folds,
    fit_model = function(z, yy) { stopifnot(inherits(yy, "Surv")); 0 },
    predict_model = function(m, z) rep(m, nrow(z)),
    score = function(yy, pred, domain) { stopifnot(inherits(yy, "Surv")); 1 })
}
