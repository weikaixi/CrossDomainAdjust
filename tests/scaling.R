library(CrossDomainAdjust)

sc_expect_error <- function(expr, pattern) {
  caught <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(caught, "error"),
            grepl(pattern, conditionMessage(caught), fixed = TRUE))
}
close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12))

# Feature means and sample SDs are estimated only from the training matrix.
x <- rbind(s1 = c(A = 1, B = 7, C = 3),
           s2 = c(A = 2, B = 7, C = 9),
           s3 = c(A = 5, B = 7, C = 4),
           s4 = c(A = 8, B = 7, C = 10))
scaler <- fit_feature_scaler(x, "zscore")
z <- scale_features(scaler, x)
stopifnot(close(colMeans(z), c(A = 0, B = 0, C = 0)),
          close(apply(z[, c("A", "C")], 2, sd), c(A = 1, C = 1)),
          scaler$scale["B"] == 1,
          identical(scaler$zero_variance_features, "B"),
          identical(dimnames(z), dimnames(x)))
before <- scaler
new <- c(C = 1000, B = 8, A = -100)
expected <- (new[colnames(x)] - colMeans(x)) / scaler$scale
got <- scale_features(scaler, new)
stopifnot(close(unname(got[1, ]), unname(expected)), identical(scaler, before),
          got[1, "B"] == 1,
          close(scale_features(scaler, x[, 3:1]), z),
          close(scale_features(scaler, cbind(extra = 99, x)), z),
          identical(scale_features(fit_feature_scaler(x), x), x))
sc_expect_error(fit_feature_scaler(x[1, , drop = FALSE], "zscore"), "at least two")
bad <- x; bad[1, 1] <- Inf
sc_expect_error(fit_feature_scaler(bad, "zscore"), "finite")
bad[1, 1] <- NA_real_
sc_expect_error(scale_features(scaler, bad), "finite")
sc_expect_error(scale_features(scaler, x[, 1:2]), "Missing required")
sc_expect_error(scale_features(list(), x), "fit_feature_scaler")

# Pipeline order: within-sample ranks, training feature scaling, then projection.
set.seed(220)
a <- matrix(rnorm(12 * 5), 12, 5,
            dimnames = list(paste0("sample", 1:12), paste0("g", 1:5)))
d <- rep(c("A", "B", "C"), each = 4)
encoder <- fit_feature_encoder(a, "rank")
ranked <- encode_features(encoder, a)
manual_scaler <- fit_feature_scaler(ranked, "zscore")
standardized <- scale_features(manual_scaler, ranked)
manual_adjuster <- fit_domain_adjuster(standardized, d, method = "project")
manual <- adjust_features(manual_adjuster, standardized, d, lambda = 0.4)
pipeline <- fit_cross_domain(a, d, representation = "rank", correction = "project",
                             lambda = 0.4, standardize = "zscore")
stopifnot(close(transform_cross_domain(pipeline, a, d), manual),
          identical(pipeline$scaler, manual_scaler),
          close(unname(transform_cross_domain(pipeline, a[1, ])),
                unname(manual[1, , drop = FALSE])))
diagnostics <- diagnose_cross_domain(pipeline, a, d)
stopifnot(diagnostics$summary$stage[1] == "standardized",
          close(diagnostics$sample_displacement,
                setNames(sqrt(rowMeans((manual - standardized)^2)), rownames(a))))

# Default scaling is numerically identical to the previous encoder-adjuster flow.
raw_adjuster <- fit_domain_adjuster(a, d, method = "project")
default <- fit_cross_domain(a, d, correction = "project", lambda = 0.4)
default_out <- transform_cross_domain(default, a, d)
stopifnot(identical(default_out, adjust_features(raw_adjuster, a, d, lambda = 0.4)))
legacy <- default
legacy$scaler <- NULL
legacy$standardize <- NULL
legacy$version <- "2.0.0"
stopifnot(identical(transform_cross_domain(legacy, a, d), default_out),
          identical(diagnose_cross_domain(legacy, a, d),
                    diagnose_cross_domain(default, a, d)))

# Fold-local fit and validation transforms are checked against independent
# training-only calculations. Large fold shifts would expose accidental refits.
b <- cbind(g1 = c(1:4, 21:24, 101:104),
           g2 = c(7, 5, 9, 8, 70, 50, 90, 80, 700, 500, 900, 800))
rownames(b) <- as.character(seq_len(nrow(b)))
folds <- rep(1:3, each = 4)
dd <- rep(c("A", "B"), 6)
seen <- 0L
tuned <- tune_cross_domain(b, dd, seq_len(nrow(b)), folds,
  fit_model = function(xx, yy) {
    mu <- colMeans(b[yy, , drop = FALSE])
    ss <- apply(b[yy, , drop = FALSE], 2, sd)
    expected <- sweep(sweep(b[yy, , drop = FALSE], 2, mu, "-"), 2, ss, "/")
    stopifnot(close(xx, expected), max(abs(colMeans(xx))) < 1e-12)
    list(mu = mu, sd = ss)
  },
  predict_model = function(model, xx) {
    ids <- as.integer(rownames(xx))
    expected <- sweep(sweep(b[ids, , drop = FALSE], 2, model$mu, "-"), 2, model$sd, "/")
    stopifnot(close(xx, expected))
    seen <<- seen + 1L
    xx[, 1]
  },
  score = function(y, prediction, domain) -mean((prediction - y)^2),
  lambda_grid = 0, representation = "raw", correction = "none", standardize = "zscore")
stopifnot(seen == 3L, close(tuned$fit$scaler$center, colMeans(b)),
          close(tuned$fit$scaler$scale, apply(b, 2, sd)))
