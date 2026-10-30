fit_cross_domain <- function(x, domain, representation = c("raw", "rank", "pair"),
                             correction = c("none", "center", "project", "hybrid"),
                             lambda = 1, frequency = c(0.2, 0.8), pairs = NULL,
                             max_pairs = NULL, max_candidates = 1000000, chunk_size = 1000,
                             anchor = "domain_mean", center = "mean", trim = 0.1,
                             n_components = NULL, standardize = c("none", "zscore")) {
  representation <- match.arg(representation)
  correction <- match.arg(correction)
  standardize <- match.arg(standardize)
  lambda <- cd_lambda(lambda, correction)
  encoder <- fit_feature_encoder(x, representation, frequency, pairs,
                                 max_pairs, max_candidates, chunk_size)
  encoded <- encode_features(encoder, x)
  scaler <- fit_feature_scaler(encoded, standardize)
  scaled <- scale_features(scaler, encoded)
  adjuster <- fit_domain_adjuster(scaled, domain, correction, anchor, center, trim, n_components)
  structure(list(encoder = encoder, scaler = scaler, adjuster = adjuster, lambda = lambda,
    representation = representation, standardize = standardize, correction = correction,
    training_samples = nrow(encoded), version = "2.1.0"), class = "cross_domain_fit")
}

transform_cross_domain <- function(object, x, domain = NULL, lambda = NULL) {
  if (!inherits(object, "cross_domain_fit")) stop("object must be a fitted cross_domain_fit.", call. = FALSE)
  encoded <- encode_features(object$encoder, x)
  # Version 2.0.0 objects had no scaler and used identity scaling.
  scaled <- if (is.null(object$scaler)) encoded else scale_features(object$scaler, encoded)
  adjust_features(object$adjuster, scaled, domain,
                  if (is.null(lambda)) object$lambda else lambda)
}

print.cross_domain_fit <- function(x, ...) {
  scaling <- if (is.null(x$scaler)) "none" else x$scaler$method
  version <- if (is.null(x$version)) "2.0" else x$version
  cat("CrossDomainAdjust ", version, ": ", x$representation, " -> ",
      scaling, " scaling -> ", x$correction, "\n", sep = "")
  cat(" Lambda:", paste(x$lambda, collapse = ", "), " Training samples:", x$training_samples, "\n")
  cat(" Output features:", length(x$adjuster$feature_names), "\n")
  invisible(x)
}

tune_cross_domain <- function(x, domain, y, folds, fit_model, predict_model, score,
                              lambda_grid = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                              representation = "raw", correction = "project", ...) {
  x <- cd_matrix(x)
  domain <- cd_domain(domain, nrow(x))
  if (!(is.atomic(y) || is.data.frame(y)) || NROW(y) != nrow(x) || anyNA(y))
    stop("y must contain one nonmissing outcome row per sample.", call. = FALSE)
  subset_y <- function(at) if (is.null(dim(y))) y[at] else y[at, , drop = FALSE]
  if (length(folds) != nrow(x) || anyNA(folds) || length(unique(folds)) < 2L)
    stop("folds must assign each sample to one of at least two validation folds.", call. = FALSE)
  if (!all(vapply(list(fit_model, predict_model, score), is.function, logical(1))))
    stop("fit_model, predict_model and score must be functions.", call. = FALSE)
  if (!is.numeric(lambda_grid) || !length(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid < 0 | lambda_grid > 1))
    stop("lambda_grid must contain finite values in [0,1].", call. = FALSE)
  grid <- sort(unique(lambda_grid))
  ids <- unique(folds)
  dots <- list(...)
  if (any(names(dots) %in% c("lambda", "x", "domain", "representation", "correction")))
    stop("Do not repeat fixed arguments through ...", call. = FALSE)
  records <- vector("list", length(ids) * length(grid))
  k <- 0L
  for (fold in ids) {
    valid <- folds == fold
    local <- do.call(fit_cross_domain, c(list(x = x[!valid, , drop = FALSE], domain = domain[!valid],
      representation = representation, correction = correction, lambda = 0), dots))
    for (strength in grid) {
      local$lambda <- cd_lambda(strength, correction)
      train_x <- transform_cross_domain(local, x[!valid, , drop = FALSE], domain[!valid])
      valid_x <- transform_cross_domain(local, x[valid, , drop = FALSE], domain[valid])
      model <- fit_model(train_x, subset_y(!valid))
      prediction <- predict_model(model, valid_x)
      value <- score(subset_y(valid), prediction, domain[valid])
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value))
        stop("score must return one finite number; larger is better.", call. = FALSE)
      k <- k + 1L
      records[[k]] <- data.frame(fold = as.character(fold), lambda = strength, score = value,
                                 features = ncol(train_x))
    }
  }
  detail <- do.call(rbind, records)
  summary <- data.frame(lambda = grid,
    mean_score = vapply(grid, function(l) mean(detail$score[detail$lambda == l]), numeric(1)),
    sd_score = vapply(grid, function(l) stats::sd(detail$score[detail$lambda == l]), numeric(1)))
  best <- summary$lambda[which.max(summary$mean_score)]
  final <- do.call(fit_cross_domain, c(list(x = x, domain = domain, representation = representation,
    correction = correction, lambda = best), dots))
  list(best_lambda = best, scores = summary, fold_scores = detail, fit = final,
       note = "Selection scores are not an unbiased performance estimate. Evaluate on an untouched outer test set.")
}

diagnose_cross_domain <- function(object, x, domain) {
  if (!inherits(object, "cross_domain_fit")) stop("Expected a cross_domain_fit.", call. = FALSE)
  before <- encode_features(object$encoder, x)
  if (!is.null(object$scaler)) before <- scale_features(object$scaler, before)
  domain <- cd_domain(domain, nrow(before))
  after <- transform_cross_domain(object, x, domain)
  lev <- unique(domain)
  center_distance <- function(z) {
    centers <- do.call(rbind, lapply(lev, function(d) colMeans(z[domain == d, , drop = FALSE])))
    if (length(lev) < 2L) return(NA_real_)
    mean(as.numeric(stats::dist(centers)) / sqrt(ncol(z)))
  }
  within <- function(z) mean(vapply(lev, function(d) {
    zz <- z[domain == d, , drop = FALSE]
    mean(sweep(zz, 2L, colMeans(zz), "-")^2)
  }, numeric(1)))
  before_stage <- if (!is.null(object$scaler) && object$scaler$method == "zscore") "standardized" else "encoded"
  list(summary = data.frame(stage = c(before_stage, "adjusted"),
    mean_centroid_distance_per_sqrt_feature = c(center_distance(before), center_distance(after)),
    mean_within_domain_squared_deviation = c(within(before), within(after))),
    sample_displacement = setNames(sqrt(rowMeans((after - before)^2)), rownames(before)),
    note = "Descriptive domain separation is not a measure of biological validity or predictive performance.")
}
