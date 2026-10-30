# Feature-wise scaling is fitted once on training samples. It is distinct from
# the within-sample rank representation implemented in features.R.
fit_feature_scaler <- function(x, method = c("none", "zscore")) {
  x <- cd_matrix(x)
  method <- match.arg(method)
  feature_names <- colnames(x)
  center <- stats::setNames(rep(0, ncol(x)), feature_names)
  scale <- stats::setNames(rep(1, ncol(x)), feature_names)
  zero_variance_features <- character()
  if (method == "zscore") {
    if (nrow(x) < 2L) {
      stop("zscore fitting requires at least two training samples to estimate sample SD.",
           call. = FALSE)
    }
    center <- colMeans(x)
    scale <- apply(x, 2L, stats::sd)
    if (any(!is.finite(center)) || any(!is.finite(scale))) {
      stop("Training means or standard deviations are nonfinite; rescale extreme input values first.",
           call. = FALSE)
    }
    constant <- scale == 0
    zero_variance_features <- feature_names[constant]
    scale[constant] <- 1
  }
  structure(list(method = method, feature_names = feature_names, center = center,
                 scale = scale, zero_variance_features = zero_variance_features,
                 training_samples = nrow(x)), class = "feature_scaler")
}

scale_features <- function(scaler, x) {
  if (!inherits(scaler, "feature_scaler")) {
    stop("scaler must be produced by fit_feature_scaler().", call. = FALSE)
  }
  x <- cd_matrix(x, scaler$feature_names)
  if (scaler$method == "none") return(x)
  if (!identical(scaler$method, "zscore")) {
    stop("Unrecognized method in scaler.", call. = FALSE)
  }
  result <- sweep(sweep(x, 2L, scaler$center, "-"), 2L, scaler$scale, "/")
  if (any(!is.finite(result))) {
    stop("Scaled feature values are nonfinite; input values exceed the fitted scale's numeric range.",
         call. = FALSE)
  }
  result
}
