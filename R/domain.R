cd_matrix <- function(x, features = NULL) {
  if (is.numeric(x) && is.null(dim(x))) {
    if (is.null(names(x))) stop("A sample vector must have feature names.", call. = FALSE)
    x <- matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
  }
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) stop("All columns must be numeric.", call. = FALSE)
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x))
    stop("x must be a nonempty numeric matrix, data frame, or named vector.", call. = FALSE)
  nm <- colnames(x)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop("Feature names must be present, nonempty and unique.", call. = FALSE)
  if (any(!is.finite(x))) stop("All feature values must be finite; impute using training data first.", call. = FALSE)
  if (!is.null(features)) {
    if (any(!features %in% nm)) stop("Missing required features: ", paste(setdiff(features, nm), collapse = ", "), call. = FALSE)
    x <- x[, features, drop = FALSE]
  }
  x
}

cd_domain <- function(domain, n) {
  if (length(domain) != n || anyNA(domain)) stop("domain must contain one label per sample.", call. = FALSE)
  domain <- as.character(domain)
  if (any(!nzchar(domain))) stop("domain labels must be nonempty.", call. = FALSE)
  domain
}

cd_lambda <- function(lambda, method) {
  if (!is.numeric(lambda) || any(!is.finite(lambda)) || any(lambda < 0 | lambda > 1) ||
      !(length(lambda) %in% if (method == "hybrid") c(1L, 2L) else 1L))
    stop("lambda must be in [0,1]; hybrid accepts one or two values (center, project).", call. = FALSE)
  if (method == "hybrid" && length(lambda) == 1L) lambda <- rep(lambda, 2L)
  unname(lambda)
}

fit_domain_adjuster <- function(x, domain, method = c("project", "center", "hybrid", "none"),
                                anchor = "domain_mean", center = c("mean", "median", "trimmed"),
                                trim = 0.1, n_components = NULL) {
  x <- cd_matrix(x)
  domain <- cd_domain(domain, nrow(x))
  method <- match.arg(method)
  center <- match.arg(center)
  lev <- unique(domain)
  if (method != "none" && length(lev) < 2L) stop("At least two training domains are required.", call. = FALSE)
  if (length(trim) != 1L || !is.numeric(trim) || !is.finite(trim) || trim < 0 || trim >= 0.5)
    stop("trim must be in [0,0.5).", call. = FALSE)
  center_fun <- switch(center, mean = function(v) mean(v), median = stats::median,
                       trimmed = function(v) mean(v, trim = trim))
  centers <- do.call(rbind, lapply(lev, function(d) apply(x[domain == d, , drop = FALSE], 2, center_fun)))
  dimnames(centers) <- list(lev, colnames(x))
  sizes <- setNames(vapply(lev, function(d) sum(domain == d), integer(1)), lev)
  if (is.character(anchor) && length(anchor) == 1L) {
    anchor_type <- match.arg(anchor, c("domain_mean", "sample_mean"))
    a <- if (anchor_type == "domain_mean") colMeans(centers) else
      colSums(centers * as.numeric(sizes / sum(sizes)))
  } else {
    if (!is.numeric(anchor) || is.null(names(anchor)) || anyDuplicated(names(anchor)) ||
        any(!colnames(x) %in% names(anchor))) stop("Custom anchor must be a uniquely named numeric vector covering all features.", call. = FALSE)
    a <- anchor[colnames(x)]
    if (any(!is.finite(a))) stop("Custom anchor must be finite.", call. = FALSE)
    anchor_type <- "custom"
  }
  names(a) <- colnames(x)
  shifts <- sweep(centers, 2L, a, "-")
  if (method %in% c("project", "hybrid")) {
    s <- svd(shifts, nu = 0L, nv = min(dim(shifts)))
    tol <- max(dim(shifts)) * .Machine$double.eps * max(s$d)
    available <- sum(s$d > tol)
    if (is.null(n_components)) n_components <- available
    if (!is.numeric(n_components) || length(n_components) != 1L || !is.finite(n_components) ||
        n_components < 0 || n_components != as.integer(n_components) || n_components > available)
      stop("n_components must be an integer between zero and the estimated rank.", call. = FALSE)
    basis <- s$v[, seq_len(n_components), drop = FALSE]
    singular_values <- s$d
  } else {
    available <- 0L
    basis <- matrix(numeric(), ncol(x), 0L)
    singular_values <- numeric()
  }
  rownames(basis) <- colnames(x)
  structure(list(method = method, feature_names = colnames(x), domain_levels = lev,
    domain_centers = centers, domain_sizes = sizes, projection_anchor = a,
    anchor = a, anchor_type = anchor_type, center = center, trim = trim,
    shift_matrix = shifts, basis = basis, rank = ncol(basis), available_rank = available,
    singular_values = singular_values), class = "domain_adjuster")
}

adjust_features <- function(adjuster, x, domain = NULL, lambda = 1) {
  if (!inherits(adjuster, "domain_adjuster")) stop("adjuster must be a fitted domain_adjuster.", call. = FALSE)
  x <- cd_matrix(x, adjuster$feature_names)
  sample_names <- rownames(x)
  lambda <- cd_lambda(lambda, adjuster$method)
  if (!is.null(domain)) domain <- cd_domain(domain, nrow(x))
  if (adjuster$method == "none" || all(lambda == 0)) return(x)
  if (adjuster$method %in% c("center", "hybrid") && lambda[1L] > 0) {
    if (is.null(domain)) stop("Center adjustment requires domain labels.", call. = FALSE)
    if (any(!domain %in% adjuster$domain_levels))
      stop("Unknown domain: frozen centering requires an estimated training-domain center; no target center is inferred.", call. = FALSE)
    x <- x - lambda[1L] * adjuster$shift_matrix[domain, , drop = FALSE]
  }
  if (adjuster$method %in% c("project", "hybrid") && adjuster$rank > 0) {
    strength <- if (adjuster$method == "hybrid") lambda[2L] else lambda[1L]
    offsets <- sweep(x, 2L, adjuster$projection_anchor, "-")
    # Use O(p*r) storage, never construct a dense p-by-p projection matrix.
    x <- x - strength * tcrossprod(offsets %*% adjuster$basis, adjuster$basis)
  }
  dimnames(x) <- list(sample_names, adjuster$feature_names)
  x
}

fit_domain_projector <- function(x, domain, anchor = "domain_mean") {
  out <- fit_domain_adjuster(x, domain, method = "project", anchor = anchor)
  class(out) <- c("domain_projector", "domain_adjuster")
  out
}

project_features <- function(projector, x, lambda = 1) {
  if (!inherits(projector, "domain_projector")) stop("Expected a domain_projector.", call. = FALSE)
  adjust_features(projector, x, lambda = lambda)
}

project_new_sample <- function(projector, sample, domain = NULL, lambda = 1) {
  x <- cd_matrix(sample)
  if (nrow(x) != 1L) stop("sample must contain exactly one row.", call. = FALSE)
  if (!is.null(domain)) domain <- cd_domain(domain, 1L)
  result <- project_features(projector, x, lambda)
  out <- setNames(as.numeric(result[1L, ]), colnames(result))
  attr(out, "domain") <- domain
  attr(out, "lambda") <- lambda
  out
}

print.domain_adjuster <- function(x, ...) {
  cat("Cross-domain adjuster:", x$method, "\n")
  cat(" Features:", length(x$feature_names), " Domains:", length(x$domain_levels), " Projection rank:", x$rank, "\n")
  cat(" Frozen", x$center, "centers;", x$anchor_type, "anchor\n")
  invisible(x)
}

print.domain_projector <- function(x, ...) print.domain_adjuster(x, ...)
