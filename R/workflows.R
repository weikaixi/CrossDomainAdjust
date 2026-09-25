# Convenient presets retain the general pipeline as the single implementation.
fit_domain_method <- function(x, domain, method = c("rank", "pair", "center", "project"),
                              lambda = 1, ...) {
  method <- match.arg(method)
  dots <- list(...)
  if (any(names(dots) %in% c("representation", "correction", "standardize")))
    stop("Presets fix representation, correction and standardize; use fit_cross_domain for a custom pipeline.", call. = FALSE)
  encoding <- method %in% c("rank", "pair")
  do.call(fit_cross_domain, c(list(x = x, domain = domain,
    representation = if (encoding) method else "raw",
    correction = if (encoding) "none" else method,
    standardize = if (encoding) "zscore" else "none",
    lambda = if (encoding) 0 else lambda), dots))
}

lambda_diagnostics <- function(object, x, domain,
                               lambda_grid = c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
  if (!inherits(object, "cross_domain_fit") ||
      !object$correction %in% c("center", "project"))
    stop("object must be a center or project cross_domain_fit.", call. = FALSE)
  if (!is.numeric(lambda_grid) || !length(lambda_grid) || any(!is.finite(lambda_grid)) ||
      any(lambda_grid < 0 | lambda_grid > 1)) stop("lambda_grid must be finite and in [0,1].", call. = FALSE)
  before <- transform_cross_domain(object, x, domain, lambda = 0)
  domain <- cd_domain(domain, nrow(before))
  if (length(unique(domain)) < 2L || nrow(before) < 2L || ncol(before) < 2L)
    stop("Diagnostics require at least two domains, samples and features.", call. = FALSE)
  pca <- stats::prcomp(before, center = TRUE, scale. = FALSE, rank. = 2)
  if (sum(pca$sdev^2) == 0) stop("PCA requires nonzero total variation.", call. = FALSE)
  ids <- rownames(before)
  if (is.null(ids)) ids <- paste0("sample_", seq_len(nrow(before)))
  summary <- coordinates <- list()
  for (strength in sort(unique(lambda_grid))) {
    z <- transform_cross_domain(object, x, domain, lambda = strength)
    lev <- unique(domain)
    centers <- do.call(rbind, lapply(lev, function(d) colMeans(z[domain == d, , drop = FALSE])))
    sizes <- vapply(lev, function(d) sum(domain == d), numeric(1))
    delta <- sweep(centers, 2, colMeans(z), "-")
    between <- sum(sizes * rowSums(delta^2))
    total <- sum(sweep(z, 2, colMeans(z), "-")^2)
    xy <- sweep(z, 2, pca$center, "-") %*% pca$rotation[, 1:2, drop = FALSE]
    summary[[length(summary) + 1L]] <- data.frame(lambda = strength,
      domain_r2 = if (total > 0) min(1, max(0, between / total)) else NA_real_,
      mean_centroid_distance = mean(as.numeric(stats::dist(centers))),
      mean_displacement = mean(sqrt(rowSums((z - before)^2))))
    coordinates[[length(coordinates) + 1L]] <- data.frame(sample = ids, domain = domain,
      lambda = strength, PC1 = xy[, 1], PC2 = xy[, 2], row.names = NULL)
  }
  list(summary = do.call(rbind, summary), coordinates = do.call(rbind, coordinates),
    pca_center = pca$center, pca_rotation = pca$rotation[, 1:2, drop = FALSE],
    baseline_variance_explained = pca$sdev[1:2]^2 / sum(pca$sdev^2),
    note = "Descriptive fitting-data geometry, not lambda selection. Equal centroids do not imply identical distributions or preserved prognosis.")
}

simulate_domain_data <- function(n_train = 180, n_test = 90, p = 20, domains = 3, seed = 42) {
  for (nm in c("n_train", "n_test", "p", "domains")) fe_positive_integer(get(nm), nm)
  if (p < 4 || domains < 2 || min(n_train, n_test) < domains)
    stop("Use p >= 4, domains >= 2, and at least one sample per domain in each split.", call. = FALSE)
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max)
    stop("seed must be a finite nonnegative integer in R's seed range.", call. = FALSE)
  if (seed != floor(seed)) stop("seed must be an integer.", call. = FALSE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  shifts <- matrix(stats::rnorm(domains * p, sd = 0.7), domains, p)
  lev <- paste0("D", seq_len(domains))
  draw <- function(n, prefix) {
    domain <- rep(lev, length.out = n)
    latent <- matrix(stats::rnorm(n * p), n, p)
    biology <- (match(domain, lev) - 1) / domains
    latent[, 1] <- latent[, 1] + biology
    eta <- 0.55 * latent[, 1] - 0.45 * latent[, 2] + 0.25 * latent[, 3]
    x <- latent + shifts[match(domain, lev), , drop = FALSE]
    dimnames(x) <- list(paste0(prefix, seq_len(n)), paste0("gene", seq_len(p)))
    death <- stats::rexp(n, rate = exp(eta) / 24)
    censor <- stats::rexp(n, rate = 1 / 40)
    list(x = x, domain = domain, regression = eta + stats::rnorm(n, sd = .5),
      binary = stats::rbinom(n, 1, stats::plogis(eta)),
      survival = data.frame(time = pmin(death, censor), event = as.integer(death <= censor)))
  }
  list(train = draw(n_train, "train"), test = draw(n_test, "test"),
    description = "Synthetic educational data with shared technical offsets and domain-associated biology; not a clinical dataset or benchmark of superiority.")
}
