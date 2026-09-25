library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

# Raw input receives no additional z-score scaling.
fits <- lapply(c("rank", "pair", "center", "project"), function(method)
  fit_domain_method(tr$x, tr$domain, method, lambda = .4, max_pairs = 20))
names(fits) <- c("rank", "pair", "center", "project")
for (method in names(fits)) {
  z <- transform_cross_domain(fits[[method]], te$x, te$domain)
  cat(method, ":", nrow(z), "samples x", ncol(z), "features
")
}
