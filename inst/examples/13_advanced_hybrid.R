library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

# Advanced composition is separate from the four study presets.
fit <- fit_cross_domain(tr$x, tr$domain, representation = "raw",
  correction = "hybrid", lambda = c(.2, .4), standardize = "none")
z <- transform_cross_domain(fit, te$x, te$domain)
stopifnot(all(is.finite(z)))
print(diagnose_cross_domain(fit, tr$x, tr$domain)$summary)
# The operations overlap; a hybrid is not assumed superior.
