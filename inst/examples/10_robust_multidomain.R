library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

# Robust centers are optional; their equality does not imply equal arithmetic means.
robust <- fit_domain_method(tr$x, tr$domain, "project", lambda = .4,
                            center = "trimmed", trim = .1, anchor = "sample_mean")
print(robust)
print(diagnose_cross_domain(robust, tr$x, tr$domain)$summary)
# Limit projection to the leading estimated domain direction.
leading <- fit_domain_method(tr$x, tr$domain, "project", lambda = .4, n_components = 1)
stopifnot(ncol(leading$adjuster$basis) == 1)
