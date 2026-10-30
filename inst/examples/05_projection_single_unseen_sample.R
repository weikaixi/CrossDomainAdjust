library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

fit <- fit_domain_method(tr$x, tr$domain, "project", lambda = .6)
# No target label or target-cohort refit is required for projection.
single <- transform_cross_domain(fit, te$x[1, ])
batch <- transform_cross_domain(fit, te$x, rep("new_site", nrow(te$x)))
stopifnot(isTRUE(all.equal(as.numeric(single), as.numeric(batch[1, ]))))
print(single)
# This is deployment feasibility, not evidence that an unseen shift is removed.
