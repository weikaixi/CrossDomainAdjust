library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

fit <- fit_domain_method(tr$x, tr$domain, "center", lambda = .6)
adjusted <- transform_cross_domain(fit, te$x, te$domain)
stopifnot(identical(dim(adjusted), dim(te$x)))
# Domain means came only from training data; no target-cohort mean is used.
print(diagnose_cross_domain(fit, tr$x, tr$domain)$summary)
unknown <- tryCatch(transform_cross_domain(fit, te$x[1, ], "new_domain"),
                    error = function(e) conditionMessage(e))
cat("Expected validation message:", unknown, "
")
stopifnot(is.character(unknown))
