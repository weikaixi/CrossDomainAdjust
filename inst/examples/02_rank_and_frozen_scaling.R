library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

fit <- fit_domain_method(tr$x, tr$domain, "rank")
# Strictly increasing transformations within a sample preserve its ordering.
a <- transform_cross_domain(fit, te$x)
b <- transform_cross_domain(fit, exp(te$x / 10))
stopifnot(isTRUE(all.equal(a, b)))
# The single-sample result equals the same row transformed in a whole cohort.
stopifnot(isTRUE(all.equal(as.numeric(a[1, ]),
                          as.numeric(transform_cross_domain(fit, te$x[1, ])))))
# Fractional ranks are [0,1]; subsequent training-fitted z-scores need not be.
print(range(encode_features(fit$encoder, te$x)))
print(range(a))
print(head(fit$scaler))
