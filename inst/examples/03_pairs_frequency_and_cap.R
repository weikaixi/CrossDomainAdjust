library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

fit <- fit_domain_method(tr$x, tr$domain, "pair", frequency = c(.2, .8),
                         max_pairs = 12, chunk_size = 5)
print(fit$encoder$pairs)
print(fit$encoder$frequencies)
stopifnot(all(fit$encoder$frequencies >= .2 & fit$encoder$frequencies <= .8),
          ncol(transform_cross_domain(fit, te$x)) <= 12)
# Explicit candidates avoid enumerating all p*(p-1)/2 possible pairs.
candidate <- rbind(c("gene1", "gene2"), c("gene3", "gene4"), c("gene5", "gene6"))
selected <- fit_domain_method(tr$x, tr$domain, "pair", pairs = candidate,
                              frequency = c(0, 1))
print(encode_features(selected$encoder, te$x[1:3, , drop = FALSE]))
