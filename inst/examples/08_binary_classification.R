library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

fit <- fit_domain_method(tr$x, tr$domain, "rank")
z <- transform_cross_domain(fit, tr$x)[, -ncol(tr$x), drop = FALSE]
# Fractional ranks have a fixed sum. Omit one column to avoid exact collinearity
# with the intercept in an unpenalized logistic model (ridge models need no such drop).
model <- stats::glm.fit(cbind(1, z), tr$binary, family = stats::binomial())
stopifnot(model$converged, all(is.finite(model$coefficients)))
test <- transform_cross_domain(fit, te$x)[, colnames(z), drop = FALSE]
prob <- stats::plogis(as.numeric(cbind(1, test) %*% model$coefficients))
cat("Test Brier score:", mean((te$binary - prob)^2), "
")
cat("Test accuracy at a fixed 0.5 threshold:", mean((prob >= .5) == te$binary), "
")
