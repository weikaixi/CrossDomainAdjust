library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

if (!requireNamespace("survival", quietly = TRUE)) stop("Install the suggested package survival.")
library(survival)
fit_cox <- function(x, y) {
  z <- x; time <- y$time; event <- y$event
  m <- coxph(Surv(time, event) ~ ridge(z, theta = 10, scale = FALSE), ties = "efron")
  beta <- coef(m)
  stopifnot(all(is.finite(beta)))
  beta
}
predict_cox <- function(model, x) as.numeric(x %*% model)
cindex <- function(y, prediction, domain) {
  time <- y$time; event <- y$event; risk <- prediction
  as.numeric(concordance(Surv(time, event) ~ risk, reverse = TRUE)$concordance)
}

set.seed(9)
folds <- sample(rep(1:3, length.out = nrow(tr$x)))
tuned <- tune_cross_domain(tr$x, tr$domain, tr$survival, folds,
  fit_cox, predict_cox, cindex, correction = "center")
model <- fit_cox(transform_cross_domain(tuned$fit, tr$x, tr$domain), tr$survival)
pred <- predict_cox(model, transform_cross_domain(tuned$fit, te$x, te$domain))
print(tuned$scores)
cat("Pooled test C-index:", cindex(te$survival, pred, te$domain), "
")
# All patient replicates must share a fold in a real study.
