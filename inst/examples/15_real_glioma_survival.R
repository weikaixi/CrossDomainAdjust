library(CrossDomainAdjust)

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

# Public-source, frozen 60-gene worked data. These examples do not reselect genes.
# Fixed illustrative strengths and penalty: NOT a reproduction of paper estimates.
d <- readRDS(system.file("extdata", "glioma_worked_data.rds", package = "CrossDomainAdjust"))
tr <- d$train; te <- d$test
print(d$provenance)
out <- list()
for (method in c("raw", "rank", "pair", "center", "project")) {
  fit <- if (method == "raw") fit_cross_domain(tr$x, tr$domain, correction = "none") else
    fit_domain_method(tr$x, tr$domain, method, lambda = .4, max_pairs = 128)
  model <- fit_cox(transform_cross_domain(fit, tr$x, tr$domain), tr$survival)
  pred <- predict_cox(model, transform_cross_domain(fit, te$x, te$domain))
  out[[method]] <- data.frame(method = method, cohort = d$test_cohort,
    n = nrow(te$x), pooled_C_index = cindex(te$survival, pred, te$domain))
}
print(do.call(rbind, out), row.names = FALSE)
# Test data are transformed with training parameters. No test outcome enters a fit.
# This previously inspected cohort is a worked example, not a new confirmatory test.
