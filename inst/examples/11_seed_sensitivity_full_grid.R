library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test
fit_ridge <- function(x, y) {
  z <- cbind(intercept = 1, x)
  solve(crossprod(z) + diag(c(0, rep(1, ncol(x)))), crossprod(z, y))
}
predict_ridge <- function(model, x) as.numeric(cbind(1, x) %*% model)
neg_mse <- function(y, prediction, domain) -mean((y - prediction)^2)

records <- list(); chosen <- list()
for (seed in c(11, 22, 33)) {
  set.seed(seed)
  folds <- sample(rep(1:3, length.out = nrow(tr$x)))
  for (method in c("center", "project")) {
    tuned <- tune_cross_domain(tr$x, tr$domain, tr$regression, folds,
      fit_ridge, predict_ridge, neg_mse, correction = method)
    records[[length(records) + 1]] <- data.frame(seed = seed, method = method, tuned$scores)
    chosen[[length(chosen) + 1]] <- data.frame(seed = seed, method = method,
      lambda = tuned$best_lambda, interior = tuned$best_lambda > 0 & tuned$best_lambda < 1)
  }
}
full_grid <- do.call(rbind, records); selected <- do.call(rbind, chosen)
print(selected)
stopifnot(nrow(full_grid) == 36, nrow(selected) == 6)
# Retain endpoint and interior results alike; seed searching is exploratory.
utils::write.csv(full_grid, file.path(tempdir(), "CrossDomainAdjust_all_seeds.csv"), row.names = FALSE)
