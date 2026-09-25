library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test
fit_ridge <- function(x, y) {
  z <- cbind(intercept = 1, x)
  solve(crossprod(z) + diag(c(0, rep(1, ncol(x)))), crossprod(z, y))
}
predict_ridge <- function(model, x) as.numeric(cbind(1, x) %*% model)
neg_mse <- function(y, prediction, domain) -mean((y - prediction)^2)

set.seed(7)
folds <- sample(rep(1:3, length.out = nrow(tr$x)))
tuned <- tune_cross_domain(tr$x, tr$domain, tr$regression, folds,
  fit_ridge, predict_ridge, neg_mse, correction = "project")
print(tuned$scores)
z <- transform_cross_domain(tuned$fit, tr$x, tr$domain)
model <- fit_ridge(z, tr$regression)
pred <- predict_ridge(model, transform_cross_domain(tuned$fit, te$x))
cat("Chosen lambda:", tuned$best_lambda, "Untouched test RMSE:",
    sqrt(mean((te$regression - pred)^2)), "
")
