library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

fit <- fit_domain_method(tr$x, tr$domain, "project", lambda = .4)
diag <- lambda_diagnostics(fit, tr$x, tr$domain)
print(diag$summary)
xy <- diag$coordinates
file <- file.path(tempdir(), "CrossDomainAdjust_lambda_grid.pdf")
grDevices::pdf(file, width = 10, height = 6)
op <- graphics::par(mfrow = c(2, 3), mar = c(3.5, 3.5, 2, 1))
for (strength in diag$summary$lambda) {
  z <- xy[xy$lambda == strength, ]
  graphics::plot(z$PC1, z$PC2, col = as.integer(factor(z$domain)), pch = 16, cex = .5,
    xlim = range(xy$PC1), ylim = range(xy$PC2), xlab = "Fixed PC1", ylab = "Fixed PC2",
    main = paste("Lambda =", strength))
}
graphics::par(op); grDevices::dev.off()
cat("Fixed-axis PDF:", file, "
")
