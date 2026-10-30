library(CrossDomainAdjust)
d <- simulate_domain_data(p = 8, seed = 42)
tr <- d$train; te <- d$test

fit <- fit_domain_method(tr$x, tr$domain, "rank")
file <- tempfile(fileext = ".rds")
saveRDS(fit, file)
restored <- readRDS(file)
# Shuffled named columns are restored to the original training order.
original <- transform_cross_domain(fit, te$x)
shuffled <- transform_cross_domain(restored, te$x[, rev(colnames(te$x))])
stopifnot(isTRUE(all.equal(original, shuffled)))
missing <- tryCatch(transform_cross_domain(restored, te$x[, -1]),
                     error = function(e) conditionMessage(e))
cat("Expected missing-gene error:", missing, "
")
stopifnot(is.character(missing))
unlink(file)
