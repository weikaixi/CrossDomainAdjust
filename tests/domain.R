library(CrossDomainAdjust)

near <- function(a, b, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance)))
}
must_error <- function(expr, pattern = NULL) {
  got <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(got, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(got)))
}

x <- rbind(c(1, 2, 5), c(3, 4, 7), c(5, 2, 8), c(7, 4, 10),
           c(3, 7, 5), c(5, 9, 7))
dimnames(x) <- list(paste0("s", seq_len(nrow(x))), paste0("g", seq_len(ncol(x))))
domain <- rep(c("A", "B", "C"), each = 2L)
centers <- rbind(A = c(2, 3, 6), B = c(6, 3, 9), C = c(4, 8, 6))
colnames(centers) <- colnames(x)
anchor <- colMeans(centers)

center_fit <- fit_domain_adjuster(x, domain, "center")
near(center_fit$domain_centers, centers)
near(center_fit$anchor, anchor)
near(adjust_features(center_fit, x, domain, 0.4),
     x - 0.4 * sweep(centers[domain, , drop = FALSE], 2L, anchor))
centered <- adjust_features(center_fit, x, domain, 1)
for (d in unique(domain)) near(colMeans(centered[domain == d, , drop = FALSE]), anchor)
stopifnot(identical(dimnames(centered), dimnames(x)))
# Centers remain frozen for future observations, even when their batch mean changes.
new_x <- x[1:2, , drop = FALSE] + 100
near(adjust_features(center_fit, new_x, c("A", "A")),
     new_x - matrix(centers["A", ] - anchor, 2L, 3L, byrow = TRUE))
must_error(adjust_features(center_fit, new_x), "requires domain")
must_error(adjust_features(center_fit, new_x, c("new", "new")), "Unknown domain")
near(adjust_features(center_fit, new_x, c("new", "new"), lambda = 0), new_x)

projection <- fit_domain_adjuster(x, domain, "project")
B <- projection$basis
stopifnot(ncol(B) == 2L, identical(dim(B), c(3L, 2L)))
near(crossprod(B), diag(2L))
near(adjust_features(projection, x, lambda = 0.35),
     x - 0.35 * sweep(x, 2L, anchor) %*% B %*% t(B))
projected <- adjust_features(projection, x, lambda = 1)
near(adjust_features(projection, projected, lambda = 1), projected)
for (d in unique(domain)) near(colMeans(projected[domain == d, , drop = FALSE]), anchor)
near(adjust_features(projection, x, rep("unseen", nrow(x))), projected)
near(adjust_features(projection, x, lambda = 0), x)
stopifnot(identical(dimnames(projected), dimnames(x)))
# New input columns are reordered to the fitted feature schema.
near(adjust_features(projection, x[, 3:1, drop = FALSE]), projected)
near(adjust_features(projection, x[1L, ]), projected[1L, , drop = FALSE])

custom_anchor <- c(g3 = -2, g1 = 8, g2 = 4)
custom <- fit_domain_adjuster(x, domain, "project", anchor = custom_anchor)
custom_result <- adjust_features(custom, x, lambda = 0.2)
near(custom_result, x - 0.2 * sweep(x, 2L, custom_anchor[colnames(x)]) %*%
       custom$basis %*% t(custom$basis))
near(adjust_features(custom, custom_anchor), custom_anchor[colnames(x)])

hybrid <- fit_domain_adjuster(x, domain, "hybrid")
first <- adjust_features(center_fit, x, domain, lambda = 0.3)
near(adjust_features(hybrid, x, domain, c(0.3, 0.6)),
     adjust_features(projection, first, lambda = 0.6))
near(adjust_features(hybrid, x, domain, 0.4),
     adjust_features(hybrid, x, domain, c(0.4, 0.4)))
near(adjust_features(hybrid, x, rep("unseen", nrow(x)), c(0, 0.6)),
     adjust_features(projection, x, lambda = 0.6))
must_error(adjust_features(hybrid, x, rep("unseen", nrow(x)), c(0.2, 0)), "Unknown domain")

# Rank-zero and single-feature matrices are valid and retain matrix dimensions.
same <- rbind(x[1:2, ], x[1:2, ])
zero <- fit_domain_adjuster(same, rep(c("A", "B"), each = 2L))
stopifnot(zero$rank == 0L, ncol(zero$basis) == 0L)
near(adjust_features(zero, same), same)
one <- matrix(c(1, 3, 7, 9), ncol = 1L,
              dimnames = list(paste0("r", 1:4), "only"))
onefit <- fit_domain_adjuster(one, c("A", "A", "B", "B"))
near(adjust_features(onefit, one), rep(5, 4L))
stopifnot(identical(dim(adjust_features(onefit, one)), dim(one)))
onecenter <- fit_domain_adjuster(one, c("A", "A", "B", "B"), "center")
near(adjust_features(onecenter, one, c("A", "A", "B", "B")), c(4, 6, 4, 6))
none <- fit_domain_adjuster(one, rep("A", nrow(one)), "none")
near(adjust_features(none, one), one)

# Unequal domain sizes distinguish equal-domain and sample-weighted anchors.
unequal <- rbind(x, x[1L, , drop = FALSE])
unequal_domain <- c(domain, "A")
weighted <- fit_domain_adjuster(unequal, unequal_domain, anchor = "sample_mean")
near(weighted$anchor, colMeans(unequal))
robust <- fit_domain_adjuster(x, domain, center = "median")
near(robust$domain_centers, centers)
trimmed <- fit_domain_adjuster(x, domain, center = "trimmed", trim = 0.2)
near(trimmed$domain_centers, centers)

limited <- fit_domain_adjuster(x, domain, n_components = 1L)
stopifnot(limited$rank == 1L)
limited_zero <- fit_domain_adjuster(x, domain, n_components = 0L)
near(adjust_features(limited_zero, x), x)
must_error(fit_domain_adjuster(x, domain, n_components = 4L), "n_components")

# Low-rank storage remains linear in p times number of domains.
set.seed(12)
wide <- matrix(rnorm(24L * 1500L), 24L, 1500L,
               dimnames = list(NULL, paste0("g", seq_len(1500L))))
wide_fit <- fit_domain_adjuster(wide, rep(c("A", "B", "C"), each = 8L))
stopifnot(ncol(wide_fit$basis) <= 2L)
stopifnot(!any(vapply(wide_fit, function(z) is.matrix(z) &&
                       identical(dim(z), c(1500L, 1500L)), logical(1))))

bad_names <- x
colnames(bad_names)[2L] <- colnames(bad_names)[1L]
must_error(fit_domain_adjuster(bad_names, domain), "unique")
bad_values <- x
bad_values[1L, 1L] <- Inf
must_error(fit_domain_adjuster(bad_values, domain), "finite")
bad_values[1L, 1L] <- NA_real_
must_error(fit_domain_adjuster(bad_values, domain), "finite")
must_error(adjust_features(projection, x[, 1:2, drop = FALSE]), "Missing")
must_error(fit_domain_adjuster(x, rep("A", nrow(x))), "two training domains")
must_error(fit_domain_adjuster(x, c(domain[-1L], NA_character_)), "domain")
must_error(adjust_features(projection, x, lambda = -0.1), "lambda")
must_error(adjust_features(projection, x, lambda = c(0.2, 0.3)), "lambda")
must_error(adjust_features(projection, x, lambda = NaN), "lambda")

# Legacy entry points accept a named sample vector and optional unseen domain.
legacy <- fit_domain_projector(x, domain)
near(project_features(legacy, x, lambda = 0.5),
     adjust_features(projection, x, lambda = 0.5))
sample_result <- project_new_sample(legacy, x[1L, ], domain = "unseen", lambda = 0.5)
near(sample_result, project_features(legacy, x[1L, , drop = FALSE], lambda = 0.5))
stopifnot(identical(names(sample_result), colnames(x)),
          identical(attr(sample_result, "domain"), "unseen"),
          identical(attr(sample_result, "lambda"), 0.5))
must_error(project_new_sample(legacy, x[1:2, , drop = FALSE]), "exactly one row")
