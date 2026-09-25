# CrossDomainAdjust 2.1.0 validation

Local release validation, 2026-09-24:

- Built and installed from source on Windows 11, R 4.6.1 (x86_64, UCRT).
- `R CMD check --no-manual --no-build-vignettes`: **Status: OK**, no errors, warnings or notes.
- All five test scripts passed: domain adjustment, feature encoding, scaling,
  fold-local pipelines, and the new preset/diagnostic/structured-outcome workflows.
- All 15 numbered worked examples completed, including the two public-source survival examples.
- Ten frozen preprocessing objects from the earlier 2.0.1 lung/glioma analyses
  (raw, rank, pair, center, projection in each) produced exactly identical matrices
  under 2.1.0 on the checked external samples: maximum absolute difference **0**.
- Rank/pair behavior, training-only scaling, pair screening within folds,
  fixed-axis lambda diagnostics, structured outcome alignment, single-sample
  behavior and caller RNG preservation were included in verification.

The manual PDF was not built during `R CMD check`; full English/Chinese Markdown guides
and R help pages are included. The check result is a software validation result, not
evidence of clinical validity or superiority of a method.

The GitHub workflow independently builds/checks the package and runs all 15 examples
on Ubuntu. Its current status is shown under the repository's Actions tab; this local
validation statement does not assume that a remote workflow has already completed.
