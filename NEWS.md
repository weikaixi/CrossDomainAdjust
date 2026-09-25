# CrossDomainAdjust 2.1.0

- Added four study-aligned presets through `fit_domain_method()`.
- Added fixed-axis PCA and full-feature lambda diagnostics without mutating a fitted object.
- Added a reproducible simulation helper that restores the caller's random-number state.
- Extended fold-local tuning to matrix, data-frame and `Surv` outcomes.
- Standardized the default lambda grid to 0, 0.2, 0.4, 0.6, 0.8, 1.
- Added 15 independently runnable examples, including two public-source survival workflows,
  and detailed English and Chinese guides.
- Preserved the transformation behavior of frozen 2.0.1 models; this release does not refit
  or replace the manuscript's frozen analyses.

# CrossDomainAdjust 2.0.1

- Added optional training-fitted feature standardization for rank and pair representations.
- Added frozen raw, rank, pair, center, projection and advanced hybrid pipelines.
- Added training-only pair frequency screening and fold-local lambda tuning.
