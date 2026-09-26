# Tests and reproducibility

Supported test environment: macOS and R 4.5.2 with the 70 package versions in `renv.lock`. Restore dependencies as described in [README](README.md). [VALIDATION_SOURCE.json](VALIDATION_SOURCE.json) identifies the tested source and dependency lock by SHA-256.

Run from the repository root:

```sh
Rscript --vanilla tests/test_engine.R
Rscript --vanilla tests/test_server.R
Rscript --vanilla tests/test_metric_support.R
Rscript --vanilla tests/test_plotting.R
```

The four suites contain 139 assertions:

- Engine: 44 assertions covering mappings, interval boundaries, four window methods, support holes, probability derivation, export provenance and PNG rendering.
- Server: 34 assertions covering CSV/RDS upload, Update and invalidation, failed mappings, group selection, missing probabilities and recovery from failed GAM fits.
- Metric support: 49 assertions covering separate score/probability populations, class minima, independent missingness, moderator geometry, export denominators and GAM weights.
- Plotting: 12 assertions covering the minimum supported-cell requirement, unusable weights, fit failure and recovery, and PNG output.

Tests generate inputs during execution and write temporary products outside the repository. No participant data or saved result snapshots are required.

For an interface check, start `shiny::runApp()` and exercise both dimensions with each window method. Check that changed settings invalidate results until Update, empty group selections disable downloads and explain the missing selection, and selecting groups again restores computation. Check CSV denominators against the displayed metric support. Unsupported or failed GAM fits must retain Raw/CSV access while disabling GAM PNG export; reducing the basis size on sufficiently supported data must allow recovery.

Window geometry uses finite moderators and known outcomes; each metric applies its own data requirements and class minima. In particular, 2D threshold metrics can use valid probabilities when scores are missing. These software tests do not establish clinical calibration, independent subgroup inference or upstream model validity.
