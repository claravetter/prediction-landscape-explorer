# Prediction Landscape Explorer

A standalone Shiny app for descriptive exploration of how classifier performance varies across continuous moderators. It computes local AUC, balanced accuracy, sensitivity and specificity from generated examples or an explicitly uploaded CSV/RDS data frame. The companion [PRONIA analysis workflow](https://github.com/claravetter/pronia-mri-transition-performance-landscapes) implements the study's inference and sensitivity analyses.

## Installation and launch

Use R 4.5.2. From this repository's root, restore the pinned R dependencies and launch:

```sh
Rscript --vanilla scripts/restore_environment.R
export R_LIBS_USER="$PWD/.Rlibrary"
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", launch.browser=TRUE)'
```

Restoration requires network access for missing packages and may need a compiler/system libraries. `renv.lock` records all 70 packages in the supported app/test runtime, including mgcv, nlme, Matrix and lattice. The restore script verifies their exact versions even when R supplies recommended packages. `.Rprofile` provides standard renv activation for interactive renv use; the commands above use the explicitly restored library. Direct dependencies include shiny, shinyjs, shinycssloaders, bslib, thematic, dplyr, tibble, ggplot2, readr, purrr, pROC, patchwork, mgcv, digest, jsonlite and renv. GAM smoothing is available in both 1D and 2D.

## Data and operation

The default example is generated in memory with a fixed seed (`R/demo_data.R`): 1,200 simulated records in five fictional groups. No data file, participant record or saved research result is bundled. Age and cognitive score affect the generating signal strength, with sampling noise; local AUC is not necessarily monotonic. Other moderators have no imposed signal in the generator. Its fitted demonstration probabilities do not represent clinical risk estimates.

For uploads, choose a CSV or RDS data frame, then verify the suggested mapping. Select the two outcome values explicitly (negative=0, positive=1), a numeric score, optional probability/group columns, and numeric moderators. Higher scores must indicate the positive outcome. Supplied probabilities must lie in [0,1]; invalid values are rejected. Non-finite scores (`Inf`, `-Inf`, `NaN`) are treated as missing. Without a supplied probability the app uses a logistic transform of standardized finite scores; a rank fallback uses only finite scores if standardization is unavailable. Missing/non-finite scores never produce substitute probabilities. Separately supplied valid probabilities remain usable independently of missing scores. This transformation is not a calibration procedure and cannot establish absolute clinical risk. Threshold metrics depend on that choice; AUC uses the score directly.

Internal moderator names `EXP_LABEL`, `Mean_Score`, `decision_z`, `prob_work` and `cohort` are reserved; if a moderator has one of these names, rename it to a non-reserved name before uploading. A probability source column with such a name is mapped safely without replacing another internal field. Uploaded column names must be unique. An upload should contain one row per intended unit; the app cannot establish unique participant identity without an identity mapping.

Probability is used only when selected in the column mapping. A column with a probability-like name selected solely as a moderator does not override the documented score transform.

After selecting data, group filters, moderators and settings, click **Update** on the relevant tab. Changes to data, mappings or settings invalidate previous plots/tables/downloads until another Update. A failed upload or invalid mapping clears prior data. Computation and presentation use the same captured settings. CSV exports include source indication, software version, computation time and full applied settings; PNG captions identify the source and principal settings. Exports contain aggregates; uploaded rows and local file paths are not included in export metadata. Group metadata records the effective filter: without a group mapping, `group` and `selected_groups` are null and `group_filter_applied` is false. With a group mapping, only selected groups present in that upload are applied and exported; clearing the selection excludes all rows.

## Methods and interpretation

The default in both dimensions is a fixed-range window: width 20% and step 5% of the selected moderator range, untrimmed, minimum seven events and seven non-events. Intervals are left-closed/right-open except the last interval, which includes the maximum. In 2D the intervals are crossed into boxes. These settings correspond to [`config/study.yml`](https://github.com/claravetter/pronia-mri-transition-performance-landscapes/blob/main/config/study.yml) in the companion analysis; the app is descriptive and does not perform the study's permutation tests, correction families or bootstrap inference.

Alternatives are fixed participant count in 1D (default N=100, 50% overlap) and K-nearest-neighbour windows on standardized coordinates in 2D. AUC is threshold-free; other metrics use the selected probability threshold. Changing range trimming, group selection or moderator changes the window population.

Window geometry uses rows with known outcomes and finite selected moderators. AUC additionally requires a finite score; BACC, sensitivity and specificity require a finite probability. Each metric must independently meet both class minima after these exclusions. Unsupported metric values are `NA`, even if the overall window is large. CSV exports report `auc_n_obs`, `auc_n_events`, `auc_n_nonevents`, `auc_supported` and the corresponding `threshold_*` columns. The on-screen table shows the count/support columns for the selected metric. The 1D methods omit windows whose overall event/non-event counts already miss the minima; 2D retains these candidate cells with missing metrics. A retained window can still have unsupported metrics after metric-specific exclusions. Overall `n_obs` and event counts describe the window before metric exclusions. Plot point sizes and 2D GAM weights use the selected metric's eligible counts.

Overlapping windows reuse observations. Marginal 2D bars and summed grid counts describe **memberships**, not unique people/events. Summaries separately report input rows/events eligible for the selected metric before windowing, and overall versus metric-specific memberships; these are not the union represented by supported cells. Local trends alone do not establish moderation, and flat trends do not establish its absence. Smoother bands do not account for dependence between overlapping estimates or replace inferential uncertainty. Unsupported cells and holes remain excluded from GAM displays. Raw grids preserve uneven centre spacing, and GAM support is matched by numeric grid positions rather than formatted coordinate strings. A 2D GAM display requires at least 20 usable metric cells and a successful model fit. Status, plot and PNG download share the same prepared plot; failed fits disable PNG export and provide guidance. CSV export remains available when a table exists. Choose Raw grid to display supported cells when a smooth surface cannot be fitted.

## Tests and files

```sh
Rscript --vanilla tests/test_engine.R
Rscript --vanilla tests/test_server.R
Rscript --vanilla tests/test_metric_support.R
Rscript --vanilla tests/test_plotting.R
```

`app.R` defines the interface and server. `R/engine.R` contains window/metric/plot functions; `R/results.R` binds computations to applied settings and exports; `R/demo_data.R` generates the demonstration. Tests cover numerical boundaries, upload mappings/collisions, all four window methods, finite-score probability transforms, stale-result invalidation, effective group metadata, and sparse/failed GAM display handling. See [validation](VALIDATION.md) for test coverage and execution instructions. Each test script prints PASS lines and a final assertion count.

## Licence and citation

MIT; see [LICENSE](LICENSE). This permits academic and commercial software reuse with the required notices and grants no rights to any separately uploaded data. Dependencies retain their own licences.

Cite Clara Vetter, *Prediction Landscape Explorer*, version 0.3.2, with the exact commit used and https://github.com/claravetter/prediction-landscape-explorer. A final study citation/DOI is not specified here; use the authors' confirmed paper/preprint citation when available. This app is a research exploration tool and does not establish clinical utility.
