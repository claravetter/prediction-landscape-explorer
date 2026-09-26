options(stringsAsFactors = FALSE)

# ---------------------------
# User-editable defaults
# ---------------------------
APP_VERSION <- "0.3.2"
DEFAULT_THRESHOLD <- 0.5
MIN_GAM_CELLS <- 20L

# Pretty labels for known column names. Anything not listed here gets
# displayed by its raw column name. Includes both the generic names
# used by the generated synthetic demo and study-specific names from the
# accompanying analysis pipeline (so real-data uploads with those
# columns also get readable labels).
MODERATOR_LABELS <- c(
  # Generic names (generated synthetic demo)
  "age"                   = "Age",
  "cognitive_score"       = "Cognitive score",
  "positive_symptoms"     = "Positive symptoms",
  "negative_symptoms"     = "Negative symptoms",
  "psychosocial_function" = "Psychosocial function",
  "brain_age_delta"       = "Brain-age delta",
  # Study-specific names (accompanying analysis pipeline)
  "age_z"                 = "Age (z)",
  "BrainAGE_corr"         = "BrainAGE (global)",
  "brainage_z"            = "BrainAGE (z)",
  "BrainAGE_SBC_corr"     = "BrainAGE_SBC (global)",
  "brainage_SBC_z"        = "BrainAGE_SBC (z)",
  "Psychosoz_aequiv"      = "Psychosocial functioning",
  "Psychosoz_aequiv_z"    = "Psychosocial functioning (z)",
  "psychosoz_z"           = "Psychosocial functioning (z)",
  "COGDIS_score"          = "COGDIS",
  "COGDIS_score_z"        = "COGDIS (z)",
  "cogdis_z"              = "COGDIS (z)",
  "SIPS_Positiv_Gesamt"   = "SIPS Positive",
  "SIPS_Positiv_Gesamt_z" = "SIPS Positive (z)",
  "sips_p_z"              = "SIPS Positive (z)",
  "SIPS_Negativ_Gesamt"   = "SIPS Negative",
  "SIPS_Negativ_Gesamt_z" = "SIPS Negative (z)",
  "sips_n_z"              = "SIPS Negative (z)"
)

# Alternative fixed-count 1D defaults (50% overlap)
DEFAULT_WINDOW_SIZE <- 100
DEFAULT_STEP_FRAC   <- 0.50
DEFAULT_MIN_EV      <- 7
DEFAULT_MIN_NONEV   <- 7
DEFAULT_HIST_MODE   <- "subject_bins"   # subject_bins | window_events
DEFAULT_HIST_BINS   <- 15

# Fixed-range window width/step, expressed as a fraction of each
# moderator's range for both 1D and 2D fixed-range methods.
# Companion configuration: config/study.yml in
# https://github.com/claravetter/pronia-mri-transition-performance-landscapes
DEFAULT_WIDTH_PCT <- 0.20
DEFAULT_STEP_PCT  <- 0.05

# 2D defaults
DEFAULT_GRID_RES     <- 15
DEFAULT_GRID_MODE    <- "quantile"  # quantile | range
# Use the complete finite moderator range by default.
DEFAULT_QTRIM_LO     <- 0
DEFAULT_QTRIM_HI     <- 1

DEFAULT_CELL_MODE    <- "fraction"  # fraction | fixed
DEFAULT_CELL_FRAC    <- 0.075
DEFAULT_CELL_N       <- 100

DEFAULT_SURFACE_MODE <- "gam"       # raw | gam
DEFAULT_GAM_FAMILY   <- "betar"     # betar | gaussian
DEFAULT_GAM_K        <- 30
DEFAULT_PRED_STEP    <- 0.25

# ---------------------------
# Data prep helpers
# ---------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

# Explicit two-class mapping: requires the caller to name exactly which raw
# value is the negative (0) and positive (1) class, and errors on anything
# else (multiclass, fractional codes, or values that don't match either
# name) rather than silently coercing with as.integer().
binary_mapping <- function(x, negative = "0", positive = "1") {
  x <- as.character(x)
  negative <- as.character(negative); positive <- as.character(positive)
  if (length(negative) != 1L || length(positive) != 1L ||
      is.na(negative) || is.na(positive) || negative == positive ||
      !nzchar(negative) || !nzchar(positive)) {
    stop("Choose two distinct outcome values to map to 0/1.")
  }
  observed <- unique(x[!is.na(x)])
  if (!setequal(observed, c(negative, positive))) {
    stop("Outcome column must contain exactly the two mapped values ('",
         negative, "' and '", positive, "'); found: ",
         paste(observed, collapse = ", "))
  }
  ifelse(is.na(x), NA_integer_, as.integer(x == positive))
}

make_prob_work <- function(df, prob_col = NULL, score_col = "Mean_Score") {
  # Probabilities read directly from a data column (rather than derived
  # in-app from decision_z/rank) are validated, not clamped: silently
  # clipping out-of-range values can mask a mis-mapped column.
  validate_prob <- function(p, source) {
    bad <- !is.na(p) & (!is.finite(p) | p < 0 | p > 1)
    if (any(bad)) {
      stop(source, " contains ", sum(bad),
           " value(s) outside [0, 1]; probabilities are not clipped — ",
           "fix the source data or choose a different column.")
    }
    p
  }
  p <- NULL
  if (!is.null(prob_col) && nzchar(prob_col)) {
    if (!prob_col %in% names(df) || !is.numeric(df[[prob_col]]))
      stop("Probability column must be numeric and present in the data.")
    p <- validate_prob(df[[prob_col]], paste0("Probability column '", prob_col, "'"))
  } else {
    scores <- if (!is.null(score_col) && score_col %in% names(df)) df[[score_col]] else rep(NA_real_, nrow(df))
    finite <- is.finite(scores)
    p <- rep(NA_real_, nrow(df))
    z_ok <- if ("decision_z" %in% names(df)) finite & is.finite(df$decision_z) else rep(FALSE, nrow(df))
    if (any(z_ok)) {
      p[z_ok] <- plogis(df$decision_z[z_ok])
    } else if (any(finite)) {
      ranks <- rank(scores[finite], ties.method = "average")
      p[finite] <- (ranks - 0.5) / sum(finite)
    }
  }
  df$prob_work <- p
  df
}

# Prepare an uploaded / generated dataset for the rest of the app.
#
# col_map: list with elements
#   outcome    – column name of a two-class outcome (required); mapped to
#                0/1 via `negative`/`positive`
#   negative   – raw value in `outcome` mapped to 0 (default "0")
#   positive   – raw value in `outcome` mapped to 1 (default "1")
#   score      – column name of the continuous classifier score (required)
#   prob       – column name of an optional pre-computed probability
#                in [0, 1] (e.g. Platt-calibrated). NULL/empty falls
#                back to a logistic standardized-score transform, or to
#                score ranks when no finite standardized scores exist.
#   cohort     – column name of an optional cohort/group factor.
#                NULL/empty disables cohort filtering.
#   moderators – character vector of column names to expose as
#                moderators (numeric). May be empty.
#
# Returns a tibble with the canonical names the compute layer expects:
#   EXP_LABEL, Mean_Score, decision_z, prob_work, cohort (optional).
# The user's chosen moderator columns are passed through under their
# original names, plus a z-scored sibling (`<col>_z`) is added when
# missing. Returns an attribute `moderator_cols` with the moderator
# names that survived prep (filtered to those numeric in df).
prep_data_for_app <- function(df, col_map) {
  stopifnot(is.data.frame(df), is.list(col_map))
  outcome <- col_map$outcome
  score   <- col_map$score
  prob    <- col_map$prob
  cohort  <- col_map$cohort
  mods    <- col_map$moderators %||% character()

  if (is.null(outcome) || !nzchar(outcome) || !outcome %in% names(df)) {
    stop("Outcome column not found in data.")
  }
  if (is.null(score) || !nzchar(score) || !score %in% names(df)) {
    stop("Score column not found in data.")
  }

  out <- tibble(
    EXP_LABEL = binary_mapping(df[[outcome]],
                               negative = col_map$negative %||% "0",
                               positive = col_map$positive %||% "1"),
    Mean_Score = suppressWarnings(as.numeric(df[[score]]))
  )
  out$Mean_Score[!is.finite(out$Mean_Score)] <- NA_real_
  if (!is.null(cohort) && nzchar(cohort) && cohort %in% names(df)) {
    out$cohort <- factor(df[[cohort]])
  }
  out$decision_z <- as.numeric(scale(out$Mean_Score))

  if (anyDuplicated(names(df))) stop("Column names must be unique.")
  mods <- unique(mods[mods %in% names(df)])
  mods <- mods[vapply(df[mods], is.numeric, logical(1))]
  reserved <- c("EXP_LABEL", "Mean_Score", "decision_z", "prob_work", "cohort")
  if (length(intersect(mods, reserved))) {
    stop("Rename moderator columns reserved for internal fields: ",
         paste(intersect(mods, reserved), collapse = ", "))
  }
  # Copy raw moderators first. A supplied *_z moderator is never silently
  # overwritten by a generated sibling, regardless of mapping order.
  for (m in mods) out[[m]] <- as.numeric(df[[m]])
  for (m in mods) {
    z_name <- paste0(m, "_z")
    if (!z_name %in% names(out)) out[[z_name]] <- as.numeric(scale(out[[m]]))
  }
  # Probability names stay in the upload namespace. They cannot overwrite
  # the mapped outcome, score, moderators, grouping or derived fields.
  if (!is.null(prob) && nzchar(prob)) {
    if (!prob %in% names(df) || !is.numeric(df[[prob]]))
      stop("Probability column must be numeric and present in the data.")
    out$prob_work <- df[[prob]]
    out <- make_prob_work(out, prob_col = "prob_work")
    attr(out, "probability_method") <- "supplied_probability"
  } else {
    out <- make_prob_work(out, score_col = "Mean_Score")
    attr(out, "probability_method") <- if (any(!is.na(out$decision_z)))
      "logistic_standardized_score (uncalibrated)" else "rank_score_transform (uncalibrated)"
  }

  attr(out, "moderator_cols") <- mods
  out
}

# Heuristics for auto-detecting candidate columns from an uploaded file
# (used by the column-mapping UI).
candidate_outcome_cols <- function(df) {
  ok <- vapply(df, function(x) {
    v <- as.character(x)
    nz <- v[!is.na(v)]
    length(nz) > 0 && length(unique(nz)) == 2
  }, logical(1))
  names(df)[ok]
}

# Default negative/positive class values for a two-valued column. Prefers
# the natural "0"/"1" ordering when present, otherwise sorts alphabetically.
default_negative_positive <- function(df, outcome) {
  if (is.null(outcome) || !nzchar(outcome) || !outcome %in% names(df)) {
    return(list(negative = "", positive = ""))
  }
  vals <- sort(unique(as.character(df[[outcome]])[!is.na(df[[outcome]])]))
  if (length(vals) != 2) return(list(negative = "", positive = ""))
  if (setequal(vals, c("0", "1"))) {
    list(negative = "0", positive = "1")
  } else {
    list(negative = vals[1], positive = vals[2])
  }
}

candidate_score_cols <- function(df) {
  names(df)[vapply(df, is.numeric, logical(1))]
}

candidate_prob_cols <- function(df) {
  num <- names(df)[vapply(df, is.numeric, logical(1))]
  num[vapply(num, function(n) {
    x <- df[[n]]; x <- x[!is.na(x)]
    length(x) > 0 && min(x) >= 0 && max(x) <= 1
  }, logical(1))]
}

candidate_cohort_cols <- function(df) {
  names(df)[vapply(df, function(x) {
    is.character(x) || is.factor(x) ||
      (is.numeric(x) && length(unique(x[!is.na(x)])) <= 20)
  }, logical(1))]
}

candidate_moderator_cols <- function(df, exclude = character()) {
  num <- names(df)[vapply(df, is.numeric, logical(1))]
  # Common pipeline byproducts that shouldn't be offered as moderators.
  pipeline_cols <- c("decision_z", "prob_work")
  setdiff(num, c(exclude, pipeline_cols))
}

# Preferred column names checked in order when auto-detecting a default
# mapping. The generated synthetic demo uses the first name in each list;
# study-specific names follow so real-data uploads with those columns
# still get a sensible default.
PREFERRED_OUTCOME_COLS <- c("event", "EXP_LABEL", "outcome", "y")
PREFERRED_SCORE_COLS   <- c("risk_score", "Mean_Score", "score")
PREFERRED_PROB_COLS    <- c("calibrated_prob", "Probs_platt", "prob")
PREFERRED_COHORT_COLS  <- c("group", "cohort", "site", "study")

# Default column mapping inferred directly from a raw data frame
# (used by auto-apply on example data, and as the seed for the
# column-mapping UI on upload). Picks columns sequentially, excluding
# anything already assigned so the same column isn't reused.
default_col_map_for <- function(df) {
  pick_preferred <- function(cands, preferred) {
    hit <- preferred[preferred %in% cands]
    if (length(hit) > 0) hit[1] else cands[1]
  }

  out_cands <- candidate_outcome_cols(df)
  outcome <- pick_preferred(out_cands, PREFERRED_OUTCOME_COLS)
  if (is.na(outcome)) outcome <- ""
  np <- default_negative_positive(df, outcome)

  score_cands <- setdiff(candidate_score_cols(df), outcome)
  score <- pick_preferred(score_cands, PREFERRED_SCORE_COLS)
  if (is.na(score)) score <- ""

  prob_cands <- setdiff(candidate_prob_cols(df), c(outcome, score))
  prob_hit <- PREFERRED_PROB_COLS[PREFERRED_PROB_COLS %in% prob_cands]
  prob <- if (length(prob_hit) > 0) prob_hit[1] else ""

  coh_cands <- setdiff(candidate_cohort_cols(df), c(outcome, score, prob))
  coh_hit <- PREFERRED_COHORT_COLS[PREFERRED_COHORT_COLS %in% coh_cands]
  cohort <- if (length(coh_hit) > 0) coh_hit[1] else ""

  excluded <- c(outcome, score,
                if (nzchar(prob)) prob else NULL,
                if (nzchar(cohort)) cohort else NULL)
  mods <- setdiff(candidate_moderator_cols(df, exclude = excluded), excluded)
  if (length(mods) > 8) mods <- mods[1:8]

  list(outcome = outcome, score = score, prob = prob,
       cohort = cohort, moderators = mods,
       negative = np$negative, positive = np$positive)
}

# ---------------------------
# Metrics
# ---------------------------
auc_from_scores <- function(y, score) {
  ok <- !is.na(y) & is.finite(score)
  y <- y[ok]; score <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  roc_obj <- pROC::roc(y, score, quiet = TRUE, direction = "<")
  as.numeric(roc_obj$auc)
}

conf_metrics <- function(y, prob, thr = 0.5) {
  ok <- !is.na(y) & is.finite(prob)
  y <- y[ok]; prob <- prob[ok]
  if (length(y) == 0 || length(unique(y)) < 2) {
    return(list(bacc = NA_real_, sens = NA_real_, spec = NA_real_))
  }
  pred <- as.integer(prob >= thr)
  tp <- sum(pred == 1 & y == 1)
  tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)

  sens <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  spec <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  bacc <- if (is.na(sens) || is.na(spec)) NA_real_ else 0.5 * (sens + spec)
  list(bacc = bacc, sens = sens, spec = spec)
}

# Window geometry uses finite moderators and a known binary outcome. Each
# metric then applies its own finite-value mask and support thresholds.
window_input <- function(df, moderators) {
  ok <- !is.na(df$EXP_LABEL) & df$EXP_LABEL %in% c(0, 1)
  for (name in moderators) ok <- ok & is.finite(df[[name]])
  df[ok, , drop = FALSE]
}
metric_input_mask <- function(df, metric) {
  ok <- !is.na(df$EXP_LABEL) & df$EXP_LABEL %in% c(0, 1)
  if (metric == "auc") ok & is.finite(df$Mean_Score)
  else ok & is.finite(df$prob_work) & df$prob_work >= 0 & df$prob_work <= 1
}
metric_count_column <- function(data, metric, events = FALSE) {
  prefix <- if (metric == "auc") "auc" else "threshold"
  name <- paste0(prefix, if (events) "_n_events" else "_n_obs")
  if (name %in% names(data)) return(name)
  if (!events) "n_obs" else if ("n_transition" %in% names(data)) "n_transition" else "n_events"
}
window_metrics <- function(sl, min_ev, min_nonev, thr) {
  a <- metric_input_mask(sl, "auc"); p <- metric_input_mask(sl, "bacc")
  ae <- sum(sl$EXP_LABEL[a] == 1); an <- sum(sl$EXP_LABEL[a] == 0)
  pe <- sum(sl$EXP_LABEL[p] == 1); pn <- sum(sl$EXP_LABEL[p] == 0)
  a_supported <- ae >= max(1, min_ev) && an >= max(1, min_nonev)
  p_supported <- pe >= max(1, min_ev) && pn >= max(1, min_nonev)
  cm <- if (p_supported) conf_metrics(sl$EXP_LABEL[p], sl$prob_work[p], thr)
        else list(bacc = NA_real_, sens = NA_real_, spec = NA_real_)
  tibble(auc = if (a_supported) auc_from_scores(sl$EXP_LABEL[a], sl$Mean_Score[a]) else NA_real_,
         bacc = cm$bacc, sens = cm$sens, spec = cm$spec,
         auc_n_obs = sum(a), auc_n_events = ae, auc_n_nonevents = an, auc_supported = a_supported,
         threshold_n_obs = sum(p), threshold_n_events = pe, threshold_n_nonevents = pn,
         threshold_supported = p_supported)
}

# ================================================================
# 1D sliding windows (step = fraction of window size)
# ================================================================
make_windows <- function(n_total, window_size = 100, step_size = 50) {
  if (n_total < window_size) return(list())
  # unique(c(..., last_start)) guarantees the final window always reaches
  # n_total, even when step_size doesn't evenly divide the usable range
  # (otherwise the last few sorted subjects are silently dropped).
  last_start <- n_total - window_size + 1
  starts <- unique(c(seq(1, last_start, by = step_size), last_start))
  lapply(starts, function(s) s:(s + window_size - 1))
}

calc_metrics_trajectory_1d <- function(df, sort_col,
                                       window_size = 100, step_frac = 0.5,
                                       min_ev = 7, min_nonev = 7,
                                       thr = 0.5) {
  if (!sort_col %in% names(df)) return(tibble())
  d <- window_input(df, sort_col) %>% arrange(.data[[sort_col]])
  n_total <- nrow(d)
  if (n_total < window_size) return(tibble())

  step_size <- max(1L, as.integer(round(step_frac * window_size)))
  step_size <- min(step_size, window_size - 1L)

  wins <- make_windows(n_total, window_size, step_size)
  if (length(wins) == 0) return(tibble())

  purrr::imap_dfr(wins, function(idx, i) {
    sl <- d[idx, , drop = FALSE]
    n_ev <- sum(sl$EXP_LABEL == 1, na.rm = TRUE)
    n_nonev <- sum(sl$EXP_LABEL == 0, na.rm = TRUE)
    if (n_ev < min_ev || n_nonev < min_nonev) return(NULL)

    metrics <- window_metrics(sl, min_ev, min_nonev, thr)

    bind_cols(tibble(
      window_id = i,
      center = mean(sl[[sort_col]], na.rm = TRUE),
      span = max(sl[[sort_col]], na.rm = TRUE) - min(sl[[sort_col]], na.rm = TRUE),
      n_obs = nrow(sl),
      n_events = n_ev,
      n_nonevents = n_nonev
    ), metrics)
  }) %>%
    mutate(method = "fixedn",
           step_size = step_size,
           window_size = window_size,
           overlap = 1 - step_size / window_size)
}

# ----------------------------------------------------------------
# Fixed-range 1D sliding windows
#
# Width and step are fractions of the moderator range, using the same
# per-axis windows as calc_grid_2d_fixedrange. See window_grid() and
# in_window() in the companion repository's R/lib/windows.R.
# ----------------------------------------------------------------
calc_metrics_fixedrange_1d <- function(df, sort_col,
                                       width_pct = DEFAULT_WIDTH_PCT,
                                       step_pct = DEFAULT_STEP_PCT,
                                       min_ev = 7, min_nonev = 7,
                                       thr = 0.5) {
  if (!sort_col %in% names(df)) return(tibble())
  d <- window_input(df, sort_col)
  x_raw <- d[[sort_col]]
  if (length(x_raw) == 0) return(tibble())

  x_min <- min(x_raw, na.rm = TRUE)
  x_max <- max(x_raw, na.rm = TRUE)
  rg <- x_max - x_min
  width <- rg * width_pct
  step  <- rg * step_pct

  win_grid <- .fixedrange_windows_1d(x_raw, width = width, step = step,
                                     x_min = x_min, x_max = x_max)
  if (nrow(win_grid) == 0) return(tibble())

  purrr::map_dfr(seq_len(nrow(win_grid)), function(i) {
    xl <- win_grid$x_lo[i]; xh <- win_grid$x_hi[i]
    idx <- if (win_grid$right_closed[i]) {
      which(x_raw >= xl & x_raw <= xh)
    } else {
      which(x_raw >= xl & x_raw < xh)
    }
    sl <- d[idx, , drop = FALSE]
    n_ev <- sum(sl$EXP_LABEL == 1, na.rm = TRUE)
    n_nonev <- sum(sl$EXP_LABEL == 0, na.rm = TRUE)
    if (n_ev < min_ev || n_nonev < min_nonev) return(NULL)

    metrics <- window_metrics(sl, min_ev, min_nonev, thr)

    bind_cols(tibble(
      window_id = i,
      center = win_grid$x_center[i],
      span = xh - xl,
      n_obs = nrow(sl),
      n_events = n_ev,
      n_nonevents = n_nonev
    ), metrics)
  }) %>%
    mutate(method = "fixedrange",
           width_pct = width_pct,
           step_pct = step_pct)
}

# Map a histogram quantity choice to the column to plot and its axis label.
# `mean_col` is required (and used) only when quantity == "mean_col".
hist_quantity_spec <- function(quantity = c("transitions", "subjects", "rate", "mean_col"),
                               mean_col = NULL) {
  quantity <- match.arg(quantity)
  switch(quantity,
    transitions = list(col = "transitions", lab = "N Transitions"),
    subjects    = list(col = "n",           lab = "N Subjects"),
    rate        = list(col = "rate",        lab = "Transition rate"),
    mean_col    = list(col = "mean_y",
                       lab = if (is.null(mean_col) || !nzchar(mean_col))
                               "Mean (no column selected)"
                             else paste("Mean", mean_col)))
}

# 1D histograms
hist_top_subject_bins <- function(df, sort_col, bins = 15,
                                  quantity = "transitions",
                                  mean_col = NULL) {
  d <- df %>% filter(!is.na(.data[[sort_col]]))
  if (nrow(d) == 0) return(ggplot() + theme_void())

  use_mean <- identical(quantity, "mean_col") &&
              !is.null(mean_col) && nzchar(mean_col) &&
              mean_col %in% names(d) && is.numeric(d[[mean_col]])

  brks <- pretty(range(d[[sort_col]], na.rm = TRUE), n = bins)
  d <- d %>%
    mutate(bin = cut(.data[[sort_col]], breaks = brks, include.lowest = TRUE))

  if (use_mean) {
    d <- d %>%
      group_by(bin) %>%
      summarise(
        center      = mean(.data[[sort_col]], na.rm = TRUE),
        transitions = sum(EXP_LABEL == 1, na.rm = TRUE),
        n           = dplyr::n(),
        mean_y      = mean(.data[[mean_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(rate = ifelse(n > 0, transitions / n, NA_real_))
  } else {
    d <- d %>%
      group_by(bin) %>%
      summarise(
        center      = mean(.data[[sort_col]], na.rm = TRUE),
        transitions = sum(EXP_LABEL == 1, na.rm = TRUE),
        n           = dplyr::n(),
        .groups = "drop"
      ) %>%
      mutate(rate = ifelse(n > 0, transitions / n, NA_real_))
  }

  spec <- hist_quantity_spec(quantity, mean_col = mean_col)
  ggplot(d, aes(x = center, y = .data[[spec$col]])) +
    geom_col(fill = "darkorange", alpha = 0.85) +
    labs(x = NULL, y = spec$lab) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_text(size = 9),
          plot.margin = margin(5, 5, 0, 5))
}

hist_top_window_events <- function(traj, quantity = "transitions") {
  if (nrow(traj) == 0) return(ggplot() + theme_void())
  if (identical(quantity, "mean_col")) quantity <- "transitions"  # not supported here
  d <- traj %>%
    transmute(
      center,
      transitions = n_events,
      n           = n_obs,
      rate        = ifelse(n_obs > 0, n_events / n_obs, NA_real_)
    )
  spec <- hist_quantity_spec(quantity)
  ggplot(d, aes(x = center, y = .data[[spec$col]])) +
    geom_col(fill = "darkorange", alpha = 0.85, width = 0.95) +
    labs(x = NULL, y = spec$lab) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_text(size = 9),
          plot.margin = margin(5, 5, 0, 5))
}

plot_1d_with_top_hist <- function(df, traj, sort_label, metric,
                                  hist_mode = "subject_bins", hist_bins = 15,
                                  hist_quantity = "transitions",
                                  hist_mean_col = NULL,
                                  zoom_auc = TRUE,
                                  base_size = 14,
                                  smooth_method = "loess") {

  if (nrow(traj) == 0) return(ggplot() + theme_void())

  # metric label
  ylab <- toupper(metric)
  if (metric == "auc")  ylab <- "AUC"
  if (metric == "bacc") ylab <- "Balanced accuracy"
  if (metric == "sens") ylab <- "Sensitivity"
  if (metric == "spec") ylab <- "Specificity"

  # y limits
  ylim_use <- c(0, 1)
  if (metric == "auc" && isTRUE(zoom_auc)) ylim_use <- c(0.4, 0.8)

  sort_col <- attr(traj, "sort_col")
  p_top <- if (hist_mode == "window_events") {
    hist_top_window_events(traj, quantity = hist_quantity)
  } else {
    hist_top_subject_bins(df, sort_col = sort_col, bins = hist_bins,
                          quantity = hist_quantity, mean_col = hist_mean_col)
  }

  # smoothing layer
  smooth_layer <- NULL
  if (smooth_method == "none") {
    smooth_layer <- NULL
  } else if (smooth_method == "lm") {
    smooth_layer <- geom_smooth(method = "lm", se = TRUE, alpha = 0.12)
  } else if (smooth_method == "loess") {
    smooth_layer <- geom_smooth(method = "loess", se = TRUE, span = 0.75, alpha = 0.12)
  } else if (smooth_method == "gam") {
    smooth_layer <- geom_smooth(method = "gam", formula = y ~ s(x, k = 10), se = TRUE, alpha = 0.12)
  } else if (smooth_method == "poly2") {
    smooth_layer <- geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, alpha = 0.12)
  }

  caption_txt <- if (identical(traj$method[1], "fixedrange")) {
    paste0(
      "fixed-range windows | width=", sprintf("%.0f%%", 100 * traj$width_pct[1]),
      " | step=", sprintf("%.0f%%", 100 * traj$step_pct[1]),
      " | n_windows=", nrow(traj)
    )
  } else {
    paste0(
      "window_size=", traj$window_size[1],
      " | step_frac=", sprintf("%.2f", traj$step_size[1] / traj$window_size[1]),
      " | step_size=", traj$step_size[1],
      " | overlap=", sprintf("%.0f%%", 100 * traj$overlap[1]),
      " | n_windows=", nrow(traj)
    )
  }

  count_col <- metric_count_column(traj, metric)
  p_main <- ggplot(traj, aes(x = center, y = .data[[metric]])) +
    geom_point(aes(size = .data[[count_col]]), alpha = 0.6) +
    smooth_layer +
    scale_size_area(max_size = 4, name = "Eligible N") +
    labs(
      title = paste0(sort_label, " — Sliding-window ", ylab),
      subtitle = caption_txt,
      x = sort_label,
      y = ylab
    ) +
    coord_cartesian(ylim = ylim_use) +
    theme_minimal(base_size = base_size) +
    theme(legend.position = "bottom",
          panel.grid = element_blank())

  if (metric %in% c("auc","bacc")) {
    p_main <- p_main + geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60")
  }

  p_top / p_main + plot_layout(heights = c(1, 3))
}

# ================================================================
# Moderators registry (for UI choices)
#
# Builds a tibble (label, sort_col) over whatever moderator columns
# the user selected in the column-mapping UI (or, on first load,
# whatever was attached as attr(df, "moderator_cols") by
# prep_data_for_app). The MODERATOR_LABELS lookup gives nice labels
# for recognised study columns; unknown columns display as their raw
# name.
# ================================================================
build_moderator_registry <- function(df, moderator_cols = NULL,
                                     label_map = MODERATOR_LABELS) {
  if (is.null(moderator_cols)) {
    moderator_cols <- attr(df, "moderator_cols") %||% character()
  }
  moderator_cols <- moderator_cols[moderator_cols %in% names(df)]
  if (length(moderator_cols) == 0) return(tibble(label = character(), sort_col = character()))

  labels <- vapply(moderator_cols, function(col) {
    if (col %in% names(label_map)) unname(label_map[[col]]) else col
  }, character(1))

  tibble(label = unname(labels), sort_col = moderator_cols)
}

# ================================================================
# 2D grid (KNN local windows) + plotting
# ================================================================
safe_scale_params <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sx <- sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) sx <- 1
  list(mu = mu, sd = sx, z = (x - mu) / sx)
}

calc_grid_2d_knn <- function(df, x_col, y_col,
                             grid_res = 15,
                             grid_mode = c("quantile","range"),
                             qtrim_lo = DEFAULT_QTRIM_LO,
                             qtrim_hi = DEFAULT_QTRIM_HI,
                             cell_n_mode = c("fraction","fixed"),
                             cell_n = 100,
                             cell_frac = 0.075,
                             min_ev = 7,
                             min_nonev = 7,
                             thr = 0.5,
                             show_progress = TRUE) {

  grid_mode <- match.arg(grid_mode)
  cell_n_mode <- match.arg(cell_n_mode)

  stopifnot(x_col %in% names(df), y_col %in% names(df))

  d <- window_input(df, c(x_col, y_col))

  if (nrow(d) < 20) stop("Not enough finite-moderator rows with a known outcome for 2D grid.")

  # choose k
  if (cell_n_mode == "fraction") {
    k <- ceiling(cell_frac * nrow(d))
  } else {
    k <- as.integer(cell_n)
  }
  k <- max(5L, k)
  k <- min(k, nrow(d))

  sx <- safe_scale_params(d[[x_col]])
  sy <- safe_scale_params(d[[y_col]])

  d <- d %>% mutate(x_norm = sx$z, y_norm = sy$z)

  x_raw <- d[[x_col]]
  y_raw <- d[[y_col]]

  if (grid_mode == "quantile") {
    probs <- seq(qtrim_lo, qtrim_hi, length.out = grid_res)
    x_seq_raw <- as.numeric(quantile(x_raw, probs = probs, na.rm = TRUE, names = FALSE))
    y_seq_raw <- as.numeric(quantile(y_raw, probs = probs, na.rm = TRUE, names = FALSE))
  } else {
    x_seq_raw <- seq(min(x_raw, na.rm = TRUE), max(x_raw, na.rm = TRUE), length.out = grid_res)
    y_seq_raw <- seq(min(y_raw, na.rm = TRUE), max(y_raw, na.rm = TRUE), length.out = grid_res)
  }

  x_seq_norm <- (x_seq_raw - sx$mu) / sx$sd
  y_seq_norm <- (y_seq_raw - sy$mu) / sy$sd

  grid_idx <- expand.grid(ix = seq_along(x_seq_norm),
                          iy = seq_along(y_seq_norm),
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)

  Xn <- d$x_norm
  Yn <- d$y_norm

  out_list <- vector("list", nrow(grid_idx))

  do_one <- function(ii) {
    ix <- grid_idx$ix[ii]
    iy <- grid_idx$iy[ii]
    gx <- x_seq_norm[ix]
    gy <- y_seq_norm[iy]

    dist2 <- (Xn - gx)^2 + (Yn - gy)^2
    idx <- order(dist2)[seq_len(k)]

    sl <- d[idx, , drop = FALSE]
    n_ev <- sum(sl$EXP_LABEL == 1, na.rm = TRUE)
    n_nonev <- sum(sl$EXP_LABEL == 0, na.rm = TRUE)

    bind_cols(tibble(
      x_center = x_seq_raw[ix],
      y_center = y_seq_raw[iy],
      n_obs = nrow(sl),
      n_transition = n_ev,
      n_nonevents = n_nonev
    ), window_metrics(sl, min_ev, min_nonev, thr))
  }

  if (show_progress) {
    withProgress(message = "Computing 2D grid (KNN windows)...", value = 0, {
      for (ii in seq_len(nrow(grid_idx))) {
        out_list[[ii]] <- do_one(ii)
        incProgress(1 / nrow(grid_idx))
      }
    })
  } else {
    for (ii in seq_len(nrow(grid_idx))) out_list[[ii]] <- do_one(ii)
  }

  g <- bind_rows(out_list) %>%
    mutate(k = k,
           grid_res = grid_res,
           grid_mode = grid_mode,
           qtrim_lo = qtrim_lo,
           qtrim_hi = qtrim_hi,
           cell_n_mode = cell_n_mode,
           cell_frac = ifelse(cell_n_mode == "fraction", cell_frac, NA_real_),
           cell_n = ifelse(cell_n_mode == "fixed", k, NA_real_))
  g
}

# ----------------------------------------------------------------
# Fixed-range 2D sliding windows
#
# Cross each axis's overlapping intervals into rectangular windows.
# Widths and steps are fractions of the quantile-bounded moderator range:
# width = range * width_pct, step = range * step_pct.
# The companion repository uses window_grid() and in_window() in
# R/lib/windows.R for the same interval convention.
# ----------------------------------------------------------------
.fixedrange_windows_1d <- function(x, width, step,
                                   x_min = NULL, x_max = NULL) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0 || !is.finite(width) || !is.finite(step) ||
      width <= 0 || step <= 0) {
    return(tibble(x_lo = numeric(0), x_hi = numeric(0), x_center = numeric(0),
                  right_closed = logical(0)))
  }
  if (is.null(x_min)) x_min <- min(x, na.rm = TRUE)
  if (is.null(x_max)) x_max <- max(x, na.rm = TRUE)
  if (!is.finite(x_min) || !is.finite(x_max)) {
    return(tibble(x_lo = numeric(0), x_hi = numeric(0), x_center = numeric(0),
                  right_closed = logical(0)))
  }
  if (x_max <= x_min) {
    ctr <- x_min
    return(tibble(x_lo = ctr - width/2, x_hi = ctr + width/2, x_center = ctr,
                  right_closed = TRUE))
  }
  if (width >= (x_max - x_min)) {
    ctr <- (x_min + x_max) / 2
    return(tibble(x_lo = ctr - width/2, x_hi = ctr + width/2, x_center = ctr,
                  right_closed = TRUE))
  }
  target <- x_max - width
  starts <- seq(from = x_min, to = target, by = step)
  if (length(starts) == 0) {
    ctr <- (x_min + x_max) / 2
    return(tibble(x_lo = ctr - width/2, x_hi = ctr + width/2, x_center = ctr,
                  right_closed = TRUE))
  }
  # Force the last window to reach x_max exactly: step doesn't generally
  # divide (target - x_min) evenly, so without this the final subjects
  # (those above the last window's x_hi) fall in an uncovered gap and are
  # silently excluded from every window.
  tolerance <- max(1, abs(x_max), abs(x_min)) * 1e-9
  if (abs(starts[length(starts)] - target) > tolerance) {
    starts <- c(starts, target)
  } else {
    starts[length(starts)] <- target
  }
  tibble(x_lo = starts, x_hi = starts + width, x_center = starts + width/2,
         right_closed = seq_along(starts) == length(starts))
}

calc_grid_2d_fixedrange <- function(df, x_col, y_col,
                                    width_pct_x = DEFAULT_WIDTH_PCT, step_pct_x = DEFAULT_STEP_PCT,
                                    width_pct_y = DEFAULT_WIDTH_PCT, step_pct_y = DEFAULT_STEP_PCT,
                                    qtrim_lo = DEFAULT_QTRIM_LO, qtrim_hi = DEFAULT_QTRIM_HI,
                                    min_ev = 7, min_nonev = 7,
                                    thr = 0.5,
                                    show_progress = TRUE) {
  stopifnot(x_col %in% names(df), y_col %in% names(df))

  d <- window_input(df, c(x_col, y_col))
  if (nrow(d) < 20) stop("Not enough finite-moderator rows with a known outcome for 2D grid.")

  x_raw <- d[[x_col]]
  y_raw <- d[[y_col]]

  # qtrim'd ranges (consistent with the KNN method's quantile-grid bounds)
  x_min <- as.numeric(quantile(x_raw, probs = qtrim_lo, na.rm = TRUE, names = FALSE))
  x_max <- as.numeric(quantile(x_raw, probs = qtrim_hi, na.rm = TRUE, names = FALSE))
  y_min <- as.numeric(quantile(y_raw, probs = qtrim_lo, na.rm = TRUE, names = FALSE))
  y_max <- as.numeric(quantile(y_raw, probs = qtrim_hi, na.rm = TRUE, names = FALSE))

  rg_x <- x_max - x_min
  rg_y <- y_max - y_min

  width_x <- rg_x * width_pct_x; step_x <- rg_x * step_pct_x
  width_y <- rg_y * width_pct_y; step_y <- rg_y * step_pct_y

  grid_x <- .fixedrange_windows_1d(x_raw, width = width_x, step = step_x,
                                   x_min = x_min, x_max = x_max)
  grid_y <- .fixedrange_windows_1d(y_raw, width = width_y, step = step_y,
                                   x_min = y_min, x_max = y_max)

  if (nrow(grid_x) == 0 || nrow(grid_y) == 0) {
    return(tibble(x_center = numeric(0), y_center = numeric(0)))
  }

  cells <- expand.grid(ix = seq_len(nrow(grid_x)),
                       iy = seq_len(nrow(grid_y)),
                       KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)

  do_one <- function(ii) {
    ix <- cells$ix[ii]; iy <- cells$iy[ii]
    xl <- grid_x$x_lo[ix]; xh <- grid_x$x_hi[ix]
    yl <- grid_y$x_lo[iy]; yh <- grid_y$x_hi[iy]
    # Half-open intervals [lo, hi), except
    # the final window along each axis, which is right-closed [lo, hi] so
    # subjects exactly at the domain maximum aren't dropped from every window.
    x_hit <- if (grid_x$right_closed[ix]) x_raw >= xl & x_raw <= xh else x_raw >= xl & x_raw < xh
    y_hit <- if (grid_y$right_closed[iy]) y_raw >= yl & y_raw <= yh else y_raw >= yl & y_raw < yh
    idx <- which(x_hit & y_hit)
    sl <- d[idx, , drop = FALSE]
    n_ev <- sum(sl$EXP_LABEL == 1, na.rm = TRUE)
    n_nonev <- sum(sl$EXP_LABEL == 0, na.rm = TRUE)

    row_base <- tibble(
      x_center = grid_x$x_center[ix], y_center = grid_y$x_center[iy],
      x_lo = xl, x_hi = xh, y_lo = yl, y_hi = yh,
      n_obs = nrow(sl), n_transition = n_ev, n_nonevents = n_nonev
    )
    bind_cols(row_base, window_metrics(sl, min_ev, min_nonev, thr))
  }

  out_list <- vector("list", nrow(cells))
  if (show_progress) {
    withProgress(message = "Computing 2D grid (fixed-range windows)...",
                 value = 0, {
      for (ii in seq_len(nrow(cells))) {
        out_list[[ii]] <- do_one(ii)
        incProgress(1 / nrow(cells))
      }
    })
  } else {
    for (ii in seq_len(nrow(cells))) out_list[[ii]] <- do_one(ii)
  }

  g <- bind_rows(out_list) %>%
    mutate(k = NA_integer_,
           grid_res = nrow(grid_x) * nrow(grid_y),
           grid_mode = "fixedrange",
           qtrim_lo = qtrim_lo, qtrim_hi = qtrim_hi,
           cell_n_mode = "fixedrange",
           cell_frac = NA_real_, cell_n = NA_real_,
           width_pct_x = width_pct_x, step_pct_x = step_pct_x,
           width_pct_y = width_pct_y, step_pct_y = step_pct_y)
  g
}

hist_top_from_grid <- function(grid, x = "x_center", quantity = "transitions") {
  if (identical(quantity, "mean_col")) quantity <- "transitions"  # not supported here
  dat <- grid %>%
    group_by(.data[[x]]) %>%
    summarise(
      transitions = sum(n_transition, na.rm = TRUE),
      n           = sum(n_obs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(rate = ifelse(n > 0, transitions / n, NA_real_))

  spec <- hist_quantity_spec(quantity)
  spec$lab <- switch(quantity, transitions = "Event memberships", subjects = "Row memberships", rate = "Event fraction of memberships")
  ggplot(dat, aes(x = .data[[x]], y = .data[[spec$col]])) +
    geom_col(fill = "darkorange", alpha = 0.85) +
    labs(x = NULL, y = spec$lab) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_text(size = 9),
          plot.margin = margin(5, 5, 0, 5))
}

hist_right_from_grid <- function(grid, y = "y_center", quantity = "transitions") {
  if (identical(quantity, "mean_col")) quantity <- "transitions"  # not supported here
  dat <- grid %>%
    group_by(.data[[y]]) %>%
    summarise(
      transitions = sum(n_transition, na.rm = TRUE),
      n           = sum(n_obs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(rate = ifelse(n > 0, transitions / n, NA_real_))

  spec <- hist_quantity_spec(quantity)
  spec$lab <- switch(quantity, transitions = "Event memberships", subjects = "Row memberships", rate = "Event fraction of memberships")
  ggplot(dat, aes(x = .data[[y]], y = .data[[spec$col]])) +
    geom_col(fill = "darkorange", alpha = 0.85) +
    coord_flip() +
    labs(x = spec$lab, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 9),
          plot.margin = margin(5, 0, 5, 5))
}

plot_grid_raw <- function(grid, metric = "auc",
                          xlab = "X", ylab = "Y",
                          overlay_points = TRUE) {

  fill_lab <- if (metric == "auc") "AUC" else if (metric == "bacc") "BACC" else toupper(metric)

  # Quantile centres (and the final fixed-range step) can be uneven. Use
  # explicit midpoint boundaries so raster regularisation cannot shift values.
  bounds <- function(x) {
    u <- sort(unique(x))
    edges <- if (length(u) == 1L) c(u - .5, u + .5) else
      c(u[1] - diff(u)[1]/2, (head(u,-1) + tail(u,-1))/2,
        tail(u,1) + tail(diff(u),1)/2)
    i <- match(x, u)
    list(lo = edges[i], hi = edges[i+1L])
  }
  bx <- bounds(grid$x_center); by <- bounds(grid$y_center)
  grid$plot_xlo <- bx$lo; grid$plot_xhi <- bx$hi
  grid$plot_ylo <- by$lo; grid$plot_yhi <- by$hi
  p <- ggplot(grid, aes(fill = .data[[metric]])) +
    geom_rect(aes(xmin = plot_xlo, xmax = plot_xhi, ymin = plot_ylo, ymax = plot_yhi)) +
    scale_fill_viridis_c(name = fill_lab, na.value = "white") +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())

  if (overlay_points) {
    g2 <- grid %>% filter(!is.na(.data[[metric]]))
    event_col <- metric_count_column(g2, metric, events = TRUE)
    p <- p + geom_point(data = g2,
                        aes(x = x_center, y = y_center, size = .data[[event_col]]),
                        inherit.aes = FALSE,
                        colour = "white", alpha = 0.55) +
      scale_size_continuous(name = "Eligible events per cell", range = c(0.2, 3)) +
      guides(size = guide_legend(override.aes = list(colour = "grey25", alpha = 1)))
  }
  p
}

gam_input <- function(grid, metric, weights_col = NULL) {
  if (is.null(weights_col)) weights_col <- metric_count_column(grid, metric)
  stopifnot(all(c(metric, "x_center", "y_center", weights_col) %in% names(grid)))
  grid %>% filter(is.finite(.data[[metric]]), is.finite(x_center), is.finite(y_center),
                  is.finite(.data[[weights_col]]), .data[[weights_col]] > 0)
}

gam_surface_from_grid <- function(grid,
                                  metric = "auc",
                                  step = 0.25,
                                  family = c("betar", "gaussian"),
                                  k = 30,
                                  weights_col = NULL,
                                  overlay_points = TRUE) {

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is not installed. Install it or choose 'Raw grid' mode.")
  }

  fam <- match.arg(family)
  if (is.null(weights_col)) weights_col <- metric_count_column(grid, metric)

  df <- gam_input(grid, metric, weights_col)
  if (nrow(df) < MIN_GAM_CELLS) stop("Too few non-missing grid cells to fit a smooth surface.")

  if (fam == "betar") {
    resp <- df[[metric]]
    eps  <- 1e-6
    resp <- pmin(1 - eps, pmax(eps, resp))
    n    <- length(resp)
    df$.__y__ <- (resp * (n - 1) + 0.5) / n
    famobj <- mgcv::betar(link = "logit")
  } else {
    df$.__y__ <- df[[metric]]
    famobj <- gaussian()
  }

  form <- as.formula(paste0(".__y__ ~ s(x_center, y_center, k = ", k, ")"))

  fit <- mgcv::gam(form,
                   family = famobj,
                   data = df,
                   weights = df[[weights_col]],
                   method = "REML")

  xr <- range(df$x_center, na.rm = TRUE)
  yr <- range(df$y_center, na.rm = TRUE)

  x_seq <- seq(xr[1], xr[2], by = step)
  y_seq <- seq(yr[1], yr[2], by = step)

  newdat <- expand.grid(x_seq, y_seq, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  names(newdat) <- c("x_center", "y_center")
  newdat$pred <- as.numeric(predict(fit, newdata = newdat, type = "response"))
  # Mask predictions outside the metric's actual spatial support so the
  # smoothed surface can't visually bridge a support hole or extrapolate
  # past the observed grid.
  newdat$pred[!supported_grid_points(grid, newdat, metric)] <- NA_real_

  fill_lab <- if (metric == "auc") "AUC (pred)" else if (metric == "bacc") "BACC (pred)" else paste0(toupper(metric), " (pred)")

  p <- ggplot(newdat, aes(x = x_center, y = y_center, fill = pred)) +
    geom_raster() +
    scale_fill_viridis_c(name = fill_lab, na.value = "grey85") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())

  if (overlay_points) {
    event_col <- metric_count_column(df, metric, events = TRUE)
    p <- p +
      geom_point(data = df,
                 aes(x = x_center, y = y_center, size = .data[[event_col]]),
                 inherit.aes = FALSE,
                 colour = "white", alpha = 0.55) +
      scale_size_continuous(name = "Eligible events per cell", range = c(0.2, 3)) +
      guides(size = guide_legend(override.aes = list(colour = "grey25", alpha = 1)))
  }

  list(model = fit, surface = newdat, plot = p)
}

# Interpolation is only trusted within rectangles whose four surrounding
# original grid corners are non-missing for `metric`; otherwise a smooth
# surface would fill unsupported/sparse regions with a colour that visually
# implies real data coverage.
supported_grid_points <- function(grid, newdat, metric) {
  xs <- sort(unique(grid$x_center[is.finite(grid$x_center)]))
  ys <- sort(unique(grid$y_center[is.finite(grid$y_center)]))
  if (!length(xs) || !length(ys)) return(rep(FALSE, nrow(newdat)))
  # Integer positions avoid vector-dependent decimal formatting/padding.
  # Repeated quantile centres are supported only if all their cells are finite.
  key <- match(grid$x_center, xs) + (match(grid$y_center, ys)-1L) * length(xs)
  flags <- tapply(is.finite(grid[[metric]]), key, all)
  supported <- rep(FALSE, length(xs)*length(ys))
  supported[as.integer(names(flags))] <- flags
  vapply(seq_len(nrow(newdat)), function(i) {
    x <- newdat$x_center[i]; y <- newdat$y_center[i]
    if (!is.finite(x) || !is.finite(y) || x < min(xs) || x > max(xs) || y < min(ys) || y > max(ys)) return(FALSE)
    xl <- max(which(xs <= x)); xh <- min(which(xs >= x))
    yl <- max(which(ys <= y)); yh <- min(which(ys >= y))
    corners <- expand.grid(x = unique(c(xl, xh)), y = unique(c(yl, yh)))
    all(supported[corners$x + (corners$y-1L)*length(xs)])
  }, logical(1))
}

combine_2d_with_marginals <- function(main_plot, top_hist, right_hist) {
  (top_hist + plot_spacer()) /
    (main_plot + right_hist) +
    plot_layout(widths = c(3, 1), heights = c(1, 3), guides = "collect") &
    theme(legend.position = "right")
}
