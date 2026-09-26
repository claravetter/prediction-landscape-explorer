# A result binds computed values, data and applied settings. Only aggregate
# values/settings are exported; uploaded records and internal hashes stay local.
apply_group_filter <- function(df, selected = NULL) {
  if (!"cohort" %in% names(df))
    return(list(data = df, selected_groups = NULL, applied = FALSE))
  used <- intersect(as.character(selected), as.character(unique(df$cohort[!is.na(df$cohort)])))
  list(data = droplevels(df %>% filter(cohort %in% used)), selected_groups = used, applied = TRUE)
}

no_groups_selected <- function(request) {
  mapping <- request$mapping
  isTRUE(mapping$group_filter_applied) && !length(mapping$selected_groups)
}

result_defaults <- function(dim) {
  if (dim == "1d") list(moderator = "", method = "fixedrange", width_pct = .2, step_pct = .05,
    window_size = 100, step_frac = .5, min_ev = 7, min_nonev = 7, threshold = .5,
    metric = "auc", hist_mode = "subject_bins", hist_bins = 15, hist_quantity = "transitions",
    hist_mean_col = "Mean_Score", zoom_auc = TRUE, base_size = 14, smooth_method = "none")
  else list(x_var = "", y_var = "", method = "fixedrange", width_pct_x = .2, step_pct_x = .05,
    width_pct_y = .2, step_pct_y = .05, qtrim_lo = 0, qtrim_hi = 1, min_ev = 7, min_nonev = 7,
    threshold = .5, grid_res = 15, grid_mode = "quantile", cell_mode = "fraction", cell_n = 100,
    cell_frac = .075, metric = "auc", surface_mode = "gam", pred_step = .25,
    gam_family = "betar", gam_k = 30, overlay_points = TRUE, hist_quantity = "transitions")
}
compute_result <- function(request, dim, show_progress = TRUE) {
  df <- request$data; s <- request$settings
  if (no_groups_selected(request)) {
    tab <- tibble()
    moderators <- if (dim == "1d") s$moderator else c(s$x_var, s$y_var)
  } else if (dim == "1d") {
    tab <- if (s$method == "fixedrange")
      calc_metrics_fixedrange_1d(df, s$moderator, s$width_pct, s$step_pct, s$min_ev, s$min_nonev, s$threshold)
    else calc_metrics_trajectory_1d(df, s$moderator, s$window_size, s$step_frac, s$min_ev, s$min_nonev, s$threshold)
    attr(tab, "sort_col") <- s$moderator
    moderators <- s$moderator
  } else {
    if (s$x_var == s$y_var) stop("Choose two different variables for X and Y.")
    tab <- if (s$method == "fixedrange")
      calc_grid_2d_fixedrange(df, s$x_var, s$y_var, s$width_pct_x, s$step_pct_x,
        s$width_pct_y, s$step_pct_y, s$qtrim_lo, s$qtrim_hi, s$min_ev, s$min_nonev,
        s$threshold, show_progress = show_progress)
    else calc_grid_2d_knn(df, s$x_var, s$y_var, s$grid_res, s$grid_mode, s$qtrim_lo,
        s$qtrim_hi, s$cell_mode, s$cell_n, s$cell_frac, s$min_ev, s$min_nonev,
        s$threshold, show_progress = show_progress)
    moderators <- c(s$x_var, s$y_var)
  }
  windowing <- window_input(df, moderators)
  eligible <- windowing[metric_input_mask(windowing, s$metric), , drop = FALSE]
  list(request = request, table = tab, dimension = dim, version = APP_VERSION,
       created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
       windowing_rows = nrow(windowing), windowing_events = sum(windowing$EXP_LABEL == 1),
       eligible_rows = nrow(eligible), eligible_events = sum(eligible$EXP_LABEL == 1))
}
result_metadata <- function(r) {
  s <- r$request$settings
  # Mapping provenance contains variable names, not uploaded values or paths.
  list(source = r$request$source, software_version = r$version, computed_utc = r$created,
       dimension = r$dimension, settings = s, mapping = r$request$mapping,
       probability_method = attr(r$request$data, "probability_method"),
       windowing_rows_before_metric_exclusions = r$windowing_rows,
       windowing_events_before_metric_exclusions = r$windowing_events,
       eligible_rows_before_windowing = r$eligible_rows,
       eligible_events_before_windowing = r$eligible_events,
       support_rule = "AUC: finite score and outcome; threshold metrics: finite probability and outcome. Each metric needs the configured event and non-event minima after exclusions.",
       counts = "n_obs/event counts describe window memberships before metric exclusions; auc_* and threshold_* describe metric-specific eligible counts/support. Eligible input rows refer to the selected metric and count records once, not necessarily unique people.")
}
result_export <- function(r) {
  if (!nrow(r$table)) stop("No windows to export.")
  out <- r$table
  out$source <- r$request$source
  out$software_version <- r$version
  out$computed_utc <- r$created
  out$applied_settings_json <- as.character(jsonlite::toJSON(result_metadata(r), auto_unbox = TRUE, null = "null"))
  out
}
result_filename <- function(r, extension) {
  s <- r$request$settings; dim <- r$dimension
  fields <- if (dim == "1d") s$moderator else paste(s$x_var, s$y_var, sep = "_x_")
  tag <- if (s$method == "fixedrange") {
    if (dim == "1d") sprintf("fixedrange_w%g_s%g", s$width_pct, s$step_pct)
    else sprintf("fixedrange_w%g_%g_s%g_%g", s$width_pct_x, s$width_pct_y, s$step_pct_x, s$step_pct_y)
  } else if (dim == "1d") sprintf("fixedn_n%g_step%g", s$window_size, s$step_frac)
    else paste("knn", s$grid_mode, s$grid_res, s$cell_mode,
               if (s$cell_mode == "fraction") s$cell_frac else s$cell_n, sep = "_")
  stem <- paste(r$request$source, dim, fields, s$metric, tag, sep = "_")
  paste0(gsub("[^A-Za-z0-9_.-]", "_", stem), ".", extension)
}
result_summary <- function(r) {
  g <- r$table; s <- r$request$settings
  if (!nrow(g)) return(tibble())
  tibble(method = s$method, metric = s$metric,
    eligible_input_rows = r$eligible_rows, eligible_input_events = r$eligible_events,
    windows_or_cells = nrow(g), valid_metric_cells = sum(is.finite(g[[s$metric]])),
    mean_metric = mean(g[[s$metric]], na.rm = TRUE), sd_metric = sd(g[[s$metric]], na.rm = TRUE),
    row_memberships = sum(g$n_obs),
    event_memberships = sum(if (r$dimension == "1d") g$n_events else g$n_transition),
    metric_row_memberships = sum(g[[metric_count_column(g, s$metric)]]),
    metric_event_memberships = sum(g[[metric_count_column(g, s$metric, events = TRUE)]]))
}
plot_unavailable_reason <- function(r) {
  g <- r$table; s <- r$request$settings
  if (no_groups_selected(r$request))
    return("No groups selected. Select at least one group and click Update.")
  if (!nrow(g)) return("No windows to plot. Adjust settings and click Update.")
  if (!any(is.finite(g[[s$metric]])))
    return("No supported values for this metric. Check missing values or adjust settings and click Update.")
  if (r$dimension == "2d" && s$surface_mode == "gam" && nrow(gam_input(g, s$metric)) < MIN_GAM_CELLS)
    return(paste0("GAM needs at least ", MIN_GAM_CELLS, " supported cells. Choose Raw grid or adjust windows and click Update. CSV remains available."))
  NULL
}

prepare_result_plot <- function(r) {
  reason <- plot_unavailable_reason(r)
  if (!is.null(reason)) return(list(ok = FALSE, plot = NULL, message = reason))
  tryCatch(list(ok = TRUE, plot = plot_result(r), message = NULL), error = function(e) {
    message <- if (r$dimension == "2d" && r$request$settings$surface_mode == "gam")
      "Unable to fit the GAM surface. Choose Raw grid, reduce the GAM basis size, or adjust windows and click Update. CSV remains available."
    else "Unable to prepare this plot. Check the display settings and click Update. CSV remains available."
    list(ok = FALSE, plot = NULL, message = message)
  })
}

plot_result <- function(r) {
  g <- r$table; s <- r$request$settings; df <- r$request$data
  reason <- plot_unavailable_reason(r)
  if (!is.null(reason)) stop(reason)
  label <- function(n) if (n %in% names(MODERATOR_LABELS)) unname(MODERATOR_LABELS[n]) else n
  if (r$dimension == "1d") {
    p <- plot_1d_with_top_hist(window_input(df, s$moderator), g, sort_label = label(s$moderator), metric = s$metric,
      hist_mode = s$hist_mode, hist_bins = s$hist_bins, hist_quantity = s$hist_quantity,
      hist_mean_col = s$hist_mean_col, zoom_auc = s$zoom_auc, base_size = s$base_size,
      smooth_method = s$smooth_method)
  } else {
    if (s$surface_mode == "raw") {
      main <- plot_grid_raw(g, s$metric, label(s$x_var), label(s$y_var), s$overlay_points) +
        labs(title = paste("Raw grid", toupper(s$metric)))
    } else {
      main <- gam_surface_from_grid(g, metric = s$metric, step = s$pred_step,
        family = s$gam_family, k = s$gam_k, weights_col = metric_count_column(g, s$metric),
        overlay_points = s$overlay_points)$plot +
        labs(x = label(s$x_var), y = label(s$y_var), title = paste("GAM-smoothed", toupper(s$metric)))
    }
    p <- combine_2d_with_marginals(main, hist_top_from_grid(g, quantity = s$hist_quantity),
                                  hist_right_from_grid(g, quantity = s$hist_quantity))
  }
  window_keys <- if (s$method == "fixedrange") {
    if (r$dimension == "1d") c("width_pct", "step_pct")
    else c("width_pct_x", "step_pct_x", "width_pct_y", "step_pct_y")
  } else if (r$dimension == "1d") c("window_size", "step_frac")
    else c("grid_res", "grid_mode", "cell_mode", if (s$cell_mode == "fixed") "cell_n" else "cell_frac")
  display_keys <- if (r$dimension == "1d") "smooth_method" else
    c("qtrim_lo", "qtrim_hi", "surface_mode", if (s$surface_mode == "gam") c("gam_family", "gam_k", "pred_step"))
  active <- s[c("method", window_keys, "threshold", "min_ev", "min_nonev", display_keys)]
  caption <- paste(r$request$source, "| Explorer", r$version, "|", r$created,
                   "\n", paste(paste(names(active), unlist(active), sep = "="), collapse = "; "),
                   "\nOverlapping windows: memberships are not independent participants. CSV contains full applied settings.")
  p + patchwork::plot_annotation(caption = paste(strwrap(caption, width = 125), collapse = "\n"))
}
