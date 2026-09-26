#!/usr/bin/env Rscript
# Prediction Landscape Explorer. See README.md for setup, methods and limits.
# Synthetic examples are generated in memory; uploads are processed locally.

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinycssloaders)
  library(bslib)
  library(thematic)
  # Individual tidyverse packages we actually use (avoids pulling in
  # tidyr / stringr which are not needed here).
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(pROC)
  library(patchwork)
})

# Auto-sync ggplot styling (fonts, colours) to the bslib theme below.
thematic::thematic_shiny(font = "auto")

APP_THEME <- bslib::bs_theme(
  version = 5,
  bootswatch = "flatly",
  base_font = bslib::font_google("Inter"),
  heading_font = bslib::font_google("Inter")
)

source("R/engine.R", local = TRUE)
source("R/demo_data.R", local = TRUE)
source("R/results.R", local = TRUE)

# ================================================================
# UI
# ================================================================
ui <- fluidPage(
  title = "Prediction Landscape Explorer",
  useShinyjs(),
  theme = APP_THEME,
  tags$head(tags$style(HTML("
    .scroll-table {
      max-height: 520px;
      overflow: auto;
      border: 1px solid #e5e5e5;
      border-radius: 4px;
      padding: 0.25rem 0.5rem;
      background: #fff;
    }
    .scroll-table table { margin-bottom: 0; }
    .scroll-table thead th {
      position: sticky; top: 0; background: #fff; z-index: 2;
      border-bottom: 2px solid #dee2e6;
    }
    .tbl-count {
      color: #6c757d; font-size: 0.85rem;
      margin: 0 0 0.25rem 0.1rem;
    }
  "))),
  div(class = "app-header",
      style = "padding: 1rem 0 0.5rem 0; border-bottom: 1px solid #e5e5e5; margin-bottom: 1rem;",
      h2("Prediction Landscape Explorer",
         style = "margin: 0; font-weight: 600;"),
      p(style = "margin: 0.25rem 0 0 0; color: #6c757d; font-size: 0.95rem;",
        "Interactive companion to the multi-cohort transition-prediction validation study. ",
        "Explore how classifier performance (AUC, balanced accuracy, sensitivity, specificity) ",
        "varies with continuous moderators using fixed-range or fixed-count 1D windows and fixed-range or KNN 2D grids.")
  ),

  uiOutput("synthetic_banner"),

  sidebarLayout(
    sidebarPanel(
      h4("Data"),
      radioButtons("data_source", NULL,
                   choices = c("Example (synthetic)" = "example",
                               "Upload your own (.csv/.rds)" = "upload"),
                   selected = "example"),
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput("data_upload", "Upload (.csv/.rds)",
                  accept = c(".rds", ".csv"))
      ),
      uiOutput("column_mapping_ui"),
      helpText("Check the column mapping, then click Update. Changed settings invalidate previous results."),
      hr(),

      h4("Filter"),
      uiOutput("cohort_selector"),
      hr(),

      # 1D controls
      conditionalPanel(
        condition = "input.main_tab == '1D'",
        h4("1D sliding-window settings"),
        uiOutput("moderator_selector_1d"),
        selectInput("metric_1d", "Metric",
                    choices = c("AUC"="auc",
                                "Balanced accuracy"="bacc",
                                "Sensitivity"="sens",
                                "Specificity"="spec"),
                    selected = "auc"),
        actionButton("update_1d", "Update 1D plot",
                     class = "btn-primary btn-block",
                     style = "width:100%; margin: 8px 0;"),
        hr(),
        h5("Window settings"),
        numericInput("threshold_1d", "Threshold (for BACC/Sens/Spec)", value = DEFAULT_THRESHOLD, min = 0.01, max = 0.99, step = 0.01),
        radioButtons("method_1d", "1D windowing method",
                     choices = c("Fixed-range windows (default)" = "fixedrange",
                                 "Fixed participant count (exploratory)"    = "fixedn"),
                     selected = "fixedrange"),
        conditionalPanel(
          condition = "input.method_1d == 'fixedrange'",
          h5("Window width and step (fraction of moderator range)"),
          sliderInput("width_pct_1d", "Width fraction",
                      min = 0.05, max = 0.75, value = DEFAULT_WIDTH_PCT, step = 0.01),
          sliderInput("step_pct_1d", "Step fraction",
                      min = 0.01, max = 0.5,  value = DEFAULT_STEP_PCT, step = 0.01)
        ),
        conditionalPanel(
          condition = "input.method_1d == 'fixedn'",
          sliderInput("window_size_1d", "Window size (N per window)", min = 50, max = 500, value = DEFAULT_WINDOW_SIZE, step = 10),
          sliderInput("step_frac_1d", "Step fraction of window size", min = 0.05, max = 0.95, value = DEFAULT_STEP_FRAC, step = 0.05)
        ),
        hr(),
        h5("Advanced settings"),
        numericInput("min_ev_1d", "Min events per window", value = DEFAULT_MIN_EV, min = 0, step = 1),
        numericInput("min_nonev_1d", "Min non-events per window", value = DEFAULT_MIN_NONEV, min = 0, step = 1),
        checkboxInput("zoom_auc_1d", "Zoom AUC to 0.4–0.8", value = TRUE),
        sliderInput("base_size_1d", "Plot font size", min = 10, max = 28, value = 14, step = 1),
        selectInput("smooth_method_1d", "Smoother",
                    choices = c("LOESS"="loess",
                                "Linear (lm)"="lm",
                                "GAM"="gam",
                                "Polynomial (deg 2)"="poly2",
                                "None"="none"),
                    selected = "loess"),
        hr(),
        h5("Top histogram"),
        selectInput("hist_mode_1d", "Histogram mode",
                    choices = c("Subject bins (fewer bars)"="subject_bins",
                                "Per-window events (one bar per window)"="window_events"),
                    selected = DEFAULT_HIST_MODE),
        sliderInput("hist_bins_1d", "Bins (subject bins mode)", min = 5, max = 40, value = DEFAULT_HIST_BINS, step = 1),
        selectInput("hist_quantity_1d", "Bar plot quantity",
                    choices = c("N transitions" = "transitions",
                                "N subjects"    = "subjects",
                                "Transition rate" = "rate",
                                "Mean of a numeric column" = "mean_col"),
                    selected = "transitions"),
        conditionalPanel(
          condition = "input.hist_quantity_1d == 'mean_col'",
          uiOutput("hist_mean_col_selector_1d"),
          helpText("Only available in 'Subject bins' mode.")
        ),
        hr(),
        downloadButton("download_plot_1d", "Download 1D plot (PNG)"),
        downloadButton("download_csv_1d",  "Download 1D window table (CSV)")
      ),

      # 2D controls
      conditionalPanel(
        condition = "input.main_tab == '2D'",
        h4("2D surface settings"),
        uiOutput("xy_selector_2d"),
        selectInput("metric_2d", "Metric",
                    choices = c("AUC"="auc",
                                "Balanced accuracy"="bacc",
                                "Sensitivity"="sens",
                                "Specificity"="spec"),
                    selected = "bacc"),
        actionButton("update_2d", "Update 2D plot",
                     class = "btn-primary btn-block",
                     style = "width:100%; margin: 8px 0;"),
        hr(),
        h5("Window settings"),
        numericInput("threshold_2d", "Threshold (for BACC/Sens/Spec)", value = DEFAULT_THRESHOLD, min = 0.01, max = 0.99, step = 0.01),

        hr(),
        h5("2D windowing method"),
        radioButtons("method_2d", NULL,
                     choices = c("Fixed-range windows (default)" = "fixedrange",
                                 "K-nearest neighbours (exploratory)"      = "knn"),
                     selected = "fixedrange"),

        # ---- Fixed-range specific inputs ----
        conditionalPanel(
          condition = "input.method_2d == 'fixedrange'",
          h5("Window widths and steps (fractions of each moderator's qtrim'd range)"),
          sliderInput("width_pct_x_2d", "X width fraction",
                      min = 0.05, max = 0.75, value = DEFAULT_WIDTH_PCT, step = 0.01),
          sliderInput("step_pct_x_2d", "X step fraction",
                      min = 0.01, max = 0.5,  value = DEFAULT_STEP_PCT, step = 0.01),
          sliderInput("width_pct_y_2d", "Y width fraction",
                      min = 0.05, max = 0.75, value = DEFAULT_WIDTH_PCT, step = 0.01),
          sliderInput("step_pct_y_2d", "Y step fraction",
                      min = 0.01, max = 0.5,  value = DEFAULT_STEP_PCT, step = 0.01)
        ),

        # ---- KNN specific inputs ----
        conditionalPanel(
          condition = "input.method_2d == 'knn'",
          h5("KNN grid definition"),
          selectInput("grid_mode_2d", "Grid mode", choices = c("Quantile grid (trimmed)"="quantile", "Range grid (min–max)"="range"), selected = DEFAULT_GRID_MODE),
          sliderInput("grid_res_2d", "Grid resolution (per axis)", min = 5, max = 30, value = DEFAULT_GRID_RES, step = 1),
          h5("KNN cell size"),
          selectInput("cell_mode_2d", "Cell size mode", choices = c("Fraction of filtered N"="fraction", "Fixed N"="fixed"), selected = DEFAULT_CELL_MODE),
          conditionalPanel(
            condition = "input.cell_mode_2d == 'fraction'",
            sliderInput("cell_frac_2d", "Cell size fraction", min = 0.01, max = 0.25, value = DEFAULT_CELL_FRAC, step = 0.005)
          ),
          conditionalPanel(
            condition = "input.cell_mode_2d == 'fixed'",
            sliderInput("cell_n_2d", "Local window size k (subjects per cell)", min = 30, max = 400, value = DEFAULT_CELL_N, step = 5)
          )
        ),

        # ---- Shared range/cell inputs ----
        numericInput("qtrim_lo_2d", "Quantile trim (low)", value = DEFAULT_QTRIM_LO, min = 0, max = 0.49, step = 0.005),
        numericInput("qtrim_hi_2d", "Quantile trim (high)", value = DEFAULT_QTRIM_HI, min = 0.51, max = 1, step = 0.005),
        numericInput("min_ev_2d", "Min events per cell", value = DEFAULT_MIN_EV, min = 0, step = 1),
        numericInput("min_nonev_2d", "Min non-events per cell", value = DEFAULT_MIN_NONEV, min = 0, step = 1),
        hr(),
        h5("Display settings"),
        selectInput("surface_mode_2d", "Heatmap mode",
                    choices = c("GAM-smoothed surface"="gam",
                                "Raw grid (tiles)"="raw"),
                    selected = DEFAULT_SURFACE_MODE),
        conditionalPanel(
          condition = "input.surface_mode_2d == 'gam'",
          selectInput("gam_family_2d", "GAM family",
                      choices = c("Beta regression (betar)"="betar",
                                  "Gaussian"="gaussian"),
                      selected = DEFAULT_GAM_FAMILY),
          sliderInput("gam_k_2d", "GAM basis dimension k", min = 10, max = 60, value = DEFAULT_GAM_K, step = 1),
          sliderInput("pred_step_2d", "Prediction grid step", min = 0.1, max = 1.0, value = DEFAULT_PRED_STEP, step = 0.05)
        ),
        checkboxInput("overlay_points_2d", "Overlay grid points (size by eligible events)", value = TRUE),
        hr(),
        h4("Marginal bar plots"),
        selectInput("hist_quantity_2d", "Bar plot quantity",
                    choices = c("Event memberships" = "transitions",
                                "Row memberships"    = "subjects",
                                "Event fraction of memberships" = "rate"),
                    selected = "transitions"),
        hr(),
        downloadButton("download_plot_2d", "Download 2D plot (PNG)"),
        downloadButton("download_csv_2d",  "Download 2D grid table (CSV)")
      )
    ),

    mainPanel(
      tabsetPanel(id = "main_tab",
                  tabPanel("Instructions",
                           div(style = "max-width: 820px; padding: 0.5rem 0.25rem;",
                               h3("What this app does"),
                               p("This explorer estimates how a binary-outcome classifier's ",
                                 "performance changes across a continuous moderator (e.g. age or ",
                                 "symptom load). It is the interactive companion to the ",
                                 "multi-cohort transition-prediction validation study."),

                               h4("Two analysis modes"),
                               tags$ul(
                                 tags$li(tags$b("1D sliding window"), " — slide windows across a single sorted ",
                                         "moderator and compute AUC / BACC / sensitivity / specificity in each. ",
                                         "Two windowing methods are supported:"),
                                 tags$li(tags$b("2D moderation surface"), " — compute local performance on a ",
                                         "grid spanned by two moderators. Two windowing methods are supported:")
                               ),
                               tags$ul(
                                 tags$li(tags$b("Fixed-range windows (default)"),
                                         " — partition the moderator's range into overlapping intervals of fixed ",
                                         "width (expressed as a fraction of that range), stepped by a fixed ",
                                         "fraction. In 1D each interval is one window; in 2D, intervals along X ",
                                         "and Y are crossed into boxes. Default widths, steps and minimum ",
                                         "class counts follow the companion analysis configuration."),
                                 tags$li(tags$b("Fixed participant count / K-nearest neighbours (exploratory)"),
                                         " — 1D: fixed-size windows of N sorted subjects, stepped by a fraction ",
                                         "of the window size. 2D: a regular grid (quantile- or range-based) where ",
                                         "each grid point takes the k subjects nearest in z-scored moderator ",
                                         "space. These provide alternative descriptive views.")
                               ),
                               p("1D trajectories and 2D surfaces can be smoothed for descriptive display. Smoother bands do not account for dependence between overlapping windows. Without a supplied probability, threshold metrics use a logistic standardized-score transform, which is not calibrated clinical risk."),

                               h4("Data sources"),
                               tags$ol(
                                 tags$li(tags$b("Example (synthetic)"), " — loads on startup. ",
                                         "Fully synthetic, no real patient records. Built-in moderation: ",
                                         tags$code("age"), " and ", tags$code("cognitive_score"),
                                         " have imposed signal tendencies with sampling noise; other moderators have no imposed signal."),
                                 tags$li(tags$b("Upload your own (.csv / .rds)"), " — switch the ",
                                         tags$em("Data"), " radio in the sidebar, choose a file, then map ",
                                         "columns to ", tags$em("Outcome"), " (two-class; choose which value ",
                                         "is the transition/positive class), ", tags$em("Score"),
                                         " (continuous), optional ", tags$em("Probability"), " and ",
                                         tags$em("Grouping variable"), ", and one or more ", tags$em("Moderators"),
                                         " (numeric). Check the suggested mapping, then click Update. ",
                                         "Changing data or settings invalidates earlier results and downloads.")
                               ),

                               h4("How to read the plots"),
                               tags$ul(
                                 tags$li("In ", tags$b("1D"), ", the x-axis is the moderator and the ",
                                         "y-axis is the chosen metric per window. Trends are descriptive: overlapping ",
                                         "windows and sampling variability prevent interpreting a trend as a moderation test."),
                                 tags$li("In ", tags$b("2D"), ", the surface colour encodes local performance. ",
                                         "Each metric needs the selected minimum events and non-events after ",
                                         "excluding missing scores (AUC) or probabilities (threshold metrics). ",
                                         "Unsupported values are missing; CSVs retain their eligible counts.")
                               ),

                               h4("Tips"),
                               tags$ul(
                                 tags$li("For broader local estimates, increase range width in fixed-range mode, ",
                                         "window size in fixed-count mode, or neighbourhood size in KNN mode."),
                                 tags$li("Use the ", tags$b("Threshold"), " input only for BACC / Sensitivity / ",
                                         "Specificity — AUC is threshold-free."),
                                 tags$li("Switch tabs above (", tags$b("1D"), " / ", tags$b("2D"),
                                         ") to change analysis mode; the sidebar controls update accordingly."),
                                 tags$li("Both trajectory and grid tables are downloadable from the sidebar.")
                               ),

                               h4("Citation"),
                               p("Cite Prediction Landscape Explorer, Clara Vetter, version ", APP_VERSION, ". ",
                                 tags$a(href = "https://github.com/claravetter/prediction-landscape-explorer", "Software and citation"), "; ",
                                 tags$a(href = "https://github.com/claravetter/pronia-mri-transition-performance-landscapes", "companion analyses"))
                           )
                  ),
                  tabPanel("1D",
                           textOutput("result_status_1d"),
                           div(id = "plot_section_1d",
                               withSpinner(plotOutput("traj_plot_1d", height = "650px"), type = 6),
                               br(),
                               h4("1D summary"),
                               div(class = "tbl-count", textOutput("summary_tbl_1d_count", inline = TRUE)),
                               div(class = "scroll-table", withSpinner(tableOutput("summary_tbl_1d"), type = 6, size = 0.5)),
                               br(),
                               h4("1D window table"),
                               div(class = "tbl-count", textOutput("traj_tbl_1d_count", inline = TRUE)),
                               div(class = "scroll-table", withSpinner(tableOutput("traj_tbl_1d"), type = 6, size = 0.5))
                           )
                  ),
                  tabPanel("2D",
                           textOutput("result_status_2d"),
                           div(id = "plot_section_2d",
                               withSpinner(plotOutput("surface_plot_2d", height = "750px"), type = 6),
                               br(),
                               h4("2D summary"),
                               div(class = "tbl-count", textOutput("summary_tbl_2d_count", inline = TRUE)),
                               div(class = "scroll-table", withSpinner(tableOutput("summary_tbl_2d"), type = 6, size = 0.5)),
                               br(),
                               h4("2D grid table"),
                               div(class = "tbl-count", textOutput("grid_tbl_2d_count", inline = TRUE)),
                               div(class = "scroll-table", withSpinner(tableOutput("grid_tbl_2d"), type = 6, size = 0.5))
                           )
                  )
      )
    )
  )
)

# ================================================================
# Server
# ================================================================
server <- function(input, output, session) {

  # Raw (un-prepped) dataset chosen by the data-source radio.
  raw_rv <- reactiveVal(NULL)
  # Whether the currently loaded raw dataset is the generated synthetic example.
  is_example_rv <- reactiveVal(TRUE)
  # Prepped dataset (post column-mapping) used by the rest of the app.
  dat_rv <- reactiveVal(NULL)

  # ---------------------------
  # Raw loader (driven by the data_source radio + upload)
  # ---------------------------
  observe({
    src <- input$data_source %||% "example"
    # A failed/new upload must never leave the previous dataset usable.
    raw_rv(NULL)
    dat_rv(NULL)
    is_example_rv(identical(src, "example"))
    df <- tryCatch({
      if (src == "example") make_demo_data() else {
        up <- input$data_upload
        req(up)
        if (grepl("\\.rds$", up$name, ignore.case = TRUE)) readRDS(up$datapath)
        else if (grepl("\\.csv$", up$name, ignore.case = TRUE))
          readr::read_csv(up$datapath, show_col_types = FALSE, name_repair = "minimal")
        else stop("Unsupported file type. Use .csv or .rds.")
      }
    }, error = function(e) {
      if (!inherits(e, "shiny.silent.error"))
        showNotification(conditionMessage(e), type = "error", duration = NULL)
      NULL
    })
    if (!is.null(df) && (!is.data.frame(df) || anyDuplicated(names(df)))) {
      showNotification("Upload a data frame with unique column names.", type = "error")
      return()
    }
    raw_rv(df)
  })

  # ---------------------------
  # Column-mapping UI (populated from the raw dataset)
  # ---------------------------
  output$column_mapping_ui <- renderUI({
    df <- raw_rv()
    if (is.null(df)) return(helpText("Loading data…"))

    defaults <- default_col_map_for(df)
    tagList(
      h5("Column mapping"),
      selectInput("col_outcome", "Outcome (two-class)",
                  choices = candidate_outcome_cols(df),
                  selected = defaults$outcome),
      uiOutput("outcome_value_mapping_ui"),
      selectInput("col_score", "Score",
                  choices = candidate_score_cols(df),
                  selected = defaults$score),
      selectInput("col_prob", "Probability (optional)",
                  choices = c("(none)" = "", candidate_prob_cols(df)),
                  selected = defaults$prob),
      selectInput("col_cohort", "Grouping variable (optional)",
                  choices = c("(none)" = "", candidate_cohort_cols(df)),
                  selected = defaults$cohort),
      selectInput("col_moderators", "Moderators (numeric)",
                  choices = candidate_moderator_cols(df),
                  selected = defaults$moderators, multiple = TRUE)
    )
  })

  # Which raw value maps to 0 (non-transition) vs 1 (transition). Depends
  # on the currently-chosen outcome column, so it's a separate uiOutput
  # rather than baked into column_mapping_ui (which only re-renders when
  # the raw dataset itself changes).
  output$outcome_value_mapping_ui <- renderUI({
    df <- raw_rv()
    if (is.null(df) || is.null(input$col_outcome) || !input$col_outcome %in% names(df)) {
      return(NULL)
    }
    vals <- sort(unique(as.character(df[[input$col_outcome]])[!is.na(df[[input$col_outcome]])]))
    if (length(vals) != 2) {
      return(helpText("Selected outcome column does not have exactly two distinct values."))
    }
    np <- default_negative_positive(df, input$col_outcome)
    tagList(
      selectInput("col_outcome_negative", "Value = non-transition (0)",
                  choices = vals, selected = np$negative),
      selectInput("col_outcome_positive", "Value = transition (1)",
                  choices = vals, selected = np$positive)
    )
  })

  apply_col_map <- function(df, col_map) {
    prepped <- tryCatch(prep_data_for_app(df, col_map),
                        error = function(e) {
                          showNotification(paste("Mapping error:", conditionMessage(e)),
                                           type = "error", duration = NULL)
                          NULL
                        })
    dat_rv(prepped)
    if (is.null(prepped)) return()

    if ("cohort" %in% names(prepped)) {
      updateSelectInput(session, "cohorts",
                        choices = levels(prepped$cohort),
                        selected = levels(prepped$cohort))
    } else {
      updateSelectInput(session, "cohorts", choices = character(), selected = character())
    }
    mods <- build_moderator_registry(prepped)
    if (nrow(mods) > 0) {
      ch <- setNames(mods$sort_col,
                     paste0(mods$label, "  [", mods$sort_col, "]"))
      updateSelectInput(session, "moderator_1d",
                        choices = ch, selected = mods$sort_col[1])
      default_x <- mods$sort_col[1]
      default_y <- if (length(mods$sort_col) >= 2) mods$sort_col[2] else mods$sort_col[1]
      updateSelectInput(session, "x_var_2d", choices = ch, selected = default_x)
      updateSelectInput(session, "y_var_2d", choices = ch, selected = default_y)
    }
  }

  # Auto-apply: re-prep whenever the raw dataset changes OR any column
  # mapping dropdown changes. No explicit "Apply" button — the user
  # edits dropdowns; selectors refresh and existing results require Update.
  observe({
    df <- raw_rv()
    dat_rv(NULL)
    req(df)
    # The column-mapping UI must be rendered (i.e. inputs populated)
    # before we read them; outcome + score are required.
    req(input$col_outcome, input$col_score,
        input$col_outcome_negative, input$col_outcome_positive)
    col_map <- list(
      outcome    = input$col_outcome,
      score      = input$col_score,
      prob       = input$col_prob,
      cohort     = input$col_cohort,
      moderators = input$col_moderators,
      negative   = input$col_outcome_negative,
      positive   = input$col_outcome_positive
    )
    apply_col_map(df, col_map)
  })

  # Bootstrap: when the raw dataset first loads (example or upload),
  # the column-mapping UI hasn't rendered yet so input$col_* are NULL.
  # Prepare the detected mapping and selectors. Computation still waits
  # for the corresponding Update button.
  observeEvent(raw_rv(), {
    df <- raw_rv()
    req(df)
    if (is.null(input$col_outcome)) {
      apply_col_map(df, default_col_map_for(df))
    }
  })

  # ---------------------------
  # Synthetic banner (top of page)
  # ---------------------------
  output$synthetic_banner <- renderUI({
    if (!isTRUE(is_example_rv())) return(NULL)
    div(style = "background:#fff4e6;border:1px solid #f0b67f;padding:8px 12px;margin:6px 0;border-radius:4px;",
        strong("Example data — fully synthetic. "),
        "Generated in memory, independently of real participant records. ",
        "Switch to ", em("Upload your own"),
        " to analyse a different dataset.")
  })

  # ---------------------------
  # Cohort selector (optional)
  # ---------------------------
  output$cohort_selector <- renderUI({
    df <- dat_rv()
    if (is.null(df)) return(helpText("Choose a valid column mapping to filter groups."))
    if (!"cohort" %in% names(df)) return(helpText("No grouping variable mapped."))
    selectInput("cohorts", "Groups", choices = levels(df$cohort),
                selected = levels(df$cohort), multiple = TRUE)
  })

  # ---------------------------
  # 1D / 2D moderator selectors
  # ---------------------------
  output$moderator_selector_1d <- renderUI({
    df <- dat_rv()
    if (is.null(df)) return(helpText("Choose a valid column mapping to select a moderator."))
    mods <- build_moderator_registry(df)
    if (nrow(mods) == 0) return(helpText("No moderator columns mapped."))
    selectInput("moderator_1d", "Moderator (x-axis)",
                choices = setNames(mods$sort_col, paste0(mods$label, "  [", mods$sort_col, "]")),
                selected = mods$sort_col[1])
  })

  # Numeric-column selector for the "Mean of a numeric column" histogram option.
  output$hist_mean_col_selector_1d <- renderUI({
    df <- dat_rv()
    if (is.null(df)) return(helpText("Choose a valid column mapping first."))
    # Skip the binary outcome (mean over a bin = transition rate, already a separate option).
    num_cols <- candidate_moderator_cols(df, exclude = "EXP_LABEL")
    if (length(num_cols) == 0) return(helpText("No numeric columns available."))
    default_col <- if ("Mean_Score" %in% num_cols) "Mean_Score" else num_cols[1]
    selectInput("hist_mean_col_1d", "Numeric column",
                choices = num_cols, selected = default_col)
  })

  output$xy_selector_2d <- renderUI({
    df <- dat_rv()
    if (is.null(df)) return(helpText("Choose a valid column mapping to select 2D variables."))
    mods <- build_moderator_registry(df)
    if (nrow(mods) == 0) return(helpText("No moderator columns mapped."))
    choices2 <- setNames(mods$sort_col, paste0(mods$label, "  [", mods$sort_col, "]"))
    default_x <- mods$sort_col[1]
    default_y <- if (length(mods$sort_col) >= 2) mods$sort_col[2] else mods$sort_col[1]
    tagList(
      selectInput("x_var_2d", "X variable", choices = choices2, selected = default_x),
      selectInput("y_var_2d", "Y variable", choices = choices2, selected = default_y)
    )
  })

  # ---------------------------
  # Filtered dataset (cohort filter is optional)
  # ---------------------------
  grouped_data <- reactive({
    df <- dat_rv()
    req(df)
    apply_group_filter(df, input$cohorts)
  })
  data_filtered <- reactive(grouped_data()$data)

  # Capture every setting used by the calculation or presentation. A result
  # is usable only while this request still matches the applied request.
  request_for <- function(dim) {
    reactive({
      df <- data_filtered()
      setting_names <- paste0(names(result_defaults(dim)), "_", dim)
      settings <- setNames(lapply(setting_names, function(n) input[[n]]), names(result_defaults(dim)))
      defaults <- result_defaults(dim)
      for (n in names(defaults)) if (is.null(settings[[n]])) settings[[n]] <- defaults[[n]]
      req(if (dim == "1d") settings$moderator else c(settings$x_var, settings$y_var))
      mapping <- setNames(lapply(c("col_outcome", "col_score", "col_prob", "col_cohort",
                            "col_moderators", "col_outcome_negative", "col_outcome_positive",
                            "data_source"), function(n) input[[n]]),
                          c("outcome", "score", "probability", "group", "moderators", "negative", "positive", "data_source"))
      if (!"cohort" %in% names(df)) mapping["group"] <- list(NULL)
      mapping["selected_groups"] <- list(grouped_data()$selected_groups)
      mapping$group_filter_applied <- grouped_data()$applied
      list(data = df, settings = settings,
           source = if (isTRUE(is_example_rv())) "synthetic_demo" else "uploaded_data",
           # Full raw data and mapping distinguish different uploads/mappings
           # even if the currently selected prepared fields happen to agree.
           data_key = digest::digest(raw_rv(), algo = "sha256"),
           mapping = mapping)
    })
  }
  request_1d <- request_for("1d")
  request_2d <- request_for("2d")
  cached_1d <- eventReactive(input$update_1d, compute_result(request_1d(), "1d"))
  cached_2d <- eventReactive(input$update_2d, compute_result(request_2d(), "2d"))
  current_result <- function(cache, request) reactive({
    r <- cache()
    validate(need(identical(r$request, request()),
                  "Data or settings changed. Click Update to recompute."))
    r
  })
  result_1d <- current_result(cached_1d, request_1d)
  result_2d <- current_result(cached_2d, request_2d)
  prepared_plot_1d <- reactive(prepare_result_plot(result_1d()))
  prepared_plot_2d <- reactive(prepare_result_plot(result_2d()))

  bind_result_ui <- function(dim, current, prepared) {
    status <- reactive({
      tryCatch({
        r <- current()
        if (!prepared()$ok) prepared()$message
        else paste("Current results:", r$request$source, "—", r$created)
      }, error = function(e) "Data or settings changed, or no result is available. Click Update.")
    })
    observe({
      ok <- tryCatch(nrow(current()$table) > 0, error = function(e) FALSE)
      plot_ok <- tryCatch(prepared()$ok, error = function(e) FALSE)
      shinyjs::toggleState(paste0("download_plot_", dim), condition = plot_ok)
      shinyjs::toggleState(paste0("download_csv_", dim), condition = ok)
    })
    output[[paste0("result_status_", dim)]] <- renderText(status())
    output[[if (dim == "1d") "traj_plot_1d" else "surface_plot_2d"]] <- renderPlot({
      p <- prepared()
      validate(need(p$ok, p$message))
      p$plot
    })
    output[[paste0("summary_tbl_", dim)]] <- renderTable(result_summary(current()), digits = 4)
    output[[paste0("summary_tbl_", dim, "_count")]] <- renderText({current(); "1 row (summary)"})
    table_id <- if (dim == "1d") "traj_tbl_1d" else "grid_tbl_2d"
    output[[table_id]] <- renderTable({
      r <- current()
      cols <- if (dim == "1d") c("window_id", "center", "n_obs", "n_events", "n_nonevents", "auc", "bacc", "sens", "spec")
              else c("x_center", "y_center", "n_obs", "n_transition", "auc", "bacc", "sens", "spec")
      prefix <- if (r$request$settings$metric == "auc") "auc_" else "threshold_"
      cols <- c(cols, paste0(prefix, c("n_obs", "n_events", "n_nonevents", "supported")))
      r$table %>% select(any_of(cols))
    }, digits = 4)
    output[[paste0(table_id, "_count")]] <- renderText(paste(nrow(current()$table), "windows/cells"))
    output[[paste0("download_csv_", dim)]] <- downloadHandler(
      filename = function() result_filename(current(), "csv"),
      content = function(file) readr::write_csv(result_export(current()), file)
    )
    output[[paste0("download_plot_", dim)]] <- downloadHandler(
      filename = function() result_filename(current(), "png"),
      content = function(file) {
        p <- prepared()
        validate(need(p$ok, p$message))
        ggplot2::ggsave(file, p$plot, width = 11, height = 9, dpi = 300)
      }
    )
  }
  bind_result_ui("1d", result_1d, prepared_plot_1d)
  bind_result_ui("2d", result_2d, prepared_plot_2d)

}

# ================================================================
# Launch
# ================================================================
shinyApp(ui, server)
