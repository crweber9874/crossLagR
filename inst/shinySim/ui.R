library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(DiagrammeR)

## estimator_choices and dgp_choices are defined in global.R so server.R sees
## them too.

## ---- Academic-grey Bootstrap 5 theme ---------------------------------------
ua_grey_theme <- bs_theme(
  version          = 5,
  bg               = "#FAFBFC",
  fg               = "#1F2937",
  primary          = "#37474F",
  secondary        = "#607D8B",
  success          = "#2E7D32",
  info             = "#455A64",
  warning          = "#B7791F",
  danger           = "#C62828",
  base_font        = font_google("Source Sans 3"),
  heading_font     = font_google("Source Sans 3"),
  code_font        = font_google("JetBrains Mono"),
  "navbar-bg"          = "#263238",
  "navbar-light-color" = "#ECEFF1",
  "border-radius"      = "8px",
  "card-border-color"  = "#E5E7EB",
  "card-cap-bg"        = "#F1F3F5",
  "body-bg"            = "#FAFBFC",
  "body-color"         = "#1F2937"
)

extra_css <- HTML("
  .navbar-brand { font-weight: 700; letter-spacing: 0.3px; }
  .navbar-nav .nav-link { color: #ECEFF1 !important; font-weight: 500; }
  .navbar-nav .nav-link.active, .navbar-nav .nav-link:hover {
    color: #FFFFFF !important; border-bottom: 2px solid #B0BEC5;
  }
  .card { box-shadow: 0 1px 2px rgba(15,23,42,0.06); }
  .card-header { font-weight: 600; color: #1F2937; border-bottom: 1px solid #E5E7EB; }
  .btn { font-weight: 600 !important; letter-spacing: 0.2px; }
  .btn-primary, .btn-info, .btn-success, .btn-warning, .btn-danger {
    color: #FFFFFF !important;
  }
  .btn-primary:hover, .btn-info:hover, .btn-success:hover,
  .btn-warning:hover, .btn-danger:hover {
    color: #FFFFFF !important; filter: brightness(0.92);
  }
  h3, h4, h5 { color: #263238; font-weight: 600; }
  .small-note { color: #6C7280; font-size: 0.88em; font-style: italic; }
  .dag-frame {
    background: #FFFFFF; border: 1px solid #E5E7EB; border-radius: 8px;
    padding: 14px; min-height: 360px;
  }
  pre {
    background: #F1F3F5; border-left: 3px solid #37474F;
    color: #1F2937; padding: 10px 12px; border-radius: 4px;
  }
  .mc-log-box {
    background: #1F2937; color: #ECEFF1;
    font-family: 'JetBrains Mono', 'SF Mono', 'Menlo', monospace;
    font-size: 0.83em;
    padding: 10px 12px; border-radius: 6px;
    height: 200px; overflow-y: auto;
    white-space: pre-wrap; word-break: break-word;
    border: 1px solid #37474F;
  }
  .mc-log-box pre {
    background: transparent; border: 0; color: #ECEFF1; padding: 0;
  }
  /* Compact sidebar inputs */
  .sim-sidebar .form-group { margin-bottom: 8px; }
  .sim-sidebar label { font-weight: 600; font-size: 0.88em; margin-bottom: 2px; color: #455A64; }
  .sim-sidebar .form-control { font-size: 0.9em; padding: 4px 8px; height: 32px; }
  .sim-sidebar h5 {
    border-bottom: 1px solid #E5E7EB; padding-bottom: 4px;
    margin-top: 14px; margin-bottom: 8px;
    font-size: 0.92em; text-transform: uppercase;
    letter-spacing: 0.6px; color: #607D8B;
  }
  .sim-sidebar .btn { padding: 8px 14px; }
  .run-btn-row { display: flex; gap: 8px; margin-top: 10px; }
  .run-btn-row .btn { flex: 1; }
  .status-banner {
    background: #FFFFFF; border: 1px solid #E5E7EB; border-radius: 8px;
    padding: 12px 16px; margin-bottom: 12px;
  }
  .status-banner .status-text { font-weight: 500; color: #263238; }
  .dataTables_wrapper { font-size: 0.92em; }
  .nav-tabs .nav-link.active {
    border-bottom: 2px solid #37474F; color: #37474F; font-weight: 600;
  }
  /* Kill busy / 'recalculating' greying everywhere.
     Without these, every output that re-renders fades to 50% opacity
     for a few hundred ms — combined across many panels that reads as
     'the screen is permanently grey'. */
  .recalculating,
  .recalculating *,
  .shiny-output-error,
  .shiny-output-error-validation,
  .shiny-bound-output { opacity: 1 !important; filter: none !important; }
  .shiny-busy-indicators-spinner, .shiny-busy-indicators-pulse,
  html.shiny-busy::after, .shiny-busy::after,
  .shiny-busy-indicator-container { display: none !important; }
  body.shiny-busy, html.shiny-busy { cursor: default !important; }
  /* bslib >= 0.7 puts a transient overlay on cards while recalculating;
     also suppress that. */
  .bslib-page-fill .recalculating { opacity: 1 !important; }
  .card.recalculating, .card .recalculating { opacity: 1 !important; }
")

## Compact card wrapper
acad_card <- function(title, ...) {
  card(full_screen = FALSE, card_header(title), card_body(...))
}

## A two-up numeric input row (label P / label Q)
ni_row <- function(id_p, id_q, label_p, label_q, value, step = 0.05) {
  tags$div(style = "display:flex; gap:8px;",
    tags$div(style = "flex:1;",
      numericInput(id_p, label_p, value = value, step = step)),
    tags$div(style = "flex:1;",
      numericInput(id_q, label_q, value = value, step = step))
  )
}

## A "label  min  max  step" row for the dynamic-grid spec.
## Creates inputs <id>_min, <id>_max, <id>_step.
dyn_range_row <- function(id, label, min_v, max_v, step_v) {
  tags$div(style = "display:flex; gap:4px; align-items:center; margin-bottom:4px;",
    tags$div(style = "flex:0 0 65px; font-weight:600; font-size:0.85em; color:#37474F;",
             label),
    tags$div(style = "flex:1;",
      numericInput(paste0(id, "_min"),  "min",  value = min_v,  step = 0.05)),
    tags$div(style = "flex:1;",
      numericInput(paste0(id, "_max"),  "max",  value = max_v,  step = 0.05)),
    tags$div(style = "flex:1;",
      numericInput(paste0(id, "_step"), "step", value = step_v, step = 0.05, min = 0.01))
  )
}

ui <- page_navbar(
  id      = "tabs",
  title   = "CrossLagR: Simulations",
  theme   = ua_grey_theme,
  bg      = "#263238",
  inverse = TRUE,
  header  = tags$head(tags$style(extra_css)),
  selected = "Simulate",

  ## ===========================================================================
  ## MODEL — diagram + lavaan syntax
  ## ===========================================================================
  nav_panel(
    "Model",
    layout_columns(
      col_widths = c(4, 8),
      div(
        acad_card("Model selection",
          selectInput("estimator", "Estimator:",
                      choices = estimator_choices, selected = "CLPM"),
          numericInput("waves", "Waves (T):", value = 5, min = 2, max = 10),
          tags$p(class = "small-note",
                 "DAG redraws at the selected wave count (clipped to 2-8).")
        ),
        br(),
        acad_card("Quick constraints",
          checkboxInput("constrain_beta",   "Constrain AR across waves",         value = TRUE),
          checkboxInput("constrain_omega",  "Constrain CL across waves",         value = TRUE),
          checkboxInput("constrain_resvar", "Constrain residual variances",      value = TRUE),
          checkboxInput("constrain_rescov", "Constrain residual covariances",    value = TRUE),
          checkboxInput("estimate_means",   "Estimate means at wave 1",          value = TRUE),
          checkboxInput("start_values",     "Use starting values",               value = FALSE),
          conditionalPanel(
            "input.estimator == 'LCHANGE'",
            checkboxInput("lcs_constant_change",  "Estimate constant-change A", value = TRUE),
            checkboxInput("lcs_change_to_change", "Estimate change-to-change",  value = FALSE),
            selectInput("lchange_variable_type", "Variables:",
                        choices = c("bivariate", "univariate"), selected = "bivariate")
          ),
          conditionalPanel(
            "input.estimator == 'LGM' || input.estimator == 'LCMSR'",
            selectInput("lgm_variable_type", "Variables:",
                        choices = c("bivariate", "univariate"), selected = "bivariate"),
            checkboxInput("lgm_quadratic", "Include quadratic factor", value = FALSE)
          ),
          conditionalPanel(
            "input.estimator == 'BB'",
            selectInput("bb_x_effect", "X effect on Y:",
                        choices = c("lagged", "concurrent", "none"), selected = "lagged"),
            checkboxInput("bb_x_autoreg",      "X autoregression",       value = TRUE),
            checkboxInput("bb_y_on_x",         "Y effect on X",          value = TRUE),
            ## Default ON so BB's lavaan syntax uses the unified labels
            ## (ar_x, ar_y, cl_xy, cl_yx) instead of per-wave forms
            ## (ar_x2, ar_x3, ...). Without this, monteCarloLavaan can't
            ## extract BB's estimates and the estimator silently disappears
            ## from compare-mode plots.
            checkboxInput("bb_constrain_coef", "Constrain coefficients", value = TRUE)
          ),
          conditionalPanel(
            "input.estimator == 'TSO'",
            checkboxInput("tso_constrain_state_var", "Constrain state variances",  value = TRUE),
            checkboxInput("tso_constrain_state_cov", "Constrain state covariances", value = TRUE)
          )
        )
      ),
      div(
        acad_card(textOutput("dag_title"),
          div(class = "dag-frame",
              DiagrammeR::grVizOutput("model_dag", height = "440px")),
          tags$p(class = "small-note", style = "margin-top: 10px;",
                 "Solid: structural paths (AR blue, CL red). ",
                 "Dashed: covariances. Circles: latent. Boxes: observed.")
        ),
        br(),
        acad_card("lavaan syntax", verbatimTextOutput("model_syntax"))
      )
    )
  ),

  ## ===========================================================================
  ## SIMULATE — single-page workflow: inputs in sidebar, results in main pane
  ## ===========================================================================
  nav_panel(
    "Simulate",

    layout_sidebar(
      sidebar = sidebar(
        width = 340, open = "always", position = "left",
        class = "sim-sidebar",
        title = "Simulation setup",

        ## ===== Step 1 — choose run type =================================
        h5("1. Run type"),
        radioButtons("grid_mode", NULL,
                     choices = c(
                       "Run Single Parameterization (one estimator)"      = "single",
                       "Run Varying Parameterizations (one estimator)"    = "dynamic",
                       "Run Varying Parameterizations across Estimators"  = "compare"),
                     selected = "single"),

        ## ===== Step 2 — choose estimator(s) =============================
        ## Single dropdown for single/dynamic modes; multi-select for compare.
        h5("2. Estimator(s)"),
        conditionalPanel(
          "input.grid_mode != 'compare'",
          selectInput("sim_estimator", NULL,
                      choices = estimator_choices, selected = "CLPM")
        ),
        conditionalPanel(
          "input.grid_mode == 'compare'",
          tags$p(class = "small-note",
                 "All checked estimators are fit to every cell."),
          checkboxGroupInput(
            "compare_estimators", NULL,
            choices = c("CLPM"        = "CLPM",
                        "RI-CLPM"     = "RICLPM",
                        "ALT"         = "ALT",
                        "LGM"         = "LGM",
                        "LCM-SR"      = "LCMSR",
                        "LCS"         = "LCHANGE",
                        "Bollen-Brand" = "BB",
                        "TSO"         = "TSO"),
            selected = c("CLPM", "RICLPM", "ALT", "LCMSR", "BB", "LCHANGE"),
            inline   = TRUE)
        ),

        ## ===== Step 3 — DGP =============================================
        h5("3. Data-generating process"),
        selectInput("data_generation", NULL,
                    choices = dgp_choices, selected = "CLPM"),
        numericInput("sim_waves", "Waves (T):", value = 5, min = 2, max = 10),

        ## DGP-specific knobs (only the confounder DGP has any currently).
        conditionalPanel(
          "input.data_generation == 'clpmu'",
          tags$div(style = "background:#FFF8E1; border:1px solid #FFE082;
                           border-radius:6px; padding:6px 10px; margin-bottom:8px;",
            tags$strong("CLPM + unmeasured confounder knobs"),
            ni_row("confounder_p", "confounder_q",
                   "effect on X", "effect on Y", value = 0.3, step = 0.05),
            ni_row("confounder_variance", "confounder_stability",
                   "variance", "stability", value = 1, step = 0.1)
          )
        ),

        ## ===== Step 4 — parameters ======================================
        ## Single mode  → six scalar inputs (one value per parameter).
        ## Varying mode → six range rows (min/max/step), mutually exclusive
        ##                with the single-mode inputs.
        h5("4. Parameters"),
        conditionalPanel(
          "input.grid_mode == 'single'",
          ni_row("stability_p", "stability_q",
                 "ar_x (stab. P)", "ar_y (stab. Q)", value = 0.30, step = 0.05),
          ni_row("cross_p", "cross_q",
                 "cl_yx (Y→X)", "cl_xy (X→Y)",       value = 0.10, step = 0.05),
          ni_row("variance_p", "variance_q",
                 "within var P", "within var Q",      value = 0.50, step = 0.1),
          ni_row("icc_x", "icc_y",
                 "ICC X (between/total)", "ICC Y (between/total)",
                 value = 0.50, step = 0.05),
          tags$div(class = "small-note",
                   style = "margin: -4px 0 8px 0; font-size: 0.78em;",
                   textOutput("derived_between_text", inline = TRUE)),
          numericInput("cov_pq", "within-time cov(p,q):", value = 0, step = 0.05)
        ),
        conditionalPanel(
          "input.grid_mode == 'dynamic' || input.grid_mode == 'compare'",
          tags$p(class = "small-note",
                 "For each parameter set min / max / step. Set min = max to ",
                 "hold the parameter fixed. The full grid is the cross-product. ",
                 "(The single-mode scalar inputs above are still bound and ",
                 "serve as fall-back defaults when a row's min/max are blank.)"),
          dyn_range_row("dyn_icc_x", "ICC X",  0.00, 0.70, 0.10),
          dyn_range_row("dyn_icc_y", "ICC Y",  0.00, 0.70, 0.10),
          dyn_range_row("dyn_ar_x",  "ar_x",   0.30, 0.30, 0.05),
          dyn_range_row("dyn_ar_y",  "ar_y",   0.30, 0.30, 0.05),
          dyn_range_row("dyn_cl_xy", "cl_xy",  0.10, 0.10, 0.05),
          dyn_range_row("dyn_cl_yx", "cl_yx",  0.10, 0.10, 0.05),
          textOutput("dyn_grid_preview") |>
            tagAppendAttributes(class = "small-note",
                                style = "margin-top: 6px; display:block;")
        ),

        ## ===== Step 5 — sample ==========================================
        h5("5. Sample"),
        tags$div(style = "display:flex; gap:8px;",
          tags$div(style = "flex:1;",
            numericInput("trials", "Trials/cell:", value = 50, min = 1, max = 10000)),
          tags$div(style = "flex:1;",
            numericInput("sample_size", "N:", value = 1000, min = 50, max = 50000))
        ),
        tags$div(style = "display:flex; gap:8px;",
          tags$div(style = "flex:1;",
            numericInput("seed", "Seed:", value = 1, min = 1)),
          tags$div(style = "flex:1; padding-top: 22px;",
            checkboxInput("verbose", "Verbose log", value = TRUE))
        ),

        tags$hr(),
        ## Always-visible model-ready indicator. The model is auto-prepared
        ## from the selected estimator + constraints (no manual init needed).
        div(style = paste("background:#F1F8E9; border:1px solid #C5E1A5;",
                          "border-radius:6px; padding:6px 10px; margin-bottom:8px;",
                          "font-size:0.85em; color:#33691E;"),
            tags$span(style = "font-weight:600;", "✓ Model ready: "),
            textOutput("ready_text", inline = TRUE)
        ),
        actionButton("preview_model", "Preview model syntax",
                     class = "btn-secondary w-100",
                     style = "margin-bottom: 8px;"),
        div(class = "run-btn-row",
          actionButton("run_simulation", "▶ Run Monte Carlo",
                       class = "btn-primary"),
          actionButton("abort_simulation", "Stop",
                       class = "btn-danger")
        ),
        actionButton("run_single_fit", "Single fit (one dataset)",
                     class = "btn-info w-100",
                     style = "margin-top: 8px;"),
        downloadButton("download_results", "Download CSV",
                       class = "btn-success w-100",
                       style = "margin-top: 8px;"),
        textOutput("grid_size_info") |>
          tagAppendAttributes(class = "small-note",
                              style = "margin-top: 10px; display: block;")
      ),

      ## ---- Main area ---------------------------------------------------------
      div(
        ## Always-visible status banner + progress bar
        div(class = "status-banner",
          shinyWidgets::progressBar(
            id = "sim_progress", value = 0, total = 100,
            display_pct = TRUE, status = "primary", title = "Idle"
          ),
          tags$div(class = "status-text", style = "margin-top: 6px;",
                   textOutput("sim_status_text"))
        ),

        ## Tabbed main area for: Log | Summary | Bias | Distributions | Plotly | Single fit
        navset_card_tab(
          id = "result_tabs",

          nav_panel("Live log",
            tags$p(class = "small-note",
                   "Streamed from the background R worker. Updates ~700 ms."),
            div(class = "mc-log-box",
                verbatimTextOutput("mc_log", placeholder = FALSE))
          ),

          nav_panel("Summary",
            tags$p(class = "small-note",
                   "Per parameter cell: mean estimate, bias vs DGP truth, RMSE, ",
                   "mean fit indices, convergence rate."),
            DT::dataTableOutput("mc_summary")
          ),

          nav_panel("Bias plot",
            layout_columns(
              col_widths = c(3, 9),
              div(
                selectInput("bias_param", "Parameter:",
                            choices = c("ar_x", "ar_y", "cl_xy", "cl_yx"),
                            selected = "cl_xy"),
                selectInput("bias_x", "X-axis (DGP knob):", choices = NULL)
              ),
              plotOutput("bias_ggplot", height = "440px")
            )
          ),

          ## Reconstructs chapter-5 fig-sim1-bias / fig-sim1-rmse: one panel
          ## per AR/CL parameter, one line per estimator, x-axis is a swept
          ## DGP knob (default ICC X).
          nav_panel("Estimator comparison",
            layout_columns(
              col_widths = c(3, 9),
              div(
                radioButtons("cmp_view", "View:",
                             choices = c(
                               "Bias / RMSE trajectory"        = "traj",
                               "Sampling distributions (ridges)" = "ridges",
                               "Performance heatmap"           = "heat",
                               "Convergence rate"              = "conv",
                               "Side-by-side table"            = "table"),
                             selected = "traj"),
                tags$hr(style = "margin: 8px 0;"),

                conditionalPanel(
                  "input.cmp_view == 'traj'",
                  selectInput("cmp_metric", "Metric:",
                              choices = c("Bias"  = "bias",
                                          "RMSE"  = "rmse"),
                              selected = "bias"),
                  selectInput("cmp_x", "X-axis (DGP knob):", choices = NULL),
                  checkboxInput("cmp_show_ci",
                                "Show 95% Monte-Carlo intervals", value = TRUE),
                  checkboxInput("cmp_drop_outliers",
                                "Clip y-axis to ±0.3 (drop CLPM ar_x outliers)",
                                value = FALSE),
                  tags$p(class = "small-note",
                         "Facets = AR/CL parameters; lines = estimators. ",
                         "Dashed line at zero is unbiased. Whiskers are mean ± ",
                         "1.96 · SD/√n_trials.")
                ),

                conditionalPanel(
                  "input.cmp_view == 'ridges'",
                  selectInput("cmp_ridges_y", "Stack ridges by (DGP knob):",
                              choices = NULL),
                  selectInput("cmp_ridges_facet", "Facet parameter:",
                              choices = c("All four (grid)" = "all",
                                          "ar_x"  = "ar_x",
                                          "ar_y"  = "ar_y",
                                          "cl_xy" = "cl_xy",
                                          "cl_yx" = "cl_yx"),
                              selected = "all"),
                  sliderInput("cmp_ridges_clip",
                              "Clip estimates to ± of truth:",
                              min = 0.5, max = 5, value = 2, step = 0.5),
                  tags$p(class = "small-note",
                         "Rows = estimators, columns = parameters. Each ridge ",
                         "is the sampling distribution at one cell; vertical ",
                         "line marks the true value. Color encodes mean |est − truth| ",
                         "(green = on, red = off).")
                ),

                conditionalPanel(
                  "input.cmp_view == 'heat'",
                  selectInput("cmp_heat_metric", "Cell value:",
                              choices = c("Bias"      = "bias",
                                          "|Bias|"    = "abs_bias",
                                          "RMSE"      = "rmse"),
                              selected = "abs_bias"),
                  checkboxInput("cmp_heat_label", "Show numeric labels", value = TRUE),
                  tags$p(class = "small-note",
                         "Rows = estimators, columns = parameters × cells. ",
                         "Fill encodes the chosen metric; green = good, red = bad.")
                ),

                conditionalPanel(
                  "input.cmp_view == 'conv'",
                  selectInput("cmp_conv_x", "X-axis (DGP knob):", choices = NULL),
                  tags$p(class = "small-note",
                         "Convergence rate per estimator at each cell. ",
                         "Cells below 100% indicate fits failed; bias / RMSE ",
                         "from those cells are based on partial data and ",
                         "should be interpreted cautiously.")
                ),

                conditionalPanel(
                  "input.cmp_view == 'table'",
                  tags$p(class = "small-note",
                         "One row per (estimator × cell). Columns: mean, bias, ",
                         "RMSE for each unified-label parameter. Sortable & scrollable.")
                )
              ),
              div(
                conditionalPanel("input.cmp_view == 'traj'",
                                  plotOutput("compare_ggplot", height = "540px")),
                conditionalPanel("input.cmp_view == 'ridges'",
                                  plotOutput("compare_ridges", height = "640px")),
                conditionalPanel("input.cmp_view == 'heat'",
                                  plotOutput("compare_heatmap", height = "540px")),
                conditionalPanel("input.cmp_view == 'conv'",
                                  plotOutput("compare_conv", height = "440px")),
                conditionalPanel("input.cmp_view == 'table'",
                                  DT::dataTableOutput("compare_table"))
              )
            )
          ),

          nav_panel("Distributions",
            layout_columns(
              col_widths = c(3, 9),
              div(
                selectInput("dist_param", "Parameter:",
                            choices = c("ar_x", "ar_y", "cl_xy", "cl_yx"),
                            selected = "cl_xy"),
                selectInput("dist_type", "Type:",
                            choices = c("Density" = "density",
                                        "Histogram" = "hist",
                                        "Ridges (by cell)" = "ridges"),
                            selected = "density"),
                checkboxInput("dist_show_truth", "Mark true value", value = TRUE)
              ),
              plotOutput("dist_ggplot", height = "440px")
            )
          ),

          nav_panel("Interactive (plotly)",
            layout_columns(
              col_widths = c(3, 9),
              div(
                selectInput("plot_var_x", "X-axis:", choices = NULL),
                selectInput("plot_var_y", "Y-axis:", choices = NULL),
                selectInput("plot_var_color", "Color:", choices = NULL),
                selectInput("plot_var_facet", "Facet by:", choices = NULL)
              ),
              plotlyOutput("results_plot", height = "440px")
            )
          ),

          nav_panel("Single fit",
            tags$p(class = "small-note",
                   "One dataset, one lavaan fit. Diagnostic detail for the chosen estimator."),
            verbatimTextOutput("single_fit_header"),
            layout_columns(col_widths = c(6, 6),
              acad_card("Fit indices",  DT::dataTableOutput("single_fit_indices")),
              acad_card("Key estimates", DT::dataTableOutput("single_fit_keys"))
            ),
            acad_card("Full parameter table (standardized)",
                      DT::dataTableOutput("single_fit_params")),
            acad_card("Raw lavaan summary",
                      verbatimTextOutput("single_fit_lavaan_summary"))
          ),

          nav_panel("Per-trial rows",
            DT::dataTableOutput("results_table")
          )
        )
      )
    )
  ),

  ## ===========================================================================
  ## DOCUMENTATION
  ## ===========================================================================
  nav_panel(
    "Documentation",
    layout_columns(
      col_widths = 12,
      acad_card("Workflow",
        tags$ol(
          tags$li("Open the ", strong("Model"), " tab to inspect the DAG and lavaan syntax for the chosen estimator."),
          tags$li("Open ", strong("Simulate"), ". Pick an estimator (the model you'll fit), a DGP (the model that generates the data), set the true parameter values + sample size."),
          tags$li("Optionally enable ", strong("Grid sweep"), " to vary one parameter across a range — each value becomes a separate cell in the MC."),
          tags$li("Click ", strong("Run Monte Carlo"), ". Progress bar and live log update as cells complete; click ", strong("Stop"), " to abort."),
          tags$li("Results populate the ", strong("Summary"), ", ", strong("Bias plot"), ", ", strong("Distributions"), ", and ", strong("Interactive"), " sub-tabs as soon as the run finishes.")
        )
      ),
      acad_card("Estimators",
        tags$ul(
          tags$li(strong("CLPM:"), " classic cross-lagged panel"),
          tags$li(strong("RI-CLPM:"), " random-intercept CLPM"),
          tags$li(strong("ALT:"), " Autoregressive Latent Trajectory; growth + AR/CL on latent scores"),
          tags$li(strong("LGM / LCM-SR:"), " latent growth, optionally with structured AR/CL residuals"),
          tags$li(strong("LCS:"), " bivariate latent change score"),
          tags$li(strong("Bollen-Brand:"), " dynamic panel with latent fixed effects"),
          tags$li(strong("TSO:"), " Trait-State-Occasion"),
          tags$li(strong("CTSEM, RI, OLS:"), " runtime estimators (MC only, no static syntax)")
        )
      ),
      acad_card("DGPs",
        p(strong("Parameterized simulators"),
          " — simCLPM, simRICLPM, simCLPMu, simLChange — use native parameterizations."),
        p(strong("Estimator-as-DGP"),
          " populates the chosen estimator's lavaan syntax with the DGP true ",
          "ar/cl values and calls lavaan::simulateData(). Lets you fit model A ",
          "to data generated by model B and measure the resulting bias.")
      )
    )
  )
)
