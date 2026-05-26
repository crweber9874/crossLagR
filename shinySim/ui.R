library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)

ui <- dashboardPage(
  dashboardHeader(title = "Monte Carlo CLPM Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Setup",        tabName = "setup",   icon = icon("cog")),
      menuItem("Results",      tabName = "results", icon = icon("chart-line")),
      menuItem("Documentation", tabName = "docs",   icon = icon("book"))
    )
  ),

  dashboardBody(
    useShinyjs(),

    tabItems(
      ## ---- Setup ----
      tabItem(
        tabName = "setup",

        fluidRow(
          tabBox(
            id = "setup_tabs", width = 12,

            ## 1. Simulation configuration
            tabPanel(
              "1. Simulation",
              fluidRow(
                column(6,
                  selectInput("estimator", "Estimator:",
                              choices  = c("OLS", "RICLPM", "CLPM", "CTSEM",
                                           "FI", "LCHANGE", "RI"),
                              selected = "CLPM"),
                  selectInput("data_generation", "Data-generating process:",
                              choices  = list("CLPM" = "clpm",
                                              "RI-CLPM" = "riclpm"),
                              selected = "clpm"),
                  conditionalPanel(
                    "input.estimator == 'RICLPM'",
                    selectInput("riclpm_type", "RI-CLPM variant:",
                                choices = c("riclpm", "riclpm_nolag"),
                                selected = "riclpm")
                  ),
                  conditionalPanel(
                    "input.estimator == 'LCHANGE'",
                    selectInput("lchange_type", "Latent-change type:",
                                choices = c("dual_change", "constant_change",
                                            "proportional_change"),
                                selected = "dual_change")
                  )
                ),
                column(6,
                  numericInput("trials", "Trials per cell:", value = 100, min = 1, max = 10000),
                  numericInput("waves",  "Waves:",            value = 3,   min = 2, max = 10),
                  numericInput("sample_size", "Sample size:", value = 1000, min = 50, max = 10000),
                  checkboxInput("verbose", "Verbose log to console", value = FALSE)
                )
              )
            ),

            ## 2. DGP parameters
            tabPanel(
              "2. DGP parameters",
              helpText("These are the population values used to simulate data."),

              h5("Stability (autoregressive)"),
              fluidRow(
                column(6,
                  numericInput("stability_p_min",  "Stability P min:",  value = 0.2, step = 0.1, min = 0, max = 1),
                  numericInput("stability_p_max",  "Stability P max:",  value = 0.3, step = 0.1, min = 0, max = 1),
                  numericInput("stability_p_step", "Step:",             value = 0.1, step = 0.05, min = 0.01, max = 0.5)
                ),
                column(6,
                  numericInput("stability_q_min",  "Stability Q min:",  value = 0.2, step = 0.1, min = 0, max = 1),
                  numericInput("stability_q_max",  "Stability Q max:",  value = 0.3, step = 0.1, min = 0, max = 1),
                  numericInput("stability_q_step", "Step:",             value = 0.1, step = 0.05, min = 0.01, max = 0.5)
                )
              ),

              h5("Cross-lagged"),
              fluidRow(
                column(6,
                  numericInput("cross_p_min",  "Cross P min:",  value = 0.0, step = 0.1, min = -1, max = 1),
                  numericInput("cross_p_max",  "Cross P max:",  value = 0.1, step = 0.1, min = -1, max = 1),
                  numericInput("cross_p_step", "Step:",         value = 0.1, step = 0.05, min = 0.01, max = 0.5)
                ),
                column(6,
                  numericInput("cross_q_min",  "Cross Q min:",  value = 0.0, step = 0.1, min = -1, max = 1),
                  numericInput("cross_q_max",  "Cross Q max:",  value = 0.1, step = 0.1, min = -1, max = 1),
                  numericInput("cross_q_step", "Step:",         value = 0.1, step = 0.05, min = 0.01, max = 0.5)
                )
              ),

              h5("Variances and covariance"),
              fluidRow(
                column(6,
                  numericInput("variance_p", "Within-var P:", value = 0.5, step = 0.1, min = 0.01, max = 5),
                  numericInput("variance_q", "Within-var Q:", value = 0.5, step = 0.1, min = 0.01, max = 5)
                ),
                column(6,
                  numericInput("variance_between_x", "Between-var X:", value = 0.8, step = 0.1, min = 0.01, max = 5),
                  numericInput("variance_between_y", "Between-var Y:", value = 0.8, step = 0.1, min = 0.01, max = 5)
                )
              ),
              numericInput("cov_pq", "Within-time covariance (p,q):", value = 0, step = 0.1, min = -2, max = 2)
            ),

            ## 3. Estimator parameters
            tabPanel(
              "3. Estimator parameters",
              helpText("Map directly to crossLagR::estimateCLPM() / estimateRICLPM() arguments."),
              fluidRow(
                column(6,
                  checkboxInput("constrain_beta",   "Constrain AR (beta) across waves",   value = TRUE),
                  checkboxInput("constrain_omega",  "Constrain CL (omega) across waves",  value = TRUE),
                  checkboxInput("constrain_resvar", "Constrain residual variances",       value = TRUE),
                  checkboxInput("constrain_rescov", "Constrain residual covariances",     value = TRUE)
                ),
                column(6,
                  checkboxInput("estimate_means", "Estimate means at wave 1", value = TRUE),
                  checkboxInput("start_values",   "Use starting values",      value = FALSE)
                )
              )
            ),

            ## 4. Code preview
            tabPanel(
              "4. Code preview",
              h5("lavaan model syntax for the chosen estimator"),
              verbatimTextOutput("model_syntax"),
              h5("Reproducible R call"),
              verbatimTextOutput("r_code")
            )
          )
        ),

        ## Grid + Run/Status row
        fluidRow(
          box(title = "Parameter grid preview", status = "success",
              solidHeader = TRUE, width = 8,
              DT::dataTableOutput("param_preview"),
              br(),
              textOutput("grid_size_info")),

          box(title = "Run simulation", status = "success",
              solidHeader = TRUE, width = 4,
              actionButton("initialize_models", "Initialize models",
                           class = "btn-warning btn-lg",
                           style = "width: 100%; margin-bottom: 12px;"),
              actionButton("run_simulation", "Run Monte Carlo simulation",
                           class = "btn-primary btn-lg",
                           style = "width: 100%; margin-bottom: 12px;"),
              shinyWidgets::progressBar(
                id = "sim_progress", value = 0, total = 100,
                display_pct = TRUE, status = "primary",
                title = "Idle"
              ),
              tags$div(id = "simulation_status",
                       style = "text-align: center; margin-top: 8px;",
                       textOutput("sim_status_text")),
              br(),
              downloadButton("download_results", "Download results",
                             class = "btn-success",
                             style = "width: 100%;"))
        )
      ),

      ## ---- Results ----
      tabItem(
        tabName = "results",
        fluidRow(
          box(title = "Simulation results", status = "primary",
              solidHeader = TRUE, width = 12,
              DT::dataTableOutput("results_table"))
        ),
        fluidRow(
          box(title = "Visualization", status = "info",
              solidHeader = TRUE, width = 12,
              selectInput("plot_var_x",     "X-axis:", choices = NULL),
              selectInput("plot_var_y",     "Y-axis:", choices = NULL),
              selectInput("plot_var_color", "Color:",  choices = NULL),
              plotlyOutput("results_plot", height = "500px"))
        )
      ),

      ## ---- Documentation ----
      tabItem(
        tabName = "docs",
        fluidRow(
          box(title = "Notes", status = "info", solidHeader = TRUE, width = 12,
              h3("Estimators"),
              tags$ul(
                tags$li(strong("OLS:"), "wide-format OLS regression"),
                tags$li(strong("CLPM:"), "lavaan CLPM via crossLagR::estimateCLPM"),
                tags$li(strong("RICLPM:"), "RI-CLPM (with riclpm_type variant)"),
                tags$li(strong("CTSEM:"), "continuous-time SEM"),
                tags$li(strong("FI / RI / LCHANGE:"), "fixed-intercept, random-intercept, and latent-change variants")
              ),
              h3("DGPs"),
              tags$ul(
                tags$li(strong("clpm:"), "data simulated from a CLPM"),
                tags$li(strong("riclpm:"), "data simulated from an RI-CLPM with stable trait factors")
              ),
              h3("Output"),
              p("Each row is one parameter cell × estimator estimate. Look for bias (estimate vs. true value), coverage, and SE."))
        )
      )
    )
  )
)
