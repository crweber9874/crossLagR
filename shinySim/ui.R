library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Monte Carlo Cross-Lagged Panel Model Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Simulation Setup", tabName = "setup", icon = icon("cog")),
      menuItem("Results", tabName = "results", icon = icon("chart-line")),
      menuItem("Documentation", tabName = "docs", icon = icon("book"))
    )
  ),

  dashboardBody(
    tabItems(
      # Setup Tab
      tabItem(tabName = "setup",
              fluidRow(
                # Simulation Configuration
                box(title = "Simulation Configuration", status = "primary", solidHeader = TRUE, width = 6,
                    selectInput("estimator", "Estimator Method:",
                                choices = c("OLS", "RICLPM", "CLPM", "CTSEM"),
                                selected = "OLS"),

                    selectInput("data_generation", "Data Generation Process:",
                                choices = list(
                                  "Cross-Lagged Panel Model" = "clpm",
                                  "Random Intercept CLPM" = "riclpm",
                                  "CLPM with Confounder" = "clpm_confounder"
                                ),
                                selected = "clpm"),

                    numericInput("trials", "Number of Trials:", value = 100, min = 1, max = 10000),
                    numericInput("waves", "Number of Waves:", value = 3, min = 2, max = 10),
                    numericInput("sample_size", "Sample Size:", value = 1000, min = 50, max = 10000),

                    checkboxInput("verbose", "Verbose Output", value = TRUE)
                ),

                # Parameter Grid Configuration
                box(title = "Parameter Configuration", status = "info", solidHeader = TRUE, width = 6,
                    h4("Core Parameters"),

                    # Stability parameters
                    div(style = "margin-bottom: 15px;",
                        h5("Stability Parameters"),
                        fluidRow(
                          column(6,
                                 numericInput("stability_p_min", "Stability P (min):", value = 0.2, step = 0.1, min = 0, max = 1),
                                 numericInput("stability_p_max", "Stability P (max):", value = 0.3, step = 0.1, min = 0, max = 1),
                                 numericInput("stability_p_step", "Step:", value = 0.1, step = 0.05, min = 0.01, max = 0.5)
                          ),
                          column(6,
                                 numericInput("stability_q_min", "Stability Q (min):", value = 0.2, step = 0.1, min = 0, max = 1),
                                 numericInput("stability_q_max", "Stability Q (max):", value = 0.3, step = 0.1, min = 0, max = 1),
                                 numericInput("stability_q_step", "Step:", value = 0.1, step = 0.05, min = 0.01, max = 0.5)
                          )
                        )
                    ),

                    # Cross-lagged parameters
                    div(style = "margin-bottom: 15px;",
                        h5("Cross-Lagged Parameters"),
                        fluidRow(
                          column(6,
                                 numericInput("cross_p_min", "Cross P (min):", value = 0.0, step = 0.1, min = -1, max = 1),
                                 numericInput("cross_p_max", "Cross P (max):", value = 0.1, step = 0.1, min = -1, max = 1),
                                 numericInput("cross_p_step", "Step:", value = 0.1, step = 0.05, min = 0.01, max = 0.5)
                          ),
                          column(6,
                                 numericInput("cross_q_min", "Cross Q (min):", value = 0.0, step = 0.1, min = -1, max = 1),
                                 numericInput("cross_q_max", "Cross Q (max):", value = 0.1, step = 0.1, min = -1, max = 1),
                                 numericInput("cross_q_step", "Step:", value = 0.1, step = 0.05, min = 0.01, max = 0.5)
                          )
                        )
                    ),

                    # Variance parameters
                    div(style = "margin-bottom: 15px;",
                        h5("Variance Parameters"),
                        fluidRow(
                          column(6,
                                 numericInput("variance_p", "Variance P:", value = 0.5, step = 0.1, min = 0.01, max = 5),
                                 numericInput("variance_q", "Variance Q:", value = 0.5, step = 0.1, min = 0.01, max = 5)
                          ),
                          column(6,
                                 numericInput("variance_between_x", "Variance Between X:", value = 0.8, step = 0.1, min = 0.01, max = 5),
                                 numericInput("variance_between_y", "Variance Between Y:", value = 0.8, step = 0.1, min = 0.01, max = 5)
                          )
                        )
                    ),

                    numericInput("cov_pq", "Covariance PQ:", value = 0, step = 0.1, min = -2, max = 2)
                )
              ),

              # Confounder Parameters (conditional)
              conditionalPanel(
                condition = "input.data_generation == 'clpm_confounder'",
                fluidRow(
                  box(title = "Confounder Parameters", status = "warning", solidHeader = TRUE, width = 12,
                      fluidRow(
                        column(3,
                               numericInput("confounder_p", "Confounder Effect on P:", value = 0.3, step = 0.1, min = -2, max = 2)
                        ),
                        column(3,
                               numericInput("confounder_q", "Confounder Effect on Q:", value = 0.3, step = 0.1, min = -2, max = 2)
                        ),
                        column(3,
                               numericInput("confounder_variance", "Confounder Variance:", value = 1, step = 0.1, min = 0.01, max = 5)
                        ),
                        column(3,
                               numericInput("confounder_stability", "Confounder Stability:", value = 0.4, step = 0.1, min = 0, max = 1)
                        )
                      )
                  )
                )
              ),

              # Parameter Preview and Run Button
              fluidRow(
                box(title = "Parameter Grid Preview", status = "success", solidHeader = TRUE, width = 8,
                    DT::dataTableOutput("param_preview"),
                    br(),
                    textOutput("grid_size_info")
                ),

                box(title = "Run Simulation", status = "success", solidHeader = TRUE, width = 4,
                    br(),
                    actionButton("run_simulation", "Run Monte Carlo Simulation",
                                 class = "btn-primary btn-lg",
                                 style = "width: 100%; margin-bottom: 20px;"),
                    br(),
                    div(id = "simulation_status",
                        style = "text-align: center; margin-top: 10px;"),
                    br(),
                    downloadButton("download_results", "Download Results",
                                   class = "btn-success",
                                   style = "width: 100%;")
                )
              )
      ),

      # Results Tab
      tabItem(tabName = "results",
              fluidRow(
                box(title = "Simulation Results", status = "primary", solidHeader = TRUE, width = 12,
                    DT::dataTableOutput("results_table")
                )
              ),

              fluidRow(
                box(title = "Results Visualization", status = "info", solidHeader = TRUE, width = 12,
                    selectInput("plot_var_x", "X-axis Variable:", choices = NULL),
                    selectInput("plot_var_y", "Y-axis Variable:", choices = NULL),
                    selectInput("plot_var_color", "Color Variable:", choices = NULL),
                    plotlyOutput("results_plot", height = "500px")
                )
              )
      ),

      # Documentation Tab
      tabItem(tabName = "docs",
              fluidRow(
                box(title = "Parameter Descriptions", status = "info", solidHeader = TRUE, width = 12,
                    h3("Estimator Methods"),
                    tags$ul(
                      tags$li(strong("OLS:"), "Ordinary Least Squares regression"),
                      tags$li(strong("RICLPM:"), "Random Intercept Cross-Lagged Panel Model"),
                      tags$li(strong("CLPM:"), "Cross-Lagged Panel Model"),
                      tags$li(strong("CTSEM:"), "Continuous Time Structural Equation Model")
                    ),

                    h3("Data Generation Processes"),
                    tags$ul(
                      tags$li(strong("CLPM:"), "Standard Cross-Lagged Panel Model"),
                      tags$li(strong("RI-CLPM:"), "Random Intercept Cross-Lagged Panel Model"),
                      tags$li(strong("CLPM with Confounder:"), "CLPM with unmeasured confounder affecting both variables")
                    ),

                    h3("Parameters"),
                    tags$ul(
                      tags$li(strong("Stability P/Q:"), "Autoregressive effects (stability of variables over time)"),
                      tags$li(strong("Cross P/Q:"), "Cross-lagged effects (causal effects between variables)"),
                      tags$li(strong("Variance P/Q:"), "Within-person variance of latent variables"),
                      tags$li(strong("Variance Between X/Y:"), "Between-person variance"),
                      tags$li(strong("Covariance PQ:"), "Covariance between variables at same time point"),
                      tags$li(strong("Confounder Parameters:"), "Effects and properties of unmeasured confounder")
                    ),

                    h3("Output Interpretation"),
                    p("The simulation results show how well each estimator recovers the true parameter values under different conditions. Look for bias (difference from true values) and efficiency (standard errors) in the estimates.")
                )
              )
      )
    )
  )
)
