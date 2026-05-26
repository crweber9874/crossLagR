library(shiny)
library(lavaan)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load crossLagR functions
devtools::load_all(here::here())
source("helpers.R", local = TRUE)

# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  withMathJax(),
  tags$head(tags$style(HTML("
    .param-section { background: #f8f9fa; padding: 12px; border-radius: 8px;
                     margin-bottom: 10px; border: 1px solid #dee2e6; }
    .param-section h5 { margin-top: 0; color: #495057; }
    .caveat-banner { background: #fff3cd; border: 1px solid #ffc107;
                     padding: 10px; border-radius: 6px; margin-bottom: 15px; }
    .btn-run { margin-top: 10px; }
    .results-panel { min-height: 400px; }
  "))),

  titlePanel(
    div(
      h2("crossLagR: Simulation Laboratory"),
      p(style = "font-size: 14px; color: grey;",
        "Simulate data from a DGP, estimate with any model, compare true vs. estimated values, and run Monte Carlo studies.")
    )
  ),

  tabsetPanel(
    id = "main_tabs",

    # ====================================================================
    # TAB 1: Single Sample
    # ====================================================================
    tabPanel("Single Sample",
      sidebarLayout(
        sidebarPanel(
          width = 3,

          # --- DGP Selection ---
          div(class = "param-section",
            h5("Data-Generating Process"),
            selectInput("dgp", "DGP:",
              choices = names(dgp_registry), selected = "RI-CLPM"),
            textOutput("dgp_desc")
          ),

          # --- Common Parameters ---
          div(class = "param-section",
            h5("Design"),
            fluidRow(
              column(4, numericInput("waves", "Waves", 5, min = 3, max = 20)),
              column(4, numericInput("sample_size", "N", 1000, min = 50, max = 20000, step = 100)),
              column(4, numericInput("seed", "Seed", 42, min = 1))
            )
          ),

          # --- DGP-Specific Parameters ---
          uiOutput("dgp_params_ui"),

          # --- Estimator ---
          div(class = "param-section",
            h5("Estimator"),
            selectInput("estimator", "Statistical Model:",
              choices = names(estimator_registry), selected = "CLPM"),
            checkboxInput("constrain_est", "Constrain coefficients across waves", FALSE)
          ),

          hr(),
          actionButton("run_single", "Simulate & Estimate",
                        class = "btn-primary btn-lg btn-block btn-run"),
          br(), br(),
          textOutput("single_status")
        ),

        mainPanel(
          width = 9,
          div(class = "results-panel",
            uiOutput("caveat_banner"),
            h4("True vs. Estimated Parameters"),
            tableOutput("single_comparison"),
            hr(),
            h4("Fit Indices"),
            tableOutput("single_fit"),
            hr(),
            h4("Generated lavaan Syntax"),
            fluidRow(
              column(6, h5("DGP Syntax"), verbatimTextOutput("dgp_syntax")),
              column(6, h5("Estimation Syntax"), verbatimTextOutput("est_syntax"))
            )
          )
        )
      )
    ),

    # ====================================================================
    # TAB 2: Monte Carlo
    # ====================================================================
    tabPanel("Monte Carlo",
      sidebarLayout(
        sidebarPanel(
          width = 3,

          div(class = "param-section",
            h5("Monte Carlo Settings"),
            helpText("Uses the DGP and parameters from the Single Sample tab."),
            sliderInput("mc_trials", "Number of Trials (R)",
                        min = 10, max = 1000, value = 100, step = 10),
            checkboxGroupInput("mc_estimators", "Estimators to Compare:",
              choices = names(estimator_registry),
              selected = c("CLPM", "RI-CLPM")),
            checkboxGroupInput("mc_waves_set", "Wave Counts to Compare:",
              choices = c("3" = 3, "5" = 5, "7" = 7, "10" = 10),
              selected = c("5"), inline = TRUE),
            checkboxGroupInput("mc_n_set", "Sample Sizes to Compare:",
              choices = c("200" = 200, "500" = 500, "1000" = 1000,
                          "2000" = 2000, "5000" = 5000),
              selected = c("1000"), inline = TRUE)
          ),

          div(class = "param-section",
            h5("Plot Options"),
            numericInput("mc_bins", "Histogram bins", 30, min = 10, max = 80),
            selectInput("mc_theme", "Theme:",
              choices = c("minimal", "bw", "classic"), selected = "minimal")
          ),

          hr(),
          actionButton("run_mc", "Run Monte Carlo",
                        class = "btn-warning btn-lg btn-block btn-run"),
          br(), br(),
          textOutput("mc_status"),
          hr(),
          downloadButton("mc_download_csv", "Download Results (CSV)"),
          downloadButton("mc_download_rdata", "Download Results (.RData)")
        ),

        mainPanel(
          width = 9,
          tabsetPanel(
            id = "mc_subtabs",
            tabPanel("Sampling Distributions",
              plotOutput("mc_dist_plot", height = "550px")
            ),
            tabPanel("Bias Summary",
              plotOutput("mc_bias_plot", height = "400px"),
              hr(),
              tableOutput("mc_summary_table")
            ),
            tabPanel("Bias by Condition",
              plotOutput("mc_bias_condition_plot", height = "450px")
            ),
            tabPanel("Convergence",
              tableOutput("mc_convergence_table"),
              plotOutput("mc_convergence_plot", height = "300px")
            )
          )
        )
      )
    )
  )
)


# ============================================================================
# Server
# ============================================================================

server <- function(input, output, session) {

  # ---- DGP description ----
  output$dgp_desc <- renderText({
    dgp_registry[[input$dgp]]$desc
  })

  # ---- Dynamic DGP parameter UI ----
  output$dgp_params_ui <- renderUI({
    dgp <- input$dgp
    defaults <- dgp_registry[[dgp]]$params

    if (dgp == "CLPM") {
      div(class = "param-section",
        h5("Structural Parameters"),
        fluidRow(
          column(6, numericInput("p_beta_x", HTML("\\(\\beta_{x}\\) AR(X)"),
                                 defaults$beta_x, min = -0.9, max = 0.9, step = 0.05)),
          column(6, numericInput("p_beta_y", HTML("\\(\\beta_{y}\\) AR(Y)"),
                                 defaults$beta_y, min = -0.9, max = 0.9, step = 0.05))
        ),
        fluidRow(
          column(6, numericInput("p_omega_xy", HTML("\\(\\omega_{xy}\\) CL(X\\(\\to\\)Y)"),
                                 defaults$omega_xy, min = -0.5, max = 0.5, step = 0.05)),
          column(6, numericInput("p_omega_yx", HTML("\\(\\omega_{yx}\\) CL(Y\\(\\to\\)X)"),
                                 defaults$omega_yx, min = -0.5, max = 0.5, step = 0.05))
        ),
        h5("Variance Parameters"),
        fluidRow(
          column(4, numericInput("p_var_x", HTML("\\(\\sigma^2_x\\)"), defaults$var_x,
                                  min = 0.1, max = 5, step = 0.1)),
          column(4, numericInput("p_var_y", HTML("\\(\\sigma^2_y\\)"), defaults$var_y,
                                  min = 0.1, max = 5, step = 0.1)),
          column(4, numericInput("p_cov_xy", HTML("Cov(x,y)"), defaults$cov_xy,
                                  min = -2, max = 2, step = 0.1))
        )
      )
    } else if (dgp == "RI-CLPM") {
      div(class = "param-section",
        h5("Within-Person (State) Parameters"),
        fluidRow(
          column(6, numericInput("p_beta_x", "AR(X) beta_x",
                                 defaults$beta_x, min = -0.9, max = 0.9, step = 0.05)),
          column(6, numericInput("p_beta_y", "AR(Y) beta_y",
                                 defaults$beta_y, min = -0.9, max = 0.9, step = 0.05))
        ),
        fluidRow(
          column(6, numericInput("p_omega_xy", "CL(X->Y) omega_xy",
                                 defaults$omega_xy, min = -0.5, max = 0.5, step = 0.05)),
          column(6, numericInput("p_omega_yx", "CL(Y->X) omega_yx",
                                 defaults$omega_yx, min = -0.5, max = 0.5, step = 0.05))
        ),
        numericInput("p_cov_pq", "Within-person Cov(p,q)", defaults$cov_pq,
                     min = -1, max = 1, step = 0.05),
        hr(),
        h5("Variance Decomposition (ICC)"),
        helpText("Set total variance and the proportion that is between-person (trait). ",
                 "Higher ICC = more trait variance = more bias when ignoring traits."),
        fluidRow(
          column(6, numericInput("p_total_var_x", "Total Var(X)", 2.5,
                                 min = 0.5, max = 10, step = 0.5)),
          column(6, numericInput("p_total_var_y", "Total Var(Y)", 2.5,
                                 min = 0.5, max = 10, step = 0.5))
        ),
        fluidRow(
          column(6, sliderInput("p_icc_x", "ICC(X)",
                                min = 0, max = 0.95, value = 0.60, step = 0.05)),
          column(6, sliderInput("p_icc_y", "ICC(Y)",
                                min = 0, max = 0.95, value = 0.60, step = 0.05))
        ),
        numericInput("p_cov_BXBY", "Trait Cov(eta_x, eta_y)",
                     defaults$cov_BXBY, min = -3, max = 5, step = 0.1),
        div(style = "background: #e9ecef; padding: 8px; border-radius: 4px; margin-top: 8px;",
          textOutput("icc_display")
        )
      )
    } else if (dgp == "Latent Change") {
      div(class = "param-section",
        h5("Proportional Effects"),
        fluidRow(
          column(6, numericInput("p_beta_x", "beta_x (X level -> X change)",
                                 defaults$beta_x, min = -0.9, max = 0.9, step = 0.05)),
          column(6, numericInput("p_beta_y", "beta_y (Y level -> Y change)",
                                 defaults$beta_y, min = -0.9, max = 0.9, step = 0.05))
        ),
        h5("Coupling Effects"),
        fluidRow(
          column(6, numericInput("p_omega_x", "gamma_x (Y level -> X change)",
                                 defaults$omega_x, min = -0.5, max = 0.5, step = 0.05)),
          column(6, numericInput("p_omega_y", "gamma_y (X level -> Y change)",
                                 defaults$omega_y, min = -0.5, max = 0.5, step = 0.05))
        ),
        h5("Variance"),
        fluidRow(
          column(4, numericInput("p_initial_var_x", "Var(X1)", defaults$initial_var_x,
                                  min = 0.1, max = 5, step = 0.1)),
          column(4, numericInput("p_initial_var_y", "Var(Y1)", defaults$initial_var_y,
                                  min = 0.1, max = 5, step = 0.1)),
          column(4, numericInput("p_indicator_var", "Indicator Var",
                                  defaults$indicator_variance_x, min = 0, max = 3, step = 0.1))
        )
      )
    }
  })

  # ---- ICC display for RI-CLPM ----
  output$icc_display <- renderText({
    req(input$dgp == "RI-CLPM")
    var_BX <- input$p_total_var_x * input$p_icc_x
    var_p  <- input$p_total_var_x * (1 - input$p_icc_x)
    var_BY <- input$p_total_var_y * input$p_icc_y
    var_q  <- input$p_total_var_y * (1 - input$p_icc_y)
    paste0("Var(trait_x) = ", round(var_BX, 2), ", Var(state_x) = ", round(var_p, 2),
           "  |  Var(trait_y) = ", round(var_BY, 2), ", Var(state_y) = ", round(var_q, 2))
  })

  # ---- Collect DGP params into a list ----
  dgp_params <- reactive({
    dgp <- input$dgp
    if (dgp == "CLPM") {
      list(beta_x = input$p_beta_x, beta_y = input$p_beta_y,
           omega_xy = input$p_omega_xy, omega_yx = input$p_omega_yx,
           var_x = input$p_var_x, var_y = input$p_var_y, cov_xy = input$p_cov_xy)
    } else if (dgp == "RI-CLPM") {
      var_BX <- input$p_total_var_x * input$p_icc_x
      var_p  <- input$p_total_var_x * (1 - input$p_icc_x)
      var_BY <- input$p_total_var_y * input$p_icc_y
      var_q  <- input$p_total_var_y * (1 - input$p_icc_y)
      list(beta_x = input$p_beta_x, beta_y = input$p_beta_y,
           omega_xy = input$p_omega_xy, omega_yx = input$p_omega_yx,
           var_p = var_p, var_q = var_q, cov_pq = input$p_cov_pq,
           var_BX = var_BX, var_BY = var_BY, cov_BXBY = input$p_cov_BXBY)
    } else if (dgp == "Latent Change") {
      list(beta_x = input$p_beta_x, beta_y = input$p_beta_y,
           omega_x = input$p_omega_x, omega_y = input$p_omega_y,
           initial_var_x = input$p_initial_var_x, initial_var_y = input$p_initial_var_y,
           indicator_variance_x = input$p_indicator_var,
           indicator_variance_y = input$p_indicator_var)
    }
  })

  # ---- Caveat banner ----
  output$caveat_banner <- renderUI({
    dgp <- input$dgp
    est <- input$estimator
    txt <- get_caveat_text(dgp, est)
    if (!is.null(txt)) {
      div(class = "caveat-banner",
        tags$strong(paste0("\u26A0 DGP: ", dgp, " | Estimator: ", est)),
        br(), br(),
        txt
      )
    }
  })

  # ==== Single Sample ====
  single_result <- reactiveVal(NULL)

  observeEvent(input$run_single, {
    req(dgp_params())
    set.seed(input$seed)

    withProgress(message = "Simulating data...", value = 0.3, {
      dat <- simulate_data(input$dgp, input$waves, input$sample_size, dgp_params())
      incProgress(0.3, detail = "Fitting model...")
      result <- fit_estimator(input$estimator, dat, input$waves)
      incProgress(0.4, detail = "Done!")
      single_result(result)
    })
  })

  output$single_status <- renderText({
    res <- single_result()
    if (is.null(res)) return("Press 'Simulate & Estimate' to start.")
    if (res$converged) "✅ Model converged." else "❌ Model did not converge."
  })

  output$single_comparison <- renderTable({
    res <- single_result(); req(res, res$converged)
    true <- get_true_values(input$dgp, dgp_params())
    matched <- match_true_values(res$params, true)
    matched |>
      dplyr::select(Parameter = label, True = true,
                     Estimate = est.std, SE = se, Bias = bias, p = pvalue)
  }, striped = TRUE, bordered = TRUE)

  output$single_fit <- renderTable({
    res <- single_result(); req(res, res$converged)
    data.frame(
      CFI = round(res$fit_measures["cfi"], 3),
      RMSEA = round(res$fit_measures["rmsea"], 3),
      SRMR = round(res$fit_measures["srmr"], 3)
    )
  }, striped = TRUE)

  output$dgp_syntax <- renderPrint({
    req(dgp_params())
    dgp <- input$dgp
    reg <- dgp_registry[[dgp]]
    fn <- match.fun(reg$sim_fn)
    args <- c(list(waves = input$waves), dgp_params())
    if (dgp == "CLPM") args$sample_size <- 10
    else if (dgp == "Latent Change") {
      args$sample.nobs <- 10; args$variable_type <- "bivariate"
    } else args$sample.nobs <- 10
    result <- tryCatch(do.call(fn, args), error = function(e) NULL)
    if (!is.null(result)) cat(result$model) else cat("(error generating syntax)")
  })

  output$est_syntax <- renderPrint({
    reg <- estimator_registry[[input$estimator]]
    fn <- match.fun(reg$est_fn)
    args <- c(list(waves = input$waves), reg$args)
    syntax <- tryCatch(do.call(fn, args), error = function(e) "(error)")
    cat(syntax)
  })

  # ==== Monte Carlo ====
  mc_results <- reactiveVal(NULL)

  observeEvent(input$run_mc, {
    req(dgp_params(), length(input$mc_estimators) > 0)

    estimators <- input$mc_estimators
    n_trials <- input$mc_trials
    waves_set <- as.integer(input$mc_waves_set)
    n_set <- as.integer(input$mc_n_set)
    params <- dgp_params()
    dgp <- input$dgp

    total <- n_trials * length(estimators) * length(waves_set) * length(n_set)
    completed <- 0

    withProgress(message = "Running Monte Carlo...", value = 0, {
      all_results <- list()

      for (w in waves_set) {
        for (n in n_set) {
          for (i in seq_len(n_trials)) {
            dat <- tryCatch(
              simulate_data(dgp, w, n, params),
              error = function(e) NULL
            )
            if (is.null(dat)) next

            for (est in estimators) {
              res <- fit_estimator(est, dat, w)
              if (res$converged && nrow(res$params) > 0) {
                row <- res$params
                row$estimator <- est
                row$trial <- i
                row$waves <- w
                row$n <- n
                row$converged <- TRUE
                all_results[[length(all_results) + 1]] <- row
              } else {
                all_results[[length(all_results) + 1]] <- data.frame(
                  label = NA, est.std = NA, se = NA, pvalue = NA,
                  ci.lower = NA, ci.upper = NA,
                  estimator = est, trial = i, waves = w, n = n,
                  converged = FALSE
                )
              }

              completed <- completed + 1
              incProgress(1 / total,
                detail = paste0("W=", w, " N=", n, " ", est, " trial ", i))
            }
          }
        }
      }

      mc_results(dplyr::bind_rows(all_results))
    })
  })

  output$mc_status <- renderText({
    res <- mc_results()
    if (is.null(res)) return("Configure and press 'Run Monte Carlo'.")
    n_conv <- sum(res$converged, na.rm = TRUE)
    paste0("Done: ", nrow(res), " runs (", n_conv, " converged)")
  })

  # ---- MC Sampling Distribution Plot ----
  output$mc_dist_plot <- renderPlot({
    res <- mc_results(); req(res)
    res_conv <- res |> filter(converged, !is.na(est.std)) |>
      mutate(base_label = strip_wave(label))
    req(nrow(res_conv) > 0)

    true <- get_true_values(input$dgp, dgp_params())
    theme_fn <- switch(input$mc_theme,
      minimal = theme_minimal, bw = theme_bw, classic = theme_classic,
      theme_minimal)

    ggplot(res_conv, aes(x = est.std, fill = estimator)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = input$mc_bins, alpha = 0.5, position = "identity") +
      geom_density(aes(color = estimator), linewidth = 0.8, fill = NA) +
      geom_vline(data = true, aes(xintercept = true),
                 linetype = "dashed", color = "red", linewidth = 1) +
      facet_wrap(~base_label, scales = "free") +
      labs(title = paste("Sampling Distributions | DGP:", input$dgp),
           subtitle = paste0("R = ", input$mc_trials, " | N = ",
                             paste(input$mc_n_set, collapse = ","),
                             " | Waves = ", paste(input$mc_waves_set, collapse = ",")),
           x = "Estimate", y = "Density", fill = "Estimator", color = "Estimator") +
      theme_fn(base_size = 14) +
      theme(legend.position = "bottom")
  })

  # ---- MC Bias Plot ----
  output$mc_bias_plot <- renderPlot({
    res <- mc_results(); req(res)
    res_conv <- res |> filter(converged, !is.na(est.std))
    req(nrow(res_conv) > 0)

    true <- get_true_values(input$dgp, dgp_params())
    summ <- mc_summary(res_conv, true)

    ggplot(summ, aes(x = base_label, y = bias, fill = estimator)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_hline(yintercept = 0, linewidth = 0.8) +
      labs(title = "Bias by Parameter and Estimator",
           x = "Parameter", y = "Bias (Mean Est. - True)", fill = "Estimator") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 30, hjust = 1))
  })

  # ---- MC Summary Table ----
  output$mc_summary_table <- renderTable({
    res <- mc_results(); req(res)
    res_conv <- res |> filter(converged, !is.na(est.std))
    req(nrow(res_conv) > 0)

    true <- get_true_values(input$dgp, dgp_params())
    summ <- mc_summary(res_conv, true)
    summ |>
      mutate(across(c(mean_est, sd_est, bias, rmse), ~ round(.x, 4)),
             true = round(true, 4)) |>
      select(Estimator = estimator, Parameter = base_label, True = true,
             `Mean Est.` = mean_est, SD = sd_est, Bias = bias, RMSE = rmse)
  }, striped = TRUE, bordered = TRUE)

  # ---- MC Bias by Condition ----
  output$mc_bias_condition_plot <- renderPlot({
    res <- mc_results(); req(res)
    res_conv <- res |> filter(converged, !is.na(est.std))
    req(nrow(res_conv) > 0)

    true <- get_true_values(input$dgp, dgp_params())

    cond_summ <- res_conv |>
      mutate(base_label = strip_wave(label)) |>
      group_by(estimator, base_label, waves, n) |>
      summarise(mean_est = mean(est.std, na.rm = TRUE), .groups = "drop") |>
      left_join(true, by = "base_label") |>
      mutate(bias = mean_est - true,
             condition = paste0("W=", waves, " N=", n))

    ggplot(cond_summ, aes(x = condition, y = bias, fill = estimator)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_hline(yintercept = 0) +
      facet_wrap(~base_label, scales = "free_y") +
      labs(title = "Bias by Condition", x = "", y = "Bias", fill = "Estimator") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # ---- MC Convergence ----
  output$mc_convergence_table <- renderTable({
    res <- mc_results(); req(res)
    res |>
      group_by(estimator, waves, n) |>
      summarise(
        Total = n(),
        Converged = sum(converged, na.rm = TRUE),
        Rate = round(mean(converged, na.rm = TRUE), 3),
        .groups = "drop"
      ) |>
      rename(Estimator = estimator, Waves = waves, N = n)
  }, striped = TRUE, bordered = TRUE)

  output$mc_convergence_plot <- renderPlot({
    res <- mc_results(); req(res)
    conv <- res |>
      group_by(estimator, waves, n) |>
      summarise(rate = mean(converged, na.rm = TRUE), .groups = "drop") |>
      mutate(condition = paste0("W=", waves, " N=", n))

    ggplot(conv, aes(x = condition, y = rate, fill = estimator)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      ylim(0, 1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      labs(title = "Convergence Rate", x = "", y = "Rate", fill = "Estimator") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })

  # ---- Downloads ----
  output$mc_download_csv <- downloadHandler(
    filename = function() paste0("mc_results_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) {
      res <- mc_results(); req(res)
      write.csv(res, file, row.names = FALSE)
    }
  )

  output$mc_download_rdata <- downloadHandler(
    filename = function() paste0("mc_results_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData"),
    content = function(file) {
      mc_sim_results <- mc_results(); req(mc_sim_results)
      save(mc_sim_results, file = file)
    }
  )
}

shinyApp(ui = ui, server = server)
