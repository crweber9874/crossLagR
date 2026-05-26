server <- function(input, output, session) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  values <- reactiveValues(
    results = NULL,
    param_grid = NULL,
    sim_status = "Initialize models before running simulation.",
    models_initialized = FALSE,
    init_config = NULL,
    init_syntax = NULL
  )

  current_model_config <- reactive({
    list(
      estimator = input$estimator,
      waves = input$waves,
      constrain_beta = isTRUE(input$constrain_beta),
      constrain_omega = isTRUE(input$constrain_omega),
      constrain_resvar = isTRUE(input$constrain_resvar),
      constrain_rescov = isTRUE(input$constrain_rescov),
      estimate_means = isTRUE(input$estimate_means),
      start_values = isTRUE(input$start_values),
      riclpm_type = input$riclpm_type %||% "riclpm",
      lchange_type = input$lchange_type %||% "dual_change"
    )
  })

  model_is_ready <- reactive({
    isTRUE(values$models_initialized) && identical(values$init_config, current_model_config())
  })

  observe({
    if (model_is_ready()) {
      shinyjs::enable("run_simulation")
    } else {
      shinyjs::disable("run_simulation")
    }
  })

  ## ---- Parameter grid ----
  param_grid <- reactive({
    expand.grid(
      stability_p        = seq(input$stability_p_min, input$stability_p_max, by = input$stability_p_step),
      stability_q        = seq(input$stability_q_min, input$stability_q_max, by = input$stability_q_step),
      cross_p            = seq(input$cross_p_min,     input$cross_p_max,     by = input$cross_p_step),
      cross_q            = seq(input$cross_q_min,     input$cross_q_max,     by = input$cross_q_step),
      variance_p         = input$variance_p,
      variance_q         = input$variance_q,
      variance_between_x = input$variance_between_x,
      variance_between_y = input$variance_between_y,
      cov_pq             = input$cov_pq
    )
  })

  output$param_preview <- DT::renderDataTable({
    grid <- param_grid()
    DT::datatable(grid,
                  options = list(scrollX = TRUE, pageLength = 5),
                  caption = "Parameter cells to be tested") |>
      DT::formatRound(columns = seq_len(ncol(grid)), digits = 3)
  })

  output$grid_size_info <- renderText({
    grid <- param_grid()
    paste0("Cells: ", nrow(grid),
           " | Total fits: ", format(nrow(grid) * input$trials, big.mark = ","))
  })

  ## ---- Code preview (lavaan model syntax + reproducible R call) ----
  output$model_syntax <- renderText({
    if (model_is_ready() && !is.null(values$init_syntax)) {
      return(values$init_syntax)
    }

    args <- list(
      waves                          = input$waves,
      constrain_beta                 = isTRUE(input$constrain_beta),
      constrain_omega                = isTRUE(input$constrain_omega),
      constrain_residual_variances   = isTRUE(input$constrain_resvar),
      constrain_residual_covariances = isTRUE(input$constrain_rescov),
      estimate_means                 = isTRUE(input$estimate_means),
      start_values                   = isTRUE(input$start_values)
    )
    tryCatch({
      switch(input$estimator,
        "CLPM"    = do.call(crossLagR::estimateCLPM,   args),
        "RICLPM"  = do.call(crossLagR::estimateRICLPM,
                            c(args["waves"],
                              list(riclpm_type = input$riclpm_type %||% "riclpm"))),
        "LCHANGE" = crossLagR::estimateLChange(
                      waves = input$waves,
                      lchange_type = input$lchange_type %||% "dual_change"),
        paste0("# Syntax preview not available for estimator: ", input$estimator)
      )
    }, error = function(e) paste("# Error generating syntax:", e$message))
  })

  output$r_code <- renderText({
    extra <- ""
    if (input$estimator == "RICLPM")
      extra <- paste0(",\n  riclpm_type     = \"",
                      input$riclpm_type %||% "riclpm", "\"")
    if (input$estimator == "LCHANGE")
      extra <- paste0(",\n  lchange_type    = \"",
                      input$lchange_type %||% "dual_change", "\"")

    paste0(
      'library(crossLagR)\n\n',
      'results <- run_mc_sims(\n',
      '  estimator       = "', input$estimator, '",', extra, ',\n',
      '  data_generation = "', input$data_generation, '",\n',
      '  trials          = ', input$trials, ',\n',
      '  waves           = ', input$waves, ',\n',
      '  sample_size     = ', input$sample_size, ',\n',
      '  param_grid      = <see Parameter grid preview>\n',
      ')'
    )
  })

  observeEvent(input$initialize_models, {
    args <- list(
      waves                          = input$waves,
      constrain_beta                 = isTRUE(input$constrain_beta),
      constrain_omega                = isTRUE(input$constrain_omega),
      constrain_residual_variances   = isTRUE(input$constrain_resvar),
      constrain_residual_covariances = isTRUE(input$constrain_rescov),
      estimate_means                 = isTRUE(input$estimate_means),
      start_values                   = isTRUE(input$start_values)
    )

    values$models_initialized <- FALSE
    values$init_config <- NULL
    values$init_syntax <- NULL
    values$sim_status <- "Initializing models..."
    output$sim_status_text <- renderText("Initializing models...")

    shinyWidgets::updateProgressBar(
      session = session, id = "sim_progress",
      value = 0, total = 100, title = "Initializing..."
    )

    init_result <- tryCatch({
      syntax <- switch(input$estimator,
        "CLPM" = do.call(crossLagR::estimateCLPM, args),
        "RICLPM" = do.call(crossLagR::estimateRICLPM,
                            c(args["waves"],
                              list(riclpm_type = input$riclpm_type %||% "riclpm"))),
        "LCHANGE" = crossLagR::estimateLChange(
          waves = input$waves,
          lchange_type = input$lchange_type %||% "dual_change"
        ),
        paste0("# No explicit lavaan syntax initialization required for estimator: ", input$estimator)
      )

      list(ok = TRUE, syntax = syntax)
    }, error = function(e) {
      list(ok = FALSE, message = e$message)
    })

    if (!isTRUE(init_result$ok)) {
      values$sim_status <- paste0("Initialization failed: ", init_result$message)
      output$sim_status_text <- renderText(values$sim_status)
      shinyWidgets::updateProgressBar(
        session = session, id = "sim_progress",
        value = 0, total = 100, title = "Initialization failed"
      )
      showNotification(values$sim_status, type = "error")
      return()
    }

    values$models_initialized <- TRUE
    values$init_config <- current_model_config()
    values$init_syntax <- init_result$syntax
    values$sim_status <- "Models initialized. Ready to run simulation."
    output$sim_status_text <- renderText(values$sim_status)

    shinyWidgets::updateProgressBar(
      session = session, id = "sim_progress",
      value = 100, total = 100, title = "Models initialized"
    )
  })

  ## ---- Run simulation, chunked, with progress updates ----
  observeEvent(input$run_simulation, {
    if (!model_is_ready()) {
      values$sim_status <- "Initialize models with current settings before running simulation."
      output$sim_status_text <- renderText(values$sim_status)
      showNotification(values$sim_status, type = "warning")
      return()
    }

    grid  <- param_grid()
    n_row <- nrow(grid)
    if (n_row == 0) return()

    shinyjs::disable("run_simulation")
    values$results <- NULL
    values$sim_status <- "Running…"
    output$sim_status_text <- renderText("Starting simulation…")

    shinyWidgets::updateProgressBar(
      session = session, id = "sim_progress",
      value = 0, total = 100, title = "Starting…"
    )

    start_time   <- Sys.time()
    results_list <- vector("list", n_row)

    for (i in seq_len(n_row)) {
      results_list[[i]] <- tryCatch(
        run_mc_sims(
          estimator       = input$estimator,
          riclpm_type     = input$riclpm_type  %||% "riclpm",
          lchange_type    = input$lchange_type %||% "dual_change",
          param_grid      = grid[i, , drop = FALSE],
          trials          = input$trials,
          waves           = input$waves,
          sample_size     = input$sample_size,
          verbose         = isTRUE(input$verbose),
          data_generation = input$data_generation
        ),
        error = function(e) {
          message("Cell ", i, " failed: ", e$message)
          NULL
        }
      )

      pct     <- round(100 * i / n_row)
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      eta     <- if (i > 0) elapsed / i * (n_row - i) else NA_real_

      shinyWidgets::updateProgressBar(
        session = session, id = "sim_progress",
        value   = pct, total = 100,
        title   = sprintf("Cell %d/%d  •  %.0fs elapsed  •  ~%.0fs remaining",
                          i, n_row, elapsed, eta)
      )
    }

    out <- dplyr::bind_rows(results_list)
    values$results <- out
    values$sim_status <- "Done"

    output$sim_status_text <- renderText({
      sprintf("Done. %s rows in %.1fs.",
              format(nrow(out), big.mark = ","),
              as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    })

    if (nrow(out) > 0) {
      numeric_cols <- names(out)[vapply(out, is.numeric, logical(1))]
      updateSelectInput(session, "plot_var_x", choices = numeric_cols,
                        selected = numeric_cols[1])
      updateSelectInput(session, "plot_var_y", choices = numeric_cols,
                        selected = if (length(numeric_cols) >= 2) numeric_cols[2]
                                   else numeric_cols[1])
      updateSelectInput(session, "plot_var_color",
                        choices = c("None" = "", names(out)),
                        selected = "")
    }

    shinyjs::enable("run_simulation")
  })

  ## ---- Results ----
  output$results_table <- DT::renderDataTable({
    req(values$results)
    DT::datatable(values$results,
                  options = list(scrollX = TRUE, pageLength = 10)) |>
      DT::formatRound(columns = which(vapply(values$results, is.numeric, logical(1))),
                      digits = 4)
  })

  output$results_plot <- renderPlotly({
    req(values$results, input$plot_var_x, input$plot_var_y)
    p <- ggplot(values$results, aes_string(x = input$plot_var_x,
                                           y = input$plot_var_y)) +
      geom_point(alpha = 0.7, size = 2) +
      theme_minimal() +
      labs(title = paste("Monte Carlo:", input$estimator),
           subtitle = paste("DGP:", input$data_generation))
    if (!is.null(input$plot_var_color) && nzchar(input$plot_var_color))
      p <- p + aes_string(color = input$plot_var_color)
    ggplotly(p)
  })

  output$download_results <- downloadHandler(
    filename = function() {
      paste0("mc_", input$estimator, "_", input$data_generation, "_",
             Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$results)
      write.csv(values$results, file, row.names = FALSE)
    }
  )
}
