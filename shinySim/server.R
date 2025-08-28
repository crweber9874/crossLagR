server <- function(input, output, session) {

  # Reactive values to store simulation results
  values <- reactiveValues(
    results = NULL,
    param_grid = NULL,
    simulation_running = FALSE
  )

  # Generate parameter grid based on inputs
  param_grid <- reactive({
    # Base parameters
    base_grid <- expand.grid(
      stability_p = seq(input$stability_p_min, input$stability_p_max, by = input$stability_p_step),
      stability_q = seq(input$stability_q_min, input$stability_q_max, by = input$stability_q_step),
      cross_p = seq(input$cross_p_min, input$cross_p_max, by = input$cross_p_step),
      cross_q = seq(input$cross_q_min, input$cross_q_max, by = input$cross_q_step),
      variance_p = input$variance_p,
      variance_q = input$variance_q,
      variance_between_x = input$variance_between_x,
      variance_between_y = input$variance_between_y,
      cov_pq = input$cov_pq
    )

    # Add confounder parameters if needed
    if (input$data_generation == "clpm_confounder") {
      base_grid$confounder_p <- input$confounder_p
      base_grid$confounder_q <- input$confounder_q
      base_grid$confounder_variance <- input$confounder_variance
      base_grid$confounder_stability <- input$confounder_stability
    }

    values$param_grid <- base_grid
    return(base_grid)
  })

  # Parameter grid preview
  output$param_preview <- DT::renderDataTable({
    grid <- param_grid()
    DT::datatable(grid,
                  options = list(scrollX = TRUE, pageLength = 5),
                  caption = "Parameter combinations to be tested") %>%
      DT::formatRound(columns = 1:ncol(grid), digits = 3)
  })

  # Grid size information
  output$grid_size_info <- renderText({
    grid <- param_grid()
    total_sims <- nrow(grid) * input$trials
    paste0("Total parameter combinations: ", nrow(grid),
           " | Total simulations: ", format(total_sims, big.mark = ","))
  })

  # Run simulation
  observeEvent(input$run_simulation, {
    values$simulation_running <- TRUE

    # Update UI to show running status
    output$simulation_status <- renderText({
      "Simulation running... Please wait."
    })

    # Disable the run button
    shinyjs::disable("run_simulation")

    # Run the simulation
    tryCatch({
      results <- run_mc_sims(
        estimator = input$estimator,
        param_grid = param_grid(),
        trials = input$trials,
        waves = input$waves,
        sample_size = input$sample_size,
        verbose = input$verbose,
        data_generation = input$data_generation
      )

      values$results <- results
      values$simulation_running <- FALSE

      # Update status
      output$simulation_status <- renderText({
        paste0("Simulation completed successfully! ", nrow(results), " results generated.")
      })

      # Update plot variable choices
      numeric_cols <- names(results)[sapply(results, is.numeric)]
      updateSelectInput(session, "plot_var_x", choices = numeric_cols, selected = numeric_cols[1])
      updateSelectInput(session, "plot_var_y", choices = numeric_cols, selected = numeric_cols[2])

      all_cols <- names(results)
      updateSelectInput(session, "plot_var_color",
                        choices = c("None" = "", all_cols),
                        selected = "")

    }, error = function(e) {
      values$simulation_running <- FALSE
      output$simulation_status <- renderText({
        paste("Error:", e$message)
      })
    })

    # Re-enable the run button
    shinyjs::enable("run_simulation")
  })

  # Results table
  output$results_table <- DT::renderDataTable({
    req(values$results)
    DT::datatable(values$results,
                  options = list(scrollX = TRUE, pageLength = 10)) %>%
      DT::formatRound(columns = which(sapply(values$results, is.numeric)), digits = 4)
  })

  # Results plot
  output$results_plot <- renderPlotly({
    req(values$results, input$plot_var_x, input$plot_var_y)

    p <- ggplot(values$results, aes_string(x = input$plot_var_x, y = input$plot_var_y))

    if (input$plot_var_color != "") {
      p <- p + aes_string(color = input$plot_var_color)
    }

    p <- p +
      geom_point(alpha = 0.7, size = 2) +
      theme_minimal() +
      labs(title = paste("Monte Carlo Results:", input$estimator, "Estimator"),
           subtitle = paste("DGP:", input$data_generation))

    ggplotly(p)
  })

  # Download handler
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("monte_carlo_results_", input$estimator, "_",
             input$data_generation, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$results)
      write.csv(values$results, file, row.names = FALSE)
    }
  )
}
