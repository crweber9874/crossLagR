server <- function(input, output, session) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  ## ---------------------------------------------------------------------------
  ## State
  ## ---------------------------------------------------------------------------
  values <- reactiveValues(
    results        = NULL,
    summary_df     = NULL,
    single_fit     = NULL,
    sim_status     = "Pick an estimator and DGP, then click Run Monte Carlo.",
    init_syntax    = NULL
  )

  mc_bg <- reactiveValues(
    handle        = NULL, log = "", running = FALSE, start = NULL,
    progress_file = NULL, log_file = NULL,
    total_cells   = 0L, done_cells = 0L
  )

  ## Estimators with a static lavaan-syntax preview
  syntax_estimators <- c("CLPM", "RICLPM", "ALT", "LGM", "LCMSR",
                         "LCHANGE", "BB", "TSO")

  ## Shared ggplot theme
  acad_theme <- function() {
    theme_minimal(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ECEFF1", linewidth = 0.3),
        axis.line = element_line(color = "#37474F", linewidth = 0.3),
        plot.title = element_text(color = "#263238", face = "bold", size = 14),
        plot.subtitle = element_text(color = "#546E7A", size = 11),
        legend.position = "bottom",
        strip.background = element_rect(fill = "#ECEFF1", color = NA),
        strip.text = element_text(color = "#263238", face = "bold")
      )
  }
  acad_palette <- c("#37474F", "#C62828", "#2E7D32", "#B7791F",
                    "#1565C0", "#6A1B9A", "#00838F", "#7F7F7F")

  ## ---------------------------------------------------------------------------
  ## Keep Model and Simulate tab inputs in sync (so the user can pick an
  ## estimator on either tab and it sticks).
  ## ---------------------------------------------------------------------------
  observeEvent(input$sim_estimator, {
    if (!identical(input$estimator, input$sim_estimator))
      updateSelectInput(session, "estimator", selected = input$sim_estimator)
  }, ignoreInit = TRUE)
  observeEvent(input$estimator, {
    if (!identical(input$sim_estimator, input$estimator))
      updateSelectInput(session, "sim_estimator", selected = input$estimator)
  }, ignoreInit = TRUE)
  observeEvent(input$sim_waves, {
    if (!identical(input$waves, input$sim_waves))
      updateNumericInput(session, "waves", value = input$sim_waves)
  }, ignoreInit = TRUE)
  observeEvent(input$waves, {
    if (!identical(input$sim_waves, input$waves))
      updateNumericInput(session, "sim_waves", value = input$waves)
  }, ignoreInit = TRUE)

  ## Whichever input has the latest value wins; both share the same value via
  ## the syncs above. Helpers below use the canonical names.
  cur_estimator <- reactive(input$sim_estimator %||% input$estimator %||% "CLPM")
  cur_waves     <- reactive(input$sim_waves     %||% input$waves     %||% 5)

  ## ---------------------------------------------------------------------------
  ## Estimator argument builder (uses Model-tab constraint checkboxes)
  ## ---------------------------------------------------------------------------
  build_estimator_args <- function() {
    e <- cur_estimator()
    args <- list(waves = cur_waves())

    common <- list(
      constrain_beta                 = isTRUE(input$constrain_beta),
      constrain_omega                = isTRUE(input$constrain_omega),
      constrain_residual_variances   = isTRUE(input$constrain_resvar),
      constrain_residual_covariances = isTRUE(input$constrain_rescov)
    )

    if (e == "CLPM") {
      args <- c(args, common,
                list(estimate_means = isTRUE(input$estimate_means),
                     start_values   = isTRUE(input$start_values)))
    } else if (e == "RICLPM") {
      args <- c(args, common, list(start_values = isTRUE(input$start_values)))
    } else if (e == "ALT") {
      args <- c(args, common)
    } else if (e == "LGM") {
      args <- c(args, list(
        variable_type                  = input$lgm_variable_type %||% "bivariate",
        constrain_residual_variances   = isTRUE(input$constrain_resvar),
        constrain_residual_covariances = isTRUE(input$constrain_rescov),
        estimate_quadratic             = isTRUE(input$lgm_quadratic),
        start_values                   = isTRUE(input$start_values)))
    } else if (e == "LCMSR") {
      args <- c(args, common, list(estimate_quadratic = isTRUE(input$lgm_quadratic)))
    } else if (e == "LCHANGE") {
      args <- c(args, list(
        variable_type             = input$lchange_variable_type %||% "bivariate",
        constrain_beta            = isTRUE(input$constrain_beta),
        constrain_omega           = isTRUE(input$constrain_omega),
        estimate_constant_change  = isTRUE(input$lcs_constant_change),
        estimate_change_to_change = isTRUE(input$lcs_change_to_change)))
    } else if (e == "BB") {
      args <- c(args, list(
        x_effect                     = input$bb_x_effect %||% "lagged",
        x_autoregression             = isTRUE(input$bb_x_autoreg),
        y_effect_on_x                = isTRUE(input$bb_y_on_x),
        constrain_coefficients       = isTRUE(input$bb_constrain_coef),
        constrain_residual_variances = isTRUE(input$constrain_resvar)))
    } else if (e == "TSO") {
      args <- c(args, list(
        constrain_state_variances   = isTRUE(input$tso_constrain_state_var),
        constrain_state_covariances = isTRUE(input$tso_constrain_state_cov)))
    }
    args
  }

  generate_syntax <- function() {
    e <- cur_estimator()
    if (!e %in% syntax_estimators) {
      return(paste0("# Static lavaan syntax not applicable for estimator: ", e,
                    "\n# (", e, " is a runtime estimator: fit per-trial only.)"))
    }
    fn <- switch(e,
      "CLPM"    = crossLagR::estimateCLPM,    "RICLPM"  = crossLagR::estimateRICLPM,
      "ALT"     = crossLagR::estimateALT,     "LGM"     = crossLagR::estimateLGM,
      "LCMSR"   = crossLagR::estimateLCMSR,   "LCHANGE" = crossLagR::estimateLChange,
      "BB"      = crossLagR::estimateBollen_and_Brand, "TSO" = crossLagR::estimateTSO)
    args <- build_estimator_args()
    args <- args[intersect(names(args), names(formals(fn)))]
    do.call(fn, args)
  }

  ## Cache the syntax so syntax preview + single fit always agree
  init_syntax <- reactive({
    tryCatch(generate_syntax(),
             error = function(e) paste("# Error generating syntax:", e$message))
  })
  observe({ values$init_syntax <- init_syntax() })

  ## ---------------------------------------------------------------------------
  ## Variances + ICC
  ## ---------------------------------------------------------------------------
  ## ICC = var_between / (var_between + var_within), so given ICC and within,
  ## var_between = within * icc / (1 - icc). We clamp icc to [0.02, 0.99]:
  ##   - upper bound avoids div-by-zero at icc = 1
  ##   - lower bound at 0.02 keeps the implied trait variance large enough
  ##     that lavaan's simulateData() can start the optimizer (anything
  ##     smaller and lav_start_check_cov throws "starting values imply a
  ##     correlation larger than 1" or "Sigma is not positive definite",
  ##     because the trait factor's implied covariance is essentially zero
  ##     and the simRICLPM implied matrix becomes singular).
  ## At icc = 0.02 the trait variance is ~2% of total — small enough to
  ## reproduce the "near-CLPM" condition without breaking the simulator.
  icc_to_between <- function(icc, within) {
    icc <- pmin(pmax(as.numeric(icc), 0.02), 0.99)
    within * icc / (1 - icc)
  }

  derived_between_x <- reactive({
    icc_to_between(input$icc_x %||% 0.5, input$variance_p %||% 0.5)
  })
  derived_between_y <- reactive({
    icc_to_between(input$icc_y %||% 0.5, input$variance_q %||% 0.5)
  })

  output$derived_between_text <- renderText({
    sprintf("→ between var X = %.3f (icc=%.2f · within=%.2f) | between var Y = %.3f (icc=%.2f · within=%.2f)",
            derived_between_x(), input$icc_x %||% 0.5, input$variance_p %||% 0.5,
            derived_between_y(), input$icc_y %||% 0.5, input$variance_q %||% 0.5)
  })

  ## ---------------------------------------------------------------------------
  ## Dynamic-spec helpers
  ## ---------------------------------------------------------------------------
  ## Build seq(min, max, by=step) defensively; if any input is NA / step ≤ 0 /
  ## max < min, fall back to a single-value vector.
  safe_seq <- function(lo, hi, by, default) {
    lo <- suppressWarnings(as.numeric(lo))
    hi <- suppressWarnings(as.numeric(hi))
    by <- suppressWarnings(as.numeric(by))
    if (is.na(lo) || is.na(hi)) return(default)
    if (isTRUE(all.equal(lo, hi))) return(lo)
    if (is.na(by) || by <= 0 || hi < lo) return(lo)
    seq(lo, hi, by = by)
  }

  ## Read a single dyn_range_row trio into a numeric vector.
  dyn_vec <- function(id, default) {
    safe_seq(input[[paste0(id, "_min")]],
             input[[paste0(id, "_max")]],
             input[[paste0(id, "_step")]],
             default = default)
  }

  ## ---------------------------------------------------------------------------
  ## Parameter grid — single / sweep / custom
  ## ---------------------------------------------------------------------------
  ## All downstream code expects the legacy column names variance_between_x /
  ## variance_between_y (used by simRICLPM / monteCarloLavaan). Convert from
  ## the ICC-parameterized sidebar values here.
  base_param_list <- reactive({
    out <- list(
      stability_p        = input$stability_p,
      stability_q        = input$stability_q,
      cross_p            = input$cross_p,
      cross_q            = input$cross_q,
      variance_p         = input$variance_p,
      variance_q         = input$variance_q,
      variance_between_x = derived_between_x(),
      variance_between_y = derived_between_y(),
      cov_pq             = input$cov_pq,
      icc_x              = input$icc_x %||% 0.5,
      icc_y              = input$icc_y %||% 0.5
    )
    ## Inject DGP-specific parameters into the param grid so each MC cell
    ## carries the right knobs. monteCarloLavaan reads e.g. params$confounder_p
    ## from these columns; without them the worker silently falls back to the
    ## hard-coded defaults and the user's sidebar inputs have no effect.
    if (isTRUE(input$data_generation == "clpmu")) {
      out$confounder_p         <- input$confounder_p         %||% 0.3
      out$confounder_q         <- input$confounder_q         %||% 0.3
      out$confounder_variance  <- input$confounder_variance  %||% 1
      out$confounder_stability <- input$confounder_stability %||% 0.4
    }
    out
  })

  param_grid <- reactive({
    mode <- input$grid_mode %||% "single"
    base <- base_param_list()

    if (mode == "single") {
      return(as.data.frame(do.call(expand.grid, base)))
    }

    if (mode == "dynamic") {
      ## Build per-parameter vectors from the dyn_range_row inputs. Each may
      ## be a scalar (min == max) or a sequence (min < max, step > 0).
      icc_x_vec <- dyn_vec("dyn_icc_x", default = input$icc_x %||% 0.5)
      icc_y_vec <- dyn_vec("dyn_icc_y", default = input$icc_y %||% 0.5)
      ar_x_vec  <- dyn_vec("dyn_ar_x",  default = input$stability_p)
      ar_y_vec  <- dyn_vec("dyn_ar_y",  default = input$stability_q)
      cl_xy_vec <- dyn_vec("dyn_cl_xy", default = input$cross_q)
      cl_yx_vec <- dyn_vec("dyn_cl_yx", default = input$cross_p)

      ## Cross-product. Note: internal names are stability_p/q (ar) and
      ## cross_q/p (cl_xy/cl_yx). variance_between_x/y are derived from ICC
      ## per row so they remain consistent.
      eg <- expand.grid(
        stability_p = ar_x_vec,
        stability_q = ar_y_vec,
        cross_q     = cl_xy_vec,   ## cl_xy (X→Y) → cross_q in legacy naming
        cross_p     = cl_yx_vec,   ## cl_yx (Y→X) → cross_p
        icc_x       = icc_x_vec,
        icc_y       = icc_y_vec,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      ## Add the held-fixed sidebar values
      eg$variance_p          <- input$variance_p
      eg$variance_q          <- input$variance_q
      eg$cov_pq              <- input$cov_pq %||% 0
      eg$variance_between_x  <- icc_to_between(eg$icc_x, eg$variance_p)
      eg$variance_between_y  <- icc_to_between(eg$icc_y, eg$variance_q)
      return(eg)
    }

    as.data.frame(do.call(expand.grid, base))
  })

  ## Live grid-size preview inside the dynamic-spec panel
  output$dyn_grid_preview <- renderText({
    if ((input$grid_mode %||% "single") != "dynamic") return("")
    spec <- list(
      `ICC X` = dyn_vec("dyn_icc_x", input$icc_x %||% 0.5),
      `ICC Y` = dyn_vec("dyn_icc_y", input$icc_y %||% 0.5),
      ar_x    = dyn_vec("dyn_ar_x",  input$stability_p),
      ar_y    = dyn_vec("dyn_ar_y",  input$stability_q),
      cl_xy   = dyn_vec("dyn_cl_xy", input$cross_q),
      cl_yx   = dyn_vec("dyn_cl_yx", input$cross_p)
    )
    counts <- vapply(spec, length, integer(1))
    total  <- prod(counts)
    paste0("Grid: ",
           paste(sprintf("%s=%d", names(counts), counts), collapse = " × "),
           sprintf(" = %s cells", format(total, big.mark = ",")))
  })

  output$grid_size_info <- renderText({
    grid <- param_grid()
    sprintf("Grid: %d cell%s × %d trials = %s total fits",
            nrow(grid), if (nrow(grid) == 1) "" else "s",
            input$trials,
            format(nrow(grid) * (input$trials %||% 1), big.mark = ","))
  })

  ## ---------------------------------------------------------------------------
  ## Model tab outputs
  ## ---------------------------------------------------------------------------
  output$dag_title <- renderText({
    nm <- names(estimator_choices)[match(cur_estimator(), estimator_choices)]
    sprintf("Path diagram: %s (T = %d)", nm %||% cur_estimator(), cur_waves())
  })
  output$model_dag <- DiagrammeR::renderGrViz({
    dot <- get_dag_dot(cur_estimator(), waves = cur_waves())
    if (is.null(dot)) {
      DiagrammeR::grViz("digraph G { label='(no DAG)'; node[shape=plaintext] x [label='']; }")
    } else {
      DiagrammeR::grViz(dot)
    }
  })
  output$model_syntax <- renderText({ values$init_syntax %||% "" })

  ## "Model ready" indicator in the Simulate-tab sidebar
  output$ready_text <- renderText({
    nm <- names(estimator_choices)[match(cur_estimator(), estimator_choices)]
    sprintf("%s, T=%d  (auto-prepared)", nm %||% cur_estimator(), cur_waves())
  })

  ## Preview-model button: jump to the Model tab so user can inspect the
  ## syntax + DAG without leaving the workflow.
  observeEvent(input$preview_model, {
    nav_select(id = "tabs", selected = "Model", session = session)
  })

  ## ---------------------------------------------------------------------------
  ## Data simulation dispatcher (Single fit only — Monte Carlo data sim happens
  ## inside the callr worker via run_mc_sims → monteCarloLavaan).
  ##
  ## Routes the consolidated DGP codes:
  ##   "clpmu" → simCLPMu (parameterized — only path with a confounder)
  ##   everything else → simFromSyntax (works for every unified-label
  ##                                    estimator; numerically robust at
  ##                                    extreme ICCs)
  ## ---------------------------------------------------------------------------
  simulate_dgp_dispatch <- function(dgp, waves, n, params) {
    if (dgp == "clpmu") {
      return(crossLagR::simCLPMu(
        waves = waves,
        stability_p = params$stability_p, stability_q = params$stability_q,
        cross_p = params$cross_p, cross_q = params$cross_q,
        variance_p = params$variance_p, variance_q = params$variance_q,
        cov_pq = params$cov_pq,
        confounder_p = input$confounder_p %||% 0.3,
        confounder_q = input$confounder_q %||% 0.3,
        confounder_variance = input$confounder_variance %||% 1,
        confounder_stability = input$confounder_stability %||% 0.4,
        sample.nobs = n)$data)
    }
    crossLagR::simFromSyntax(
      estimator = dgp, waves = waves, sample_size = n,
      ar_x = params$stability_p, ar_y = params$stability_q,
      cl_xy = params$cross_q, cl_yx = params$cross_p,
      d_var_x = params$variance_p, d_var_y = params$variance_q,
      d_cov_xy = params$cov_pq)
  }

  ## ---------------------------------------------------------------------------
  ## Single fit
  ## ---------------------------------------------------------------------------
  observeEvent(input$run_single_fit, {
    e <- cur_estimator()
    if (!e %in% syntax_estimators) {
      showNotification(paste("Single fit not available for", e,
                             "(runtime estimator)."), type = "warning")
      return()
    }
    one <- tryCatch({
      grid <- param_grid()
      params <- as.list(grid[1, , drop = FALSE])
      set.seed(input$seed %||% 1)
      data <- simulate_dgp_dispatch(input$data_generation, cur_waves(),
                                    input$sample_size, params)
      fit <- lavaan::lavaan(values$init_syntax, data = data,
                            meanstructure = TRUE)
      converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
      pe <- if (converged) lavaan::parameterEstimates(fit, standardized = TRUE) else NULL
      fm <- if (converged) lavaan::fitMeasures(fit, c(
        "chisq","df","pvalue","cfi","tli","rmsea","rmsea.ci.lower",
        "rmsea.ci.upper","srmr","aic","bic","npar")) else NULL
      n_obs  <- tryCatch(lavaan::lavInspect(fit, "nobs"),
                         error = function(e) nrow(data))
      n_miss <- sum(!stats::complete.cases(data))
      summary_text <- tryCatch(
        paste(capture.output(lavaan::summary(fit, fit.measures = TRUE,
                                             standardized = TRUE)),
              collapse = "\n"),
        error = function(e) paste("Summary unavailable:", e$message))
      list(ok = TRUE, converged = converged, n_obs = n_obs, n_miss = n_miss,
           fit_measures = fm, param_table = pe, summary_text = summary_text)
    }, error = function(e) list(ok = FALSE, message = e$message))

    if (!isTRUE(one$ok)) {
      showNotification(paste("Single fit failed:", one$message), type = "error")
      return()
    }
    values$single_fit <- one
    output$sim_status_text <- renderText(sprintf(
      "Single fit done. converged=%s, n_obs=%d, n_miss=%d.",
      one$converged, one$n_obs, one$n_miss))
    bslib::nav_select(id = "result_tabs", selected = "Single fit",
                      session = session)
  })

  output$single_fit_header <- renderText({
    sf <- values$single_fit
    if (is.null(sf)) return("Click 'Single fit' in the sidebar to run one fit.")
    sprintf("Estimator: %s   DGP: %s   T: %d   N: %d   seed: %d\nConverged: %s   n_obs: %d   n_miss: %d",
            cur_estimator(), input$data_generation, cur_waves(),
            input$sample_size, input$seed %||% 1,
            sf$converged, sf$n_obs, sf$n_miss)
  })
  output$single_fit_indices <- DT::renderDataTable({
    req(values$single_fit, values$single_fit$fit_measures)
    fm <- values$single_fit$fit_measures
    df <- data.frame(measure = names(fm), value = as.numeric(fm))
    DT::datatable(df, options = list(pageLength = 12, dom = "t"),
                  rownames = FALSE) |>
      DT::formatRound(columns = "value", digits = 4)
  })
  output$single_fit_keys <- DT::renderDataTable({
    req(values$single_fit, values$single_fit$param_table)
    pe <- values$single_fit$param_table
    keys <- pe[pe$label %in% c("ar_x","ar_y","cl_xy","cl_yx"),
               c("label","est","se","z","pvalue","ci.lower","ci.upper","std.all")]
    keys <- unique(keys)
    DT::datatable(keys, options = list(pageLength = 10, dom = "t"),
                  rownames = FALSE) |>
      DT::formatRound(columns = c("est","se","z","pvalue","ci.lower",
                                  "ci.upper","std.all"), digits = 4)
  })
  output$single_fit_params <- DT::renderDataTable({
    req(values$single_fit, values$single_fit$param_table)
    DT::datatable(values$single_fit$param_table,
                  options = list(scrollX = TRUE, pageLength = 25)) |>
      DT::formatRound(columns = which(vapply(values$single_fit$param_table,
                                             is.numeric, logical(1))),
                      digits = 4)
  })
  output$single_fit_lavaan_summary <- renderText({
    req(values$single_fit); values$single_fit$summary_text
  })

  ## ---------------------------------------------------------------------------
  ## Monte Carlo — callr background process with file-based progress + log
  ## ---------------------------------------------------------------------------
  pkg_root <- find_pkg_root() %||% getwd()

  observeEvent(input$run_simulation, {
    if (isTRUE(mc_bg$running)) {
      showNotification("A simulation is already running. Click Stop first.",
                       type = "warning"); return()
    }

    grid <- param_grid()
    if (nrow(grid) == 0) {
      showNotification("Parameter grid is empty.", type = "warning"); return()
    }

    progress_file <- tempfile("mc_progress_", fileext = ".txt")
    log_file      <- tempfile("mc_log_",      fileext = ".txt")
    file.create(progress_file); file.create(log_file)
    writeLines(sprintf("0/%d", nrow(grid)), progress_file)

    values$results       <- NULL
    values$summary_df    <- NULL
    mc_bg$log            <- ""
    mc_bg$start          <- Sys.time()
    mc_bg$running        <- TRUE
    mc_bg$progress_file  <- progress_file
    mc_bg$log_file       <- log_file
    mc_bg$total_cells    <- nrow(grid)
    mc_bg$done_cells     <- 0L

    est_args  <- build_estimator_args()
    dgp       <- input$data_generation
    trials    <- input$trials
    waves     <- cur_waves()
    n         <- input$sample_size

    ## In compare mode the worker iterates over ALL checked estimators per
    ## cell. Otherwise it uses the single sidebar estimator.
    if (isTRUE(input$grid_mode == "compare")) {
      estimators <- input$compare_estimators
      if (!length(estimators)) {
        showNotification("Pick at least one estimator to compare.",
                         type = "warning")
        mc_bg$running <- FALSE
        return()
      }
    } else {
      estimators <- cur_estimator()
    }

    output$sim_status_text <- renderText(
      sprintf("Launching… %d cell%s × %d trials × N=%d.",
              nrow(grid), if (nrow(grid)==1) "" else "s", trials, n))
    shinyWidgets::updateProgressBar(session = session, id = "sim_progress",
                                    value = 0, total = 100,
                                    title = sprintf("0 / %d cells", nrow(grid)))

    ## Pop the Live log tab so the user sees output streaming
    bslib::nav_select(id = "result_tabs", selected = "Live log", session = session)

    mc_bg$handle <- callr::r_bg(
      func = function(pkg_root, estimators, dgp, grid, trials, waves, n, est_args,
                      progress_file, log_file) {
        log_con <- file(log_file, open = "at")
        log_say <- function(msg) {
          line <- sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg)
          cat(line); flush.console()
          writeLines(line, log_con); flush(log_con)
        }
        if (requireNamespace("devtools", quietly = TRUE) &&
            file.exists(file.path(pkg_root, "DESCRIPTION"))) {
          suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
        } else { library(crossLagR) }

        N <- nrow(grid)
        log_say(sprintf("starting on %s DGP - %d cell(s) x %d trials x %d waves x N=%d  estimators: %s",
                        dgp, N, trials, waves, n,
                        paste(estimators, collapse = ", ")))

        all <- vector("list", N * length(estimators))
        k <- 0L
        for (i in seq_len(N)) {
          t0 <- Sys.time()
          per_cell <- vector("list", length(estimators))
          for (j in seq_along(estimators)) {
            est_j <- estimators[j]
            res_ij <- tryCatch(
              crossLagR::run_mc_sims(
                estimator = est_j, data_generation = dgp,
                param_grid = grid[i, , drop = FALSE],
                trials = trials, waves = waves, sample_size = n,
                verbose = FALSE, estimator_args = est_args),
              error = function(e) {
                log_say(sprintf("cell %d/%d [%s] FAILED: %s", i, N, est_j, e$message))
                data.frame(error_occurred = TRUE, error_message = e$message,
                            estimator = est_j, stringsAsFactors = FALSE)
              })
            ## Always ensure the estimator name is recorded (some legacy
            ## wrappers may overwrite it inconsistently).
            res_ij$estimator <- est_j
            per_cell[[j]] <- res_ij
          }
          res_i <- dplyr::bind_rows(per_cell)
          k <- k + 1L
          all[[k]] <- res_i
          dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

          ## Per-cell fit diagnostics: conv rate, mean CFI on converged fits,
          ## and counts of common lavaan pathologies. Surfaces "why is this
          ## cell slow / problematic" right in the log.
          summary_bits <- character(0)
          if ("converged" %in% names(res_i)) {
            n_conv <- sum(res_i$converged, na.rm = TRUE)
            summary_bits <- c(summary_bits,
                              sprintf("conv=%d/%d", n_conv, nrow(res_i)))
            if (n_conv > 0 && "cfi" %in% names(res_i)) {
              cfi_conv <- res_i$cfi[res_i$converged %in% TRUE]
              if (length(cfi_conv) && any(!is.na(cfi_conv)))
                summary_bits <- c(summary_bits,
                                  sprintf("CFI=%.2f", mean(cfi_conv, na.rm=TRUE)))
            }
          }
          if ("warn_iter" %in% names(res_i)) {
            n_iter <- sum(res_i$warn_iter, na.rm = TRUE)
            if (n_iter > 0)
              summary_bits <- c(summary_bits, sprintf("iter-limit=%d", n_iter))
          }
          if ("warn_npd" %in% names(res_i)) {
            n_npd <- sum(res_i$warn_npd, na.rm = TRUE)
            if (n_npd > 0)
              summary_bits <- c(summary_bits, sprintf("nonPD=%d", n_npd))
          }
          if ("warn_heywood" %in% names(res_i)) {
            n_hey <- sum(res_i$warn_heywood, na.rm = TRUE)
            if (n_hey > 0)
              summary_bits <- c(summary_bits, sprintf("Heywood=%d", n_hey))
          }
          if ("error_occurred" %in% names(res_i)) {
            n_err <- sum(res_i$error_occurred, na.rm = TRUE)
            if (n_err > 0)
              summary_bits <- c(summary_bits, sprintf("errors=%d", n_err))
          }

          ## Echo the current cell's DGP parameters so the user can correlate
          ## slow cells with extreme parameter values (e.g. high ICC).
          cell_params <- character(0)
          row1 <- grid[i, , drop = FALSE]
          for (nm in c("icc_x", "icc_y", "stability_p", "stability_q",
                       "cross_p", "cross_q")) {
            if (nm %in% names(row1) && is.numeric(row1[[nm]]))
              cell_params <- c(cell_params, sprintf("%s=%.2f", nm, row1[[nm]]))
          }

          log_say(sprintf("cell %d/%d  %.1fs  [%s] | %s",
                          i, N, dt,
                          paste(cell_params, collapse = ", "),
                          paste(summary_bits, collapse = " | ")))
          tmp <- paste0(progress_file, ".tmp")
          writeLines(sprintf("%d/%d", i, N), tmp); file.rename(tmp, progress_file)
        }
        log_say("combining cell results...")
        ## bind_rows tolerates mixed schemas (successful + error rows); rbind
        ## would error out and we'd lose all results.
        out <- tryCatch(dplyr::bind_rows(all),
                        error = function(e) {
                          log_say(sprintf("bind_rows failed: %s", e$message))
                          ## last-resort: keep only the largest dataframe
                          all[[which.max(vapply(all, function(x)
                            if (is.data.frame(x)) ncol(x) else 0L, integer(1)))]]
                        })
        log_say(sprintf("final result: %d rows x %d cols", nrow(out), ncol(out)))
        close(log_con)
        out
      },
      args = list(pkg_root = pkg_root, estimators = estimators, dgp = dgp,
                  grid = grid, trials = trials, waves = waves, n = n,
                  est_args = est_args,
                  progress_file = progress_file, log_file = log_file),
      stdout = "|", stderr = "2>&1", supervise = TRUE
    )

    mc_bg$log <- sprintf("[%s] launched background R PID %s\n",
                         format(Sys.time(), "%H:%M:%S"),
                         tryCatch(mc_bg$handle$get_pid(), error = function(e) "?"))
  })

  observeEvent(input$abort_simulation, {
    if (!isTRUE(mc_bg$running) || is.null(mc_bg$handle)) {
      showNotification("Nothing to stop.", type = "message"); return()
    }
    tryCatch(mc_bg$handle$kill(), error = function(e) NULL)
    mc_bg$running <- FALSE
    mc_bg$log <- paste0(mc_bg$log,
      sprintf("\n[%s] STOPPED by user.\n", format(Sys.time(), "%H:%M:%S")))
    mc_bg$handle <- NULL
    output$sim_status_text <- renderText("Stopped.")
    shinyWidgets::updateProgressBar(session = session, id = "sim_progress",
                                    value = 0, total = 100, title = "Stopped")
  })

  ## Poll worker: refresh log + progress bar; on completion, populate results
  observe({
    if (!isTRUE(mc_bg$running) || is.null(mc_bg$handle)) return()
    invalidateLater(700, session)

    bg <- mc_bg$handle

    log_text <- tryCatch(
      if (file.exists(mc_bg$log_file))
        paste(readLines(mc_bg$log_file, warn = FALSE), collapse = "\n") else "",
      error = function(e) "")
    if (nzchar(log_text)) mc_bg$log <- paste0(log_text, "\n")

    done <- mc_bg$done_cells
    if (!is.null(mc_bg$progress_file) && file.exists(mc_bg$progress_file)) {
      txt <- tryCatch(readLines(mc_bg$progress_file, warn = FALSE)[1],
                      error = function(e) NA_character_)
      if (!is.na(txt) && grepl("^\\d+/\\d+$", txt)) {
        done <- as.integer(strsplit(txt, "/", fixed = TRUE)[[1]][1])
      }
    }
    mc_bg$done_cells <- done

    total   <- max(1L, mc_bg$total_cells)
    pct     <- min(100L, as.integer(round(100 * done / total)))
    elapsed <- as.numeric(difftime(Sys.time(), mc_bg$start, units = "secs"))
    eta     <- if (done > 0) elapsed / done * (total - done) else NA_real_

    shinyWidgets::updateProgressBar(
      session = session, id = "sim_progress", value = pct, total = 100,
      title = sprintf("Cell %d / %d  •  %.0fs elapsed%s", done, total, elapsed,
                      if (is.finite(eta) && eta > 0) sprintf("  •  ~%.0fs left", eta) else ""))
    output$sim_status_text <- renderText(sprintf(
      "Running — %d / %d cells (%d%%). %.0fs elapsed. Click Stop to abort.",
      done, total, pct, elapsed))

    alive <- tryCatch(bg$is_alive(), error = function(e) FALSE)
    if (!alive) {
      tail_out <- tryCatch(bg$read_output(), error = function(e) "")
      if (length(tail_out) && nzchar(tail_out))
        mc_bg$log <- paste0(mc_bg$log, tail_out)
      result <- tryCatch(bg$get_result(),
                         error = function(e) {
                           mc_bg$log <- paste0(mc_bg$log,
                             sprintf("\n[ERROR] %s\n", e$message)); NULL })
      mc_bg$running <- FALSE
      mc_bg$handle  <- NULL

      if (is.null(result) || !is.data.frame(result) || nrow(result) == 0) {
        output$sim_status_text <- renderText("Failed — see log.")
        shinyWidgets::updateProgressBar(session = session, id = "sim_progress",
                                        value = 0, total = 100, title = "Failed")
        return()
      }

      result <- normalize_unified_columns(result)
      values$results    <- result
      summary <- tryCatch(summarize_mc_results(result),
                          error = function(e) {
                            mc_bg$log <- paste0(mc_bg$log,
                              sprintf("\n[server] summary failed: %s\n", e$message))
                            result  ## fall back to raw per-trial data
                          })
      ## Guarantee summary_df is non-NULL when results is set (else outputs
      ## that req(summary_df) stay blank).
      values$summary_df <- if (is.null(summary) || nrow(summary) == 0) result
                           else summary

      mc_bg$log <- paste0(mc_bg$log,
        sprintf("[%s] [server] results received: %d rows × %d cols; summary: %d rows × %d cols\n",
                format(Sys.time(), "%H:%M:%S"),
                nrow(result), ncol(result),
                nrow(values$summary_df), ncol(values$summary_df)))

      output$sim_status_text <- renderText({
        sprintf("Done. %s rows in %.1fs. Convergence: %.0f%%.",
                format(nrow(result), big.mark = ","),
                as.numeric(difftime(Sys.time(), mc_bg$start, units = "secs")),
                100 * mean(result$converged %||% rep(NA, nrow(result)), na.rm = TRUE))
      })
      shinyWidgets::updateProgressBar(session = session, id = "sim_progress",
                                      value = 100, total = 100, title = "Done")

      ## Populate plot input choices
      numeric_cols <- names(result)[vapply(result, is.numeric, logical(1))]
      grid_cols <- intersect(c("stability_p","stability_q","cross_p","cross_q",
                                "variance_p","variance_q","variance_between_x",
                                "variance_between_y","cov_pq",
                                "icc_x","icc_y"), names(result))
      updateSelectInput(session, "plot_var_x", choices = numeric_cols,
                        selected = if ("true_cl_xy" %in% numeric_cols) "true_cl_xy"
                                   else numeric_cols[1])
      updateSelectInput(session, "plot_var_y", choices = numeric_cols,
                        selected = if ("cl_xy" %in% numeric_cols) "cl_xy"
                                   else numeric_cols[min(2, length(numeric_cols))])
      updateSelectInput(session, "plot_var_color",
                        choices = c("None" = "", names(result)), selected = "")
      updateSelectInput(session, "plot_var_facet",
                        choices = c("None" = "", grid_cols), selected = "")
      updateSelectInput(session, "bias_x", choices = grid_cols,
                        selected = if (length(grid_cols)) grid_cols[1] else NULL)
      ## Comparison-tab X knob: prefer the parameter the user actually varied.
      varied <- character(0)
      for (gc in grid_cols) {
        if (length(unique(result[[gc]])) > 1) varied <- c(varied, gc)
      }
      cmp_choices <- if (length(varied)) varied else grid_cols
      updateSelectInput(session, "cmp_x", choices = cmp_choices,
                        selected = if (length(cmp_choices)) cmp_choices[1] else NULL)
      updateSelectInput(session, "cmp_ridges_y", choices = cmp_choices,
                        selected = if (length(cmp_choices)) cmp_choices[1] else NULL)
      updateSelectInput(session, "cmp_conv_x", choices = cmp_choices,
                        selected = if (length(cmp_choices)) cmp_choices[1] else NULL)

      ## Pop the user to the Summary tab so results aren't hidden behind logs.
      bslib::nav_select(id = "result_tabs", selected = "Summary",
                        session = session)
    }
  })

  ## Kill worker on session end
  session$onSessionEnded(function() {
    isolate({
      if (!is.null(mc_bg$handle))
        tryCatch(mc_bg$handle$kill(), error = function(e) NULL)
    })
  })

  output$mc_log <- renderText({
    if (nzchar(mc_bg$log)) mc_bg$log else
      "Log will stream here when a simulation is running."
  })

  ## ---------------------------------------------------------------------------
  ## Normalize column names across the package's three result schemas.
  ## The lavaan-generic wrapper (monteCarloLavaan) returns the unified labels
  ## ar_x / ar_y / cl_xy / cl_yx plus true_*. The legacy wrappers
  ## (monteCarloCLPM, monteCarloRICLPM, …) return xlag_x / xlag_y / ylag_x /
  ## ylag_y and no truth columns. To keep the plots and summary working for
  ## every estimator/DGP combo we map both schemas to the unified one here.
  ## ---------------------------------------------------------------------------
  normalize_unified_columns <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)

    rename_map <- c(ar_x  = "xlag_x",
                    ar_y  = "ylag_y",
                    cl_xy = "xlag_y",
                    cl_yx = "ylag_x")
    for (new_name in names(rename_map)) {
      old_name <- rename_map[[new_name]]
      if (old_name %in% names(df) && !new_name %in% names(df)) {
        df[[new_name]] <- df[[old_name]]
      }
    }

    truth_map <- c(true_ar_x  = "stability_p",
                   true_ar_y  = "stability_q",
                   true_cl_yx = "cross_p",
                   true_cl_xy = "cross_q")
    for (new_name in names(truth_map)) {
      old_name <- truth_map[[new_name]]
      if (old_name %in% names(df) && !new_name %in% names(df)) {
        df[[new_name]] <- df[[old_name]]
      }
    }
    df
  }

  ## ---------------------------------------------------------------------------
  ## Aggregate summary
  ## ---------------------------------------------------------------------------
  ## Hand-rolled aggregator. Earlier versions tried to compute bias inside
  ## `across(... ~ { ... .data[[true_col]] ... })`, but the `.data` pronoun
  ## doesn't see locals assigned inside the across() lambda, so
  ## `true_col` resolved as "not found" and the whole summarise threw,
  ## leaving values$summary_df with no bias columns. Computing each
  ## estimator's bias/RMSE in a plain loop sidesteps the tidy-eval landmine
  ## and is also easier to debug.
  summarize_mc_results <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    group_cols <- intersect(
      c("estimator", "dgp", "param_combo",
        "stability_p", "stability_q", "cross_p", "cross_q",
        "variance_p", "variance_q",
        "variance_between_x", "variance_between_y", "cov_pq",
        ## Keep ICC columns so plots can use them as the X-axis knob:
        "icc_x", "icc_y"),
      names(df))
    if (!length(group_cols)) return(df)
    est_cols <- intersect(c("ar_x", "ar_y", "cl_xy", "cl_yx"), names(df))
    true_map <- c(ar_x = "true_ar_x", ar_y = "true_ar_y",
                  cl_xy = "true_cl_xy", cl_yx = "true_cl_yx")
    fit_cols <- intersect(c("cfi","tli","rmsea","srmr","chisq","df","pvalue",
                            "aic","bic","n_obs","n_miss"), names(df))

    ## Means / SDs / fit-index means via the safe across() form (no .data).
    base <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::summarise(
        n_trials  = dplyr::n(),
        dplyr::across(dplyr::all_of(est_cols),
                      list(mean = ~mean(.x, na.rm = TRUE),
                           sd   = ~stats::sd(.x, na.rm = TRUE)),
                      .names = "{.col}_{.fn}"),
        dplyr::across(dplyr::all_of(fit_cols),
                      ~mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
        .groups = "drop"
      )

    ## Convergence rate — computed outside summarise so coercion is explicit.
    if ("converged" %in% names(df)) {
      tmp <- df
      tmp$.conv <- as.numeric(tmp$converged)
      conv_per_cell <- tmp %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
        dplyr::summarise(.conv = mean(.data$.conv, na.rm = TRUE),
                          .groups = "drop")
      base <- dplyr::left_join(base, conv_per_cell, by = group_cols)
      base$conv_rate <- base$.conv
      base$.conv <- NULL
    } else {
      base$conv_rate <- NA_real_
    }

    ## Bias / RMSE in a plain loop. For each unified-label estimator we look
    ## up its truth column and join the per-cell bias + RMSE back onto `base`.
    for (est in est_cols) {
      tcol <- true_map[[est]]
      bcol <- paste0(est, "_bias")
      rcol <- paste0(est, "_rmse")
      if (!tcol %in% names(df)) {
        base[[bcol]] <- NA_real_
        base[[rcol]] <- NA_real_
        next
      }
      diffs <- df[[est]] - df[[tcol]]
      tmp <- df
      tmp$.diff  <- diffs
      tmp$.diff2 <- diffs^2
      per_cell <- tmp %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
        dplyr::summarise(
          .bias = mean(.data$.diff,  na.rm = TRUE),
          .rmse = sqrt(mean(.data$.diff2, na.rm = TRUE)),
          .groups = "drop"
        )
      base <- dplyr::left_join(base, per_cell, by = group_cols)
      base[[bcol]] <- base$.bias
      base[[rcol]] <- base$.rmse
      base$.bias  <- NULL
      base$.rmse  <- NULL
    }
    base
  }

  output$mc_summary <- DT::renderDataTable({
    req(values$summary_df)
    DT::datatable(values$summary_df,
                  options = list(scrollX = TRUE, pageLength = 10),
                  caption = "Per-cell aggregates") |>
      DT::formatRound(columns = which(vapply(values$summary_df,
                                             is.numeric, logical(1))),
                      digits = 4)
  })
  output$results_table <- DT::renderDataTable({
    req(values$results)
    DT::datatable(values$results,
                  options = list(scrollX = TRUE, pageLength = 10)) |>
      DT::formatRound(columns = which(vapply(values$results, is.numeric, logical(1))),
                      digits = 4)
  })

  ## ---------------------------------------------------------------------------
  ## Plots
  ## ---------------------------------------------------------------------------
  ## Helper: render a textual ggplot canvas as a fallback message
  blank_plot <- function(msg) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = msg,
               color = "#37474F", size = 5, fontface = "italic") +
      xlim(0, 1) + ylim(0, 1) + theme_void()
  }

  output$bias_ggplot <- renderPlot({
    if (is.null(values$summary_df) || nrow(values$summary_df) == 0)
      return(blank_plot("Run a Monte Carlo to see bias by cell."))
    s <- values$summary_df
    if (is.null(input$bias_param) || !nzchar(input$bias_param))
      return(blank_plot("Pick a parameter."))

    bias_col <- paste0(input$bias_param, "_bias")
    rmse_col <- paste0(input$bias_param, "_rmse")

    if (!bias_col %in% names(s)) {
      avail <- grep("_bias$", names(s), value = TRUE)
      return(blank_plot(paste0("No bias column for '", input$bias_param, "'.\n",
                               if (length(avail))
                                 paste0("Available: ", paste(avail, collapse = ", "))
                               else
                                 "(No bias columns — the estimator did not return unified ar/cl labels.)")))
    }

    x_col <- input$bias_x
    if (is.null(x_col) || !nzchar(x_col) || !x_col %in% names(s)) {
      ## Fall back to param_combo so we always render *something*
      x_col <- if ("param_combo" %in% names(s)) "param_combo" else names(s)[1]
    }

    ## Single-cell case: render a one-point summary with the value annotated.
    ## Use a numeric x (not factor) so scale_x_continuous works without the
    ## "Discrete value supplied to a continuous scale" error that was leaving
    ## the panel blank.
    if (nrow(s) == 1) {
      v <- s[[bias_col]][1]
      r <- if (rmse_col %in% names(s)) s[[rmse_col]][1] else NA_real_
      pt <- data.frame(x = 1, y = v)
      p <- ggplot(pt, aes(x = .data$x, y = .data$y)) +
        geom_hline(yintercept = 0, color = "#90A4AE", linetype = "dashed") +
        geom_point(size = 5, color = "#37474F")
      if (!is.na(r)) {
        p <- p + geom_errorbar(aes(ymin = .data$y - r, ymax = .data$y + r),
                                width = 0.08, color = "#37474F", alpha = 0.6)
      }
      p <- p +
        annotate("text", x = 1.15, y = v,
                 label = sprintf("bias = %.4f%s", v,
                                 if (!is.na(r)) sprintf("\nRMSE = %.4f", r) else ""),
                 hjust = 0, color = "#263238", size = 4.5) +
        scale_x_continuous(limits = c(0.5, 2), breaks = NULL) +
        labs(title = sprintf("Bias of %s (single cell)", input$bias_param),
             subtitle = "Enable 'Sweep' in the sidebar to vary a parameter across cells.",
             x = NULL, y = "Bias  (mean est − true)") +
        acad_theme() +
        theme(axis.ticks.x = element_blank())
      return(p)
    }

    p <- ggplot(s, aes(x = .data[[x_col]], y = .data[[bias_col]])) +
      geom_hline(yintercept = 0, color = "#90A4AE", linetype = "dashed") +
      geom_point(size = 3, color = "#37474F") +
      geom_line(aes(group = 1), color = "#37474F", alpha = 0.5) +
      labs(title = sprintf("Bias of %s vs %s", input$bias_param, x_col),
           subtitle = sprintf("Per parameter cell. Whiskers = ±RMSE (%s).", rmse_col),
           x = x_col, y = "Bias  (mean est − true)") +
      acad_theme()
    if (rmse_col %in% names(s)) {
      p <- p + geom_errorbar(aes(ymin = .data[[bias_col]] - .data[[rmse_col]],
                                  ymax = .data[[bias_col]] + .data[[rmse_col]]),
                              width = 0.01, color = "#37474F", alpha = 0.5)
    }
    p
  })

  output$dist_ggplot <- renderPlot({
    if (is.null(values$results) || nrow(values$results) == 0)
      return(blank_plot("Run a Monte Carlo to see estimate distributions."))
    df <- values$results
    if (is.null(input$dist_param) || !nzchar(input$dist_param))
      return(blank_plot("Pick a parameter."))
    if (!input$dist_param %in% names(df)) {
      avail <- intersect(c("ar_x","ar_y","cl_xy","cl_yx"), names(df))
      return(blank_plot(paste0("Column '", input$dist_param, "' not in results.\n",
                               if (length(avail))
                                 paste0("Available: ", paste(avail, collapse = ", "))
                               else
                                 "(No unified ar/cl columns — try a different estimator.)")))
    }
    truth_col <- paste0("true_", input$dist_param)
    base <- ggplot(df, aes(x = .data[[input$dist_param]])) +
      labs(title = sprintf("Sampling distribution of %s", input$dist_param),
           subtitle = sprintf("Across %s trials × cells",
                              format(nrow(df), big.mark = ",")),
           x = input$dist_param, y = NULL) +
      acad_theme()
    p <- switch(input$dist_type,
      "density" = base + geom_density(fill = "#37474F", alpha = 0.4,
                                       color = "#37474F", linewidth = 0.6),
      "hist"    = base + geom_histogram(fill = "#37474F", color = "white",
                                         bins = 40, alpha = 0.85),
      "ridges"  = {
        if (!"param_combo" %in% names(df))
          base + geom_density(fill = "#37474F", alpha = 0.4)
        else if (requireNamespace("ggridges", quietly = TRUE)) {
          ggplot(df, aes(x = .data[[input$dist_param]],
                         y = factor(.data$param_combo),
                         fill = factor(.data$param_combo))) +
            ggridges::geom_density_ridges(alpha = 0.7, color = "white") +
            scale_fill_manual(values = grDevices::colorRampPalette(acad_palette)(
              length(unique(df$param_combo))), guide = "none") +
            labs(title = sprintf("Ridge plot of %s by cell", input$dist_param),
                 x = input$dist_param, y = "Cell") + acad_theme()
        } else {
          base + geom_density(aes(group = .data$param_combo,
                                  fill = factor(.data$param_combo)), alpha = 0.4) +
            scale_fill_manual(values = grDevices::colorRampPalette(acad_palette)(
              length(unique(df$param_combo))))
        }
      })
    if (isTRUE(input$dist_show_truth) && truth_col %in% names(df)) {
      truths <- unique(df[[truth_col]])
      truths <- truths[!is.na(truths)]
      if (length(truths))
        p <- p + geom_vline(xintercept = truths, linetype = "dashed",
                            color = "#C62828", linewidth = 0.6)
    }
    p
  })

  output$results_plot <- renderPlotly({
    if (is.null(values$results) || nrow(values$results) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers") |>
        plotly::layout(title = "Run a Monte Carlo to populate the interactive plot.",
                        xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE)))
    }
    if (is.null(input$plot_var_x) || !nzchar(input$plot_var_x) ||
        is.null(input$plot_var_y) || !nzchar(input$plot_var_y)) {
      return(plotly_empty() |>
        plotly::layout(title = "Pick X / Y in the sidebar above."))
    }
    df <- values$results
    p <- ggplot(df, aes(x = .data[[input$plot_var_x]],
                        y = .data[[input$plot_var_y]])) +
      geom_point(alpha = 0.7, size = 2, color = "#37474F") +
      acad_theme() +
      labs(title = paste("Monte Carlo:", cur_estimator()),
           subtitle = paste("DGP:", input$data_generation))
    if (!is.null(input$plot_var_color) && nzchar(input$plot_var_color)) {
      p <- p + aes(color = factor(.data[[input$plot_var_color]])) +
        scale_color_manual(values = acad_palette)
    }
    if (!is.null(input$plot_var_facet) && nzchar(input$plot_var_facet)) {
      p <- p + facet_wrap(stats::as.formula(paste("~", input$plot_var_facet)))
    }
    ggplotly(p)
  })

  ## ---------------------------------------------------------------------------
  ## Estimator comparison plot — reconstructs chapter-5 fig-sim1-bias /
  ## fig-sim1-rmse / fig-sim1-conv. One facet per AR/CL parameter, one
  ## colored line per estimator, x-axis is the swept DGP knob.
  ## ---------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  ## Helper: long-form per-(estimator × cell × parameter) frame with bias,
  ## RMSE, mean, SD, n_trials, and 95% Monte-Carlo intervals. Used by both
  ## the trajectory and heatmap views.
  ## ---------------------------------------------------------------------------
  compare_long <- reactive({
    s <- values$summary_df
    if (is.null(s) || nrow(s) == 0 || !"estimator" %in% names(s)) return(NULL)
    pcols <- c("ar_x", "ar_y", "cl_xy", "cl_yx")
    keep_meta <- intersect(c("estimator", "n_trials", "conv_rate",
                              "stability_p", "stability_q", "cross_p", "cross_q",
                              "variance_p", "variance_q",
                              "variance_between_x", "variance_between_y",
                              "cov_pq", "icc_x", "icc_y"), names(s))
    out_long <- NULL
    for (pp in pcols) {
      mc <- paste0(pp, "_mean"); bc <- paste0(pp, "_bias")
      rc <- paste0(pp, "_rmse"); sc <- paste0(pp, "_sd")
      if (!bc %in% names(s)) next
      row <- s[, keep_meta, drop = FALSE]
      row$parameter <- pp
      row$mean <- if (mc %in% names(s)) s[[mc]] else NA_real_
      row$bias <- s[[bc]]
      row$rmse <- if (rc %in% names(s)) s[[rc]] else NA_real_
      row$sd   <- if (sc %in% names(s)) s[[sc]] else NA_real_
      ## 95% Monte-Carlo interval for the cell mean: mean ± 1.96 · SD/sqrt(n)
      n_eff <- pmax(row$n_trials * (row$conv_rate %||% 1), 1)
      row$mc_lo_mean <- row$mean - qnorm(0.975) * row$sd / sqrt(n_eff)
      row$mc_hi_mean <- row$mean + qnorm(0.975) * row$sd / sqrt(n_eff)
      row$mc_lo_bias <- row$bias - qnorm(0.975) * row$sd / sqrt(n_eff)
      row$mc_hi_bias <- row$bias + qnorm(0.975) * row$sd / sqrt(n_eff)
      out_long <- if (is.null(out_long)) row else rbind(out_long, row)
    }
    if (!is.null(out_long))
      out_long$parameter <- factor(out_long$parameter, levels = pcols)
    out_long
  })

  output$compare_ggplot <- renderPlot({
    if (is.null(values$summary_df) || nrow(values$summary_df) == 0)
      return(blank_plot("Set grid mode to 'Compare estimators' and run."))
    s <- values$summary_df
    if (!"estimator" %in% names(s))
      return(blank_plot("Result has no 'estimator' column. Use grid mode = 'Compare estimators'."))

    x_col  <- input$cmp_x
    metric <- input$cmp_metric %||% "bias"

    if (is.null(x_col) || !nzchar(x_col) || !x_col %in% names(s))
      return(blank_plot("Pick an X-axis knob (must be a varied DGP parameter)."))

    d <- compare_long()
    if (is.null(d) || nrow(d) == 0)
      return(blank_plot("No bias/RMSE columns in summary."))
    if (!x_col %in% names(d))
      return(blank_plot(paste0("X-axis column '", x_col, "' not in long form.")))

    if (metric == "rmse") {
      d$.value <- d$rmse
      ylab <- "RMSE"
      ref_line <- NA_real_
    } else {
      d$.value <- d$bias
      d$.lo    <- d$mc_lo_bias
      d$.hi    <- d$mc_hi_bias
      ylab <- "Bias (mean est − true)"
      ref_line <- 0
    }

    ## Sort estimators for a consistent legend order
    est_levels <- unique(d$estimator)
    est_levels <- est_levels[order(match(est_levels,
                                          c("CLPM","RICLPM","ALT","LGM",
                                            "LCMSR","BB","LCHANGE","TSO")))]
    d$estimator <- factor(d$estimator, levels = est_levels)

    n_est <- length(est_levels)
    palette <- grDevices::colorRampPalette(acad_palette)(n_est)

    p <- ggplot(d, aes(x = .data[[x_col]], y = .data$.value,
                       color = .data$estimator, shape = .data$estimator,
                       group = .data$estimator)) +
      facet_wrap(~ parameter, scales = "free_y", nrow = 1) +
      scale_color_manual(values = palette) +
      scale_shape_manual(values = c(16, 17, 15, 18, 7, 8, 9, 10)[seq_len(n_est)]) +
      labs(title = switch(metric, bias = "Bias by estimator",
                                  rmse = "RMSE by estimator"),
           subtitle = sprintf("Facets = AR/CL parameters  |  X-axis: %s  |  %d cells × %d trials/cell (median)",
                              x_col, length(unique(d[[x_col]])),
                              stats::median(d$n_trials, na.rm = TRUE)),
           x = x_col, y = ylab,
           color = "Estimator", shape = "Estimator") +
      acad_theme() +
      theme(panel.spacing.x = grid::unit(1.2, "lines"),
            legend.key.width = grid::unit(2, "lines"),
            strip.text = element_text(size = 13, face = "bold"))

    if (!is.na(ref_line)) {
      p <- p + geom_hline(yintercept = ref_line, color = "#90A4AE",
                          linetype = "dashed", linewidth = 0.6)
    }
    ## 95% MC intervals (only for bias; RMSE uncertainty is harder to derive)
    if (metric == "bias" && isTRUE(input$cmp_show_ci %||% TRUE)) {
      p <- p + geom_ribbon(aes(ymin = .data$.lo, ymax = .data$.hi,
                                fill = .data$estimator),
                            alpha = 0.12, color = NA, show.legend = FALSE) +
                scale_fill_manual(values = palette)
    }
    p <- p +
      geom_line(linewidth = 1.1, alpha = 0.9) +
      geom_point(size = 3.4, stroke = 1.2)

    if (isTRUE(input$cmp_drop_outliers)) {
      p <- p + coord_cartesian(ylim = c(-0.3, 0.3))
    }
    p
  })
  outputOptions(output, "compare_ggplot", suspendWhenHidden = FALSE)

  ## ---------------------------------------------------------------------------
  ## Compare-mode RIDGE plot — chapter-5 fig-sim1-ridges style. Rows =
  ## estimator, columns = parameter, ridges stacked vertically by the swept
  ## DGP knob; fill encodes mean |estimate − truth|.
  ## ---------------------------------------------------------------------------
  output$compare_ridges <- renderPlot({
    if (is.null(values$results) || nrow(values$results) == 0)
      return(blank_plot("Run a 'Compare estimators' Monte Carlo first."))
    df <- values$results
    if (!"estimator" %in% names(df))
      return(blank_plot("No 'estimator' column — use grid mode 'Compare estimators'."))

    yvar <- input$cmp_ridges_y
    if (is.null(yvar) || !nzchar(yvar) || !yvar %in% names(df))
      return(blank_plot("Pick a Y-axis (DGP knob) for the ridges."))

    pcols    <- c("ar_x", "ar_y", "cl_xy", "cl_yx")
    true_map <- c(ar_x = "true_ar_x", ar_y = "true_ar_y",
                  cl_xy = "true_cl_xy", cl_yx = "true_cl_yx")
    pcols <- intersect(pcols, names(df))
    if (!length(pcols))
      return(blank_plot("No unified ar/cl columns in results."))

    focus <- input$cmp_ridges_facet %||% "all"
    if (focus != "all" && focus %in% pcols) pcols <- focus

    long_pieces <- lapply(pcols, function(pp) {
      tc <- true_map[[pp]]
      data.frame(
        estimator = df$estimator, parameter = pp,
        estimate  = df[[pp]],
        truth     = if (tc %in% names(df)) df[[tc]] else NA_real_,
        yval      = df[[yvar]],
        stringsAsFactors = FALSE)
    })
    long <- do.call(rbind, long_pieces)
    long <- long[stats::complete.cases(long$estimator, long$parameter, long$estimate), ]
    if (!nrow(long))
      return(blank_plot("No usable rows for the ridge plot."))

    ## Per-parameter clipping: drop estimates that fall outside the true
    ## value ± k. This removes runaway tails (e.g. LCS ar_x in the ±5
    ## range) that otherwise smear every panel's x-axis.
    k <- input$cmp_ridges_clip %||% 2
    if (is.finite(k) && k > 0) {
      keep <- with(long,
        is.na(truth) | (estimate >= truth - k & estimate <= truth + k))
      long <- long[keep, , drop = FALSE]
    }

    ## Mean |estimate − truth| per (estimator × parameter × yval) for fill
    mad <- aggregate(abs(long$estimate - long$truth),
                      by = list(estimator = long$estimator,
                                parameter = long$parameter,
                                yval      = long$yval),
                      FUN = function(x) mean(x, na.rm = TRUE))
    names(mad)[4] <- "mad"
    long <- merge(long, mad, by = c("estimator","parameter","yval"),
                  all.x = TRUE)

    long$yval_f <- factor(round(long$yval, 3))
    long$parameter <- factor(long$parameter,
                              levels = c("ar_x","ar_y","cl_xy","cl_yx"))

    ## Estimator order (consistent across views)
    est_levels <- unique(long$estimator)
    est_levels <- est_levels[order(match(est_levels,
                                          c("CLPM","RICLPM","ALT","LGM",
                                            "LCMSR","BB","LCHANGE","TSO")))]
    long$estimator <- factor(long$estimator, levels = est_levels)

    has_ridges <- requireNamespace("ggridges", quietly = TRUE)
    midpoint <- stats::quantile(long$mad, 0.5, na.rm = TRUE)
    p <- ggplot(long, aes(x = .data$estimate, y = .data$yval_f,
                           fill = .data$mad)) +
      scale_fill_gradient2(low = "#2E7D32", mid = "#FFFDE7", high = "#C62828",
                            midpoint = midpoint,
                            name = "mean\n|est − truth|") +
      labs(title = "Sampling distributions of AR/CL estimates",
           subtitle = sprintf("Rows = estimator, columns = parameter, ridges stacked by %s  |  clipped to ±%g of truth",
                              yvar, k),
           x = "Estimate", y = yvar) +
      acad_theme() +
      theme(strip.placement   = "outside",
            strip.text         = element_text(size = 12, face = "bold",
                                              color = "#263238"),
            strip.background.x = element_rect(fill = "#ECEFF1", color = NA),
            strip.background.y = element_rect(fill = "#CFD8DC", color = NA),
            panel.spacing      = grid::unit(0.6, "lines"),
            legend.position    = "right",
            legend.key.height  = grid::unit(1.6, "lines"))

    p <- if (has_ridges) {
      p + ggridges::geom_density_ridges(alpha = 0.88, color = "white",
                                         scale = 1.0, rel_min_height = 0.005)
    } else {
      p + geom_violin(alpha = 0.85, color = "white",
                       draw_quantiles = 0.5)
    }

    ## Truth lines — heavier so they're visible behind ridges
    truth_lines <- unique(long[, c("parameter", "truth")])
    truth_lines <- truth_lines[!is.na(truth_lines$truth), ]
    if (nrow(truth_lines))
      p <- p + geom_vline(data = truth_lines,
                          aes(xintercept = .data$truth),
                          color = "#263238", linewidth = 0.9,
                          linetype = "longdash")

    if (length(unique(long$parameter)) > 1) {
      p <- p + facet_grid(estimator ~ parameter,
                          scales = "free_x", switch = "y")
    } else {
      p <- p + facet_wrap(~ estimator, ncol = 1, strip.position = "right",
                          scales = "free_x")
    }
    p
  })
  outputOptions(output, "compare_ridges", suspendWhenHidden = FALSE)

  ## ---------------------------------------------------------------------------
  ## Compare-mode PERFORMANCE HEATMAP — at-a-glance view: rows = estimator,
  ## columns = (parameter × cell value). Fill = bias, |bias|, or RMSE.
  ## ---------------------------------------------------------------------------
  output$compare_heatmap <- renderPlot({
    d <- compare_long()
    if (is.null(d) || nrow(d) == 0)
      return(blank_plot("Run a 'Compare estimators' Monte Carlo first."))
    metric <- input$cmp_heat_metric %||% "abs_bias"

    ## Pick the varied DGP knob as the cell-id column
    grid_cols <- intersect(c("stability_p","stability_q","cross_p","cross_q",
                              "variance_p","variance_q","variance_between_x",
                              "variance_between_y","cov_pq","icc_x","icc_y"),
                            names(d))
    varied <- grid_cols[vapply(grid_cols, function(g)
                                 length(unique(d[[g]])) > 1, logical(1))]
    cell_col <- if (length(varied)) varied[1] else
                  if (length(grid_cols)) grid_cols[1] else "parameter"

    d$.cell <- factor(round(d[[cell_col]], 3))
    d$.val  <- switch(metric,
                       bias     = d$bias,
                       abs_bias = abs(d$bias),
                       rmse     = d$rmse)
    d$.lbl  <- sprintf("%.3f", d$.val)

    est_levels <- unique(d$estimator)
    est_levels <- est_levels[order(match(est_levels,
                                          c("CLPM","RICLPM","ALT","LGM",
                                            "LCMSR","BB","LCHANGE","TSO")))]
    d$estimator <- factor(d$estimator, levels = rev(est_levels))

    if (metric == "bias") {
      lim <- max(abs(d$.val), na.rm = TRUE)
      fill_layer <- scale_fill_gradient2(low = "#2E7D32", mid = "#FFFDE7",
                                          high = "#C62828", midpoint = 0,
                                          limits = c(-lim, lim),
                                          name = "Bias")
    } else {
      fill_layer <- scale_fill_gradient(low = "#FFFDE7", high = "#C62828",
                                         name = if (metric == "rmse") "RMSE" else "|Bias|")
    }

    p <- ggplot(d, aes(x = .data$.cell, y = .data$estimator, fill = .data$.val)) +
      geom_tile(color = "white", linewidth = 1.0) +
      facet_wrap(~ parameter, nrow = 1) +
      fill_layer +
      labs(title = sprintf("Performance heatmap (%s)",
                            switch(metric, bias = "signed bias",
                                            abs_bias = "absolute bias",
                                            rmse = "RMSE")),
           subtitle = sprintf("Rows = estimator, columns = %s  |  Facets = AR/CL parameters",
                              cell_col),
           x = cell_col, y = NULL) +
      acad_theme() +
      theme(strip.text = element_text(size = 13, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank())

    if (isTRUE(input$cmp_heat_label %||% TRUE)) {
      p <- p + geom_text(aes(label = .data$.lbl), size = 3.0,
                          color = "#263238")
    }
    p
  })
  outputOptions(output, "compare_heatmap", suspendWhenHidden = FALSE)

  ## ---------------------------------------------------------------------------
  ## Compare-mode CONVERGENCE plot — bar chart of conv_rate per (estimator ×
  ## cell), colored red below 100% to flag dropouts.
  ## ---------------------------------------------------------------------------
  output$compare_conv <- renderPlot({
    s <- values$summary_df
    if (is.null(s) || nrow(s) == 0 || !"estimator" %in% names(s) ||
        !"conv_rate" %in% names(s))
      return(blank_plot("Run a 'Compare estimators' MC; conv_rate not yet available."))

    x_col <- input$cmp_conv_x %||% input$cmp_x
    if (is.null(x_col) || !nzchar(x_col) || !x_col %in% names(s)) {
      grid_cols <- intersect(c("stability_p","stability_q","cross_p","cross_q",
                                "variance_p","variance_q","variance_between_x",
                                "variance_between_y","cov_pq","icc_x","icc_y"),
                              names(s))
      varied <- grid_cols[vapply(grid_cols, function(g)
                                   length(unique(s[[g]])) > 1, logical(1))]
      x_col <- if (length(varied)) varied[1] else
               if (length(grid_cols)) grid_cols[1] else "estimator"
    }

    est_levels <- unique(s$estimator)
    est_levels <- est_levels[order(match(est_levels,
                                          c("CLPM","RICLPM","ALT","LGM",
                                            "LCMSR","BB","LCHANGE","TSO")))]
    s$estimator <- factor(s$estimator, levels = est_levels)
    palette <- grDevices::colorRampPalette(acad_palette)(length(est_levels))

    ggplot(s, aes(x = factor(round(.data[[x_col]], 3)),
                  y = .data$conv_rate,
                  fill = .data$estimator)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      geom_hline(yintercept = 1, color = "#90A4AE",
                  linetype = "dashed", linewidth = 0.5) +
      scale_fill_manual(values = palette) +
      scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, by = 0.2)) +
      labs(title = "Convergence rate by estimator and cell",
           subtitle = sprintf("Dashed line = 100%%. X-axis: %s.", x_col),
           x = x_col, y = "Convergence rate", fill = "Estimator") +
      acad_theme() +
      theme(legend.position = "right")
  })
  outputOptions(output, "compare_conv", suspendWhenHidden = FALSE)

  ## ---------------------------------------------------------------------------
  ## Compare-mode side-by-side TABLE — one row per (estimator × cell), with
  ## bias / RMSE / mean for each unified-label parameter and conv_rate.
  ## ---------------------------------------------------------------------------
  output$compare_table <- DT::renderDataTable({
    s <- values$summary_df
    if (is.null(s) || nrow(s) == 0 || !"estimator" %in% names(s))
      return(DT::datatable(data.frame(message = "Run a 'Compare estimators' MC first."),
                           options = list(dom = "t"), rownames = FALSE))

    ## Pick the varied DGP knob to keep as a 'cell' column (else use param_combo).
    grid_cols <- intersect(c("stability_p","stability_q","cross_p","cross_q",
                              "variance_p","variance_q","variance_between_x",
                              "variance_between_y","cov_pq","icc_x","icc_y"),
                            names(s))
    varied <- grid_cols[vapply(grid_cols, function(g)
                                 length(unique(s[[g]])) > 1, logical(1))]
    cell_col <- if (length(varied)) varied[1] else
                  if ("param_combo" %in% names(s)) "param_combo" else NULL

    keep_cols <- c("estimator", cell_col, "n_trials", "conv_rate",
                   intersect(c("ar_x_mean", "ar_x_bias", "ar_x_rmse",
                               "ar_y_mean", "ar_y_bias", "ar_y_rmse",
                               "cl_xy_mean","cl_xy_bias","cl_xy_rmse",
                               "cl_yx_mean","cl_yx_bias","cl_yx_rmse"),
                              names(s)))
    keep_cols <- intersect(keep_cols, names(s))
    out <- s[, keep_cols, drop = FALSE]
    if (!is.null(cell_col))
      out <- out[order(out[[cell_col]], out$estimator), , drop = FALSE]

    num_cols <- which(vapply(out, is.numeric, logical(1)))
    bias_cols <- grep("_bias$", names(out))

    dt <- DT::datatable(
      out,
      options = list(scrollX = TRUE, pageLength = 25, dom = "tip"),
      rownames = FALSE,
      caption = sprintf("Side-by-side comparison (%d estimator%s × %d cell%s)",
                        length(unique(out$estimator)),
                        if (length(unique(out$estimator)) == 1) "" else "s",
                        if (!is.null(cell_col)) length(unique(out[[cell_col]])) else 1,
                        if (!is.null(cell_col) && length(unique(out[[cell_col]])) != 1)
                          "s" else "")
    )
    dt <- DT::formatRound(dt, columns = num_cols, digits = 4)
    if (length(bias_cols)) {
      ## Diverging color: large |bias| = red, near-zero = white
      max_b <- max(abs(unlist(out[, bias_cols], use.names = FALSE)),
                   na.rm = TRUE)
      if (is.finite(max_b) && max_b > 0) {
        brks <- seq(-max_b, max_b, length.out = 11)
        clrs <- colorRampPalette(c("#2E7D32", "#FFFFFF", "#C62828"))(12)
        dt <- DT::formatStyle(dt, columns = bias_cols,
                              backgroundColor = DT::styleInterval(brks, clrs))
      }
    }
    dt
  })
  outputOptions(output, "compare_table", suspendWhenHidden = FALSE)

  output$download_results <- downloadHandler(
    filename = function() {
      paste0("mc_", cur_estimator(), "_", input$data_generation, "_",
             Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$results)
      write.csv(values$results, file, row.names = FALSE)
    }
  )

  ## ---------------------------------------------------------------------------
  ## Force-render the result outputs even when their tab is not currently
  ## visible. bslib's navset_card_tab suspends hidden outputs by default; that
  ## means when MC completes and we navigate to Summary, the plot outputs on
  ## other tabs (bias_ggplot, dist_ggplot, results_plot) still haven't been
  ## evaluated. When the user clicks the tab they see a blank panel because
  ## the input choices (bias_x, etc.) were updated *while* the output was
  ## suspended. Eager rendering forces evaluation as soon as their deps land.
  ## ---------------------------------------------------------------------------
  outputOptions(output, "mc_summary",         suspendWhenHidden = FALSE)
  outputOptions(output, "bias_ggplot",        suspendWhenHidden = FALSE)
  outputOptions(output, "dist_ggplot",        suspendWhenHidden = FALSE)
  outputOptions(output, "results_plot",       suspendWhenHidden = FALSE)
  outputOptions(output, "results_table",      suspendWhenHidden = FALSE)
  outputOptions(output, "single_fit_indices", suspendWhenHidden = FALSE)
  outputOptions(output, "single_fit_keys",    suspendWhenHidden = FALSE)
  outputOptions(output, "single_fit_params",  suspendWhenHidden = FALSE)
}
