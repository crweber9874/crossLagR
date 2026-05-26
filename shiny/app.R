library(shiny)
library(MASS)
library(lavaan)
library(ggplot2)
library(tidyr)
library(dplyr)

# ============================================================================
# Helper Functions: Shared
# ============================================================================

# Build RI-CLPM data-generating model syntax for lavaan
build_riclpm_dgp <- function(waves, beta_x, beta_y, omega_xy, omega_yx,
                              var_p, var_q, cov_pq,
                              var_BX, var_BY, cov_BXBY) {
  m <- ""
  m <- paste0(m, "BX =~ 1*x1")
  for (w in 2:waves) m <- paste0(m, " + 1*x", w)
  m <- paste0(m, "\nBY =~ 1*y1")
  for (w in 2:waves) m <- paste0(m, " + 1*y", w)

  for (w in 1:waves) {
    m <- paste0(m, "\np", w, " =~ 1*x", w)
    m <- paste0(m, "\nq", w, " =~ 1*y", w)
  }

  m <- paste0(m, "\nBX ~ 0*1\nBY ~ 0*1")
  for (w in 1:waves) {
    m <- paste0(m, "\np", w, " ~ 0*1")
    m <- paste0(m, "\nq", w, " ~ 0*1")
  }
  for (w in 1:waves) {
    m <- paste0(m, "\nx", w, " ~ 0*1")
    m <- paste0(m, "\ny", w, " ~ 0*1")
  }

  for (w in 2:waves) {
    m <- paste0(m, "\np", w, " ~ ", beta_x, "*p", w - 1, " + ", omega_yx, "*q", w - 1)
    m <- paste0(m, "\nq", w, " ~ ", beta_y, "*q", w - 1, " + ", omega_xy, "*p", w - 1)
  }

  for (w in 1:waves) {
    m <- paste0(m, "\np", w, " ~~ ", var_p, "*p", w)
    m <- paste0(m, "\nq", w, " ~~ ", var_q, "*q", w)
  }
  for (w in 1:waves) {
    m <- paste0(m, "\np", w, " ~~ ", cov_pq, "*q", w)
  }
  for (w in 1:waves) {
    m <- paste0(m, "\nx", w, " ~~ 0*x", w)
    m <- paste0(m, "\ny", w, " ~~ 0*y", w)
  }

  m <- paste0(m, "\nBX ~~ ", var_BX, "*BX")
  m <- paste0(m, "\nBY ~~ ", var_BY, "*BY")
  m <- paste0(m, "\nBX ~~ ", cov_BXBY, "*BY")
  m <- paste0(m, "\np1 ~~ 0*BX\np1 ~~ 0*BY\nq1 ~~ 0*BX\nq1 ~~ 0*BY")

  return(m)
}

# Build standard CLPM estimation model syntax (no trait factors)
build_clpm_estimator <- function(waves) {
  m <- ""
  for (w in 1:waves) {
    m <- paste0(m, "p", w, " =~ 1*y", w, "\n")
    m <- paste0(m, "q", w, " =~ 1*x", w, "\n")
  }

  m <- paste0(m, "p1 ~ 1\nq1 ~ 1\n")
  for (w in 2:waves) {
    m <- paste0(m, "p", w, " ~ 0*1\n")
    m <- paste0(m, "q", w, " ~ 0*1\n")
  }
  for (w in 1:waves) {
    m <- paste0(m, "y", w, " ~ 0*1\n")
    m <- paste0(m, "x", w, " ~ 0*1\n")
  }

  for (w in 2:waves) {
    m <- paste0(m, "p", w, " ~ beta_y*p", w - 1, " + omega_xy*q", w - 1, "\n")
    m <- paste0(m, "q", w, " ~ beta_x*q", w - 1, " + omega_yx*p", w - 1, "\n")
  }

  m <- paste0(m, "p1 ~~ var_y1*p1\nq1 ~~ var_x1*q1\n")
  for (w in 2:waves) {
    m <- paste0(m, "p", w, " ~~ var_y*p", w, "\n")
    m <- paste0(m, "q", w, " ~~ var_x*q", w, "\n")
  }

  m <- paste0(m, "p1 ~~ cov_xy1*q1\n")
  for (w in 2:waves) {
    m <- paste0(m, "p", w, " ~~ cov_xy*q", w, "\n")
  }

  for (w in 1:waves) {
    m <- paste0(m, "y", w, " ~~ 0*y", w, "\n")
    m <- paste0(m, "x", w, " ~~ 0*x", w, "\n")
  }

  return(m)
}


# ============================================================================
# Helper Functions: Tab 1 (Trait Bias)
# ============================================================================

run_one_trial_clpm <- function(waves, sample_size,
                                beta_x, beta_y, omega_xy, omega_yx,
                                var_p, var_q, cov_pq,
                                var_BX, var_BY, cov_BXBY) {

  dgp_syntax <- build_riclpm_dgp(
    waves = waves, beta_x = beta_x, beta_y = beta_y,
    omega_xy = omega_xy, omega_yx = omega_yx,
    var_p = var_p, var_q = var_q, cov_pq = cov_pq,
    var_BX = var_BX, var_BY = var_BY, cov_BXBY = cov_BXBY
  )
  est_syntax <- build_clpm_estimator(waves = waves)

  tryCatch({
    dat <- lavaan::simulateData(dgp_syntax, sample.nobs = sample_size,
                                 int.ov.free = TRUE)
    fit <- lavaan::lavaan(est_syntax, data = dat)
    if (!lavaan::lavInspect(fit, "converged")) {
      return(data.frame(beta_x_est = NA, beta_y_est = NA,
                        omega_xy_est = NA, omega_yx_est = NA, converged = FALSE))
    }
    pt <- lavaan::parameterEstimates(fit)
    data.frame(
      beta_x_est  = pt$est[pt$label == "beta_x"],
      beta_y_est  = pt$est[pt$label == "beta_y"],
      omega_xy_est = pt$est[pt$label == "omega_xy"],
      omega_yx_est = pt$est[pt$label == "omega_yx"],
      converged = TRUE
    )
  }, error = function(e) {
    data.frame(beta_x_est = NA, beta_y_est = NA,
               omega_xy_est = NA, omega_yx_est = NA, converged = FALSE)
  })
}


# ============================================================================
# Helper Functions: Tab 2 (Nickell Bias)
# ============================================================================

sim_riclpm_first_principles <- function(waves, sample_size,
                                         beta_x, beta_y, omega_xy, omega_yx,
                                         var_p, var_q, cov_pq,
                                         var_BX, var_BY, cov_BXBY) {
  trait_sigma <- matrix(c(var_BX, cov_BXBY, cov_BXBY, var_BY), nrow = 2)
  traits <- MASS::mvrnorm(n = sample_size, mu = c(0, 0), Sigma = trait_sigma)
  eta_x <- traits[, 1]
  eta_y <- traits[, 2]

  innov_sigma <- matrix(c(var_p, cov_pq, cov_pq, var_q), nrow = 2)
  wp_x <- matrix(NA, nrow = sample_size, ncol = waves)
  wp_y <- matrix(NA, nrow = sample_size, ncol = waves)

  w1 <- MASS::mvrnorm(n = sample_size, mu = c(0, 0), Sigma = innov_sigma)
  wp_x[, 1] <- w1[, 1]
  wp_y[, 1] <- w1[, 2]

  for (t in 2:waves) {
    innov <- MASS::mvrnorm(n = sample_size, mu = c(0, 0), Sigma = innov_sigma)
    wp_x[, t] <- beta_x * wp_x[, t - 1] + omega_yx * wp_y[, t - 1] + innov[, 1]
    wp_y[, t] <- beta_y * wp_y[, t - 1] + omega_xy * wp_x[, t - 1] + innov[, 2]
  }

  obs_x <- wp_x + eta_x
  obs_y <- wp_y + eta_y

  dat <- data.frame(obs_x, obs_y)
  names(dat) <- c(paste0("x", 1:waves), paste0("y", 1:waves))
  return(dat)
}

demean_and_reshape <- function(dat_wide, waves) {
  n <- nrow(dat_wide)
  x_cols <- paste0("x", 1:waves)
  y_cols <- paste0("y", 1:waves)

  x_mat <- as.matrix(dat_wide[, x_cols])
  y_mat <- as.matrix(dat_wide[, y_cols])

  x_mean <- rowMeans(x_mat)
  y_mean <- rowMeans(y_mat)

  x_within <- x_mat - x_mean
  y_within <- y_mat - y_mean

  records <- vector("list", n)
  for (i in seq_len(n)) {
    for (t in 2:waves) {
      records[[length(records) + 1]] <- data.frame(
        id = i, wave = t,
        x = x_mat[i, t], y = y_mat[i, t],
        within_x = x_within[i, t], within_y = y_within[i, t],
        xlag = x_mat[i, t - 1], ylag = y_mat[i, t - 1],
        within_x_lag = x_within[i, t - 1], within_y_lag = y_within[i, t - 1]
      )
    }
  }
  do.call(rbind, records)
}

run_one_nickell_trial <- function(waves, sample_size,
                                   beta_x, beta_y, omega_xy, omega_yx,
                                   var_p, var_q, cov_pq,
                                   var_BX, var_BY, cov_BXBY,
                                   dgp_method = "first_principles") {
  tryCatch({
    if (dgp_method == "first_principles") {
      dat_wide <- sim_riclpm_first_principles(
        waves = waves, sample_size = sample_size,
        beta_x = beta_x, beta_y = beta_y,
        omega_xy = omega_xy, omega_yx = omega_yx,
        var_p = var_p, var_q = var_q, cov_pq = cov_pq,
        var_BX = var_BX, var_BY = var_BY, cov_BXBY = cov_BXBY
      )
    } else {
      dgp_syntax <- build_riclpm_dgp(
        waves = waves, beta_x = beta_x, beta_y = beta_y,
        omega_xy = omega_xy, omega_yx = omega_yx,
        var_p = var_p, var_q = var_q, cov_pq = cov_pq,
        var_BX = var_BX, var_BY = var_BY, cov_BXBY = cov_BXBY
      )
      dat_wide <- lavaan::simulateData(dgp_syntax, sample.nobs = sample_size,
                                        int.ov.free = TRUE)
    }

    dat_long <- demean_and_reshape(dat_wide, waves)

    ols_y <- lm(y ~ xlag + ylag, data = dat_long)
    ols_x <- lm(x ~ xlag + ylag, data = dat_long)

    fe_y <- lm(within_y ~ within_x_lag + within_y_lag, data = dat_long)
    fe_x <- lm(within_x ~ within_x_lag + within_y_lag, data = dat_long)

    data.frame(
      ols_ar_x  = coef(ols_x)[["xlag"]],
      ols_ar_y  = coef(ols_y)[["ylag"]],
      ols_cl_xy = coef(ols_y)[["xlag"]],
      ols_cl_yx = coef(ols_x)[["ylag"]],
      fe_ar_x  = coef(fe_x)[["within_x_lag"]],
      fe_ar_y  = coef(fe_y)[["within_y_lag"]],
      fe_cl_xy = coef(fe_y)[["within_x_lag"]],
      fe_cl_yx = coef(fe_x)[["within_y_lag"]],
      converged = TRUE, dgp_method = dgp_method
    )
  }, error = function(e) {
    data.frame(
      ols_ar_x = NA, ols_ar_y = NA, ols_cl_xy = NA, ols_cl_yx = NA,
      fe_ar_x = NA, fe_ar_y = NA, fe_cl_xy = NA, fe_cl_yx = NA,
      converged = FALSE, dgp_method = dgp_method
    )
  })
}


# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  withMathJax(),

  titlePanel(
    div(
      h2("Cross-Lagged Panel Model: Bias Demonstrations"),
      p(style = "font-size: 14px; color: grey; margin-top: -5px;",
        "Comparing \\(\\hat{\\beta}_{xx}\\), \\(\\hat{\\beta}_{yy}\\),
         \\(\\hat{\\omega}_{xy}\\), \\(\\hat{\\omega}_{yx}\\)
         to their true DGP values")
    )
  ),

  tabsetPanel(
    id = "main_tabs",

    # ==== TAB 1: Trait Bias ====
    tabPanel("Trait Bias (CLPM)",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("True DGP Parameters"),
          helpText("Data generated from RI-CLPM, estimated with CLPM."),
          hr(),
          h5("Structural Parameters"),
          fluidRow(
            column(6, numericInput("t1_beta_x",
              label = HTML("\\(\\beta_{xx}\\) AR(X)"),
              value = 0.2, min = -0.9, max = 0.9, step = 0.05)),
            column(6, numericInput("t1_beta_y",
              label = HTML("\\(\\beta_{yy}\\) AR(Y)"),
              value = 0.2, min = -0.9, max = 0.9, step = 0.05))
          ),
          fluidRow(
            column(6, numericInput("t1_omega_xy",
              label = HTML("\\(\\omega_{xy}\\) CL(X\\(\\to\\)Y)"),
              value = 0.0, min = -0.5, max = 0.5, step = 0.05)),
            column(6, numericInput("t1_omega_yx",
              label = HTML("\\(\\omega_{yx}\\) CL(Y\\(\\to\\)X)"),
              value = 0.0, min = -0.5, max = 0.5, step = 0.05))
          ),
          hr(),
          h5("Trait (Between-Person) Variance"),
          sliderInput("t1_var_BX", HTML("\\(\\text{Var}(\\eta_x)\\)"),
                      min = 0, max = 3, value = 0.8, step = 0.1),
          sliderInput("t1_var_BY", HTML("\\(\\text{Var}(\\eta_y)\\)"),
                      min = 0, max = 3, value = 0.8, step = 0.1),
          numericInput("t1_cov_BXBY", HTML("\\(\\text{Cov}(\\eta_x, \\eta_y)\\)"),
                       value = 0.3, min = -1, max = 2, step = 0.1),
          hr(),
          h5("State (Within-Person) Variance"),
          fluidRow(
            column(6, numericInput("t1_var_p", HTML("\\(\\sigma^2_p\\)"),
                                   value = 1, min = 0.1, max = 3, step = 0.1)),
            column(6, numericInput("t1_var_q", HTML("\\(\\sigma^2_q\\)"),
                                   value = 1, min = 0.1, max = 3, step = 0.1))
          ),
          numericInput("t1_cov_pq", HTML("\\(\\text{Cov}(p,q)\\)"),
                       value = 0.1, min = -1, max = 1, step = 0.05),
          hr(),
          h5("Simulation Settings"),
          fluidRow(
            column(6, numericInput("t1_waves", "Waves", value = 5, min = 3, max = 10, step = 1)),
            column(6, numericInput("t1_sample_size", "N", value = 1000, min = 100, max = 10000, step = 100))
          ),
          sliderInput("t1_trials", "Monte Carlo Trials", min = 10, max = 500, value = 100, step = 10),
          hr(),

          h5("Plot Aesthetics"),
          fluidRow(
            column(6, numericInput("t1_bins", "Histogram bins", value = 40, min = 10, max = 100, step = 5)),
            column(6, selectInput("t1_hist_color", "Histogram fill",
              choices = c("steelblue", "coral", "seagreen", "mediumpurple",
                          "goldenrod", "slategrey"),
              selected = "steelblue"))
          ),
          fluidRow(
            column(6, selectInput("t1_bias_pos_color", "Bias bar (+)",
              choices = c("tomato", "coral", "firebrick", "darkorange"),
              selected = "tomato")),
            column(6, selectInput("t1_bias_neg_color", "Bias bar (-)",
              choices = c("steelblue", "dodgerblue", "slategrey", "seagreen"),
              selected = "steelblue"))
          ),
          checkboxInput("t1_free_y", "Free y-axis across facets", value = TRUE),

          hr(),
          actionButton("t1_run_btn", "Run Simulation", class = "btn-primary btn-lg btn-block"),
          br(), br(),
          textOutput("t1_progress_text"),
          hr(),
          h5("Download Results"),
          downloadButton("t1_download_csv", "Download CSV"),
          downloadButton("t1_download_rdata", "Download .RData"),
          br(), br(),
          h5("Download Plots"),
          fluidRow(
            column(4, downloadButton("t1_download_dist_pdf", "Dist. (PDF)")),
            column(4, downloadButton("t1_download_dist_png", "Dist. (PNG)")),
            column(4, downloadButton("t1_download_bias_pdf", "Bias (PDF)"))
          )
        ),
        mainPanel(
          width = 9,
          wellPanel(
            h4("How to read these results"),
            p("DGP: RI-CLPM (trait + state). Estimator: CLPM (ignores traits).",
              "As trait variance increases, CLPM estimates become biased.",
              strong("Dashed lines"), "= true values.", strong("Distributions"),
              "= CLPM estimates across MC trials.")
          ),
          plotOutput("t1_sampling_dist_plot", height = "500px"),
          br(),
          plotOutput("t1_bias_summary_plot", height = "300px"),
          br(),
          fluidRow(
            column(6, tableOutput("t1_bias_table")),
            column(6, tableOutput("t1_convergence_table"))
          )
        )
      )
    ),

    # ==== TAB 2: Nickell Bias ====
    tabPanel("Nickell Bias (FE)",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Nickell Bias Demonstration"),
          helpText("The fixed effects (demeaned OLS) estimator removes trait bias
                    but introduces Nickell bias: the demeaned lag is mechanically
                    correlated with the demeaned error. This bias is \\(O(1/T)\\) and
                    shrinks as waves increase."),
          hr(),

          h5("Structural Parameters"),
          fluidRow(
            column(6, numericInput("t2_beta_x",
              label = HTML("\\(\\beta_{xx}\\) AR(X)"),
              value = 0.3, min = -0.9, max = 0.9, step = 0.05)),
            column(6, numericInput("t2_beta_y",
              label = HTML("\\(\\beta_{yy}\\) AR(Y)"),
              value = 0.3, min = -0.9, max = 0.9, step = 0.05))
          ),
          fluidRow(
            column(6, numericInput("t2_omega_xy",
              label = HTML("\\(\\omega_{xy}\\) CL(X\\(\\to\\)Y)"),
              value = 0.0, min = -0.5, max = 0.5, step = 0.05)),
            column(6, numericInput("t2_omega_yx",
              label = HTML("\\(\\omega_{yx}\\) CL(Y\\(\\to\\)X)"),
              value = 0.0, min = -0.5, max = 0.5, step = 0.05))
          ),

          hr(),
          h5("Trait Variance"),
          sliderInput("t2_var_BX", HTML("\\(\\text{Var}(\\eta_x)\\)"),
                      min = 0, max = 3, value = 1.0, step = 0.1),
          sliderInput("t2_var_BY", HTML("\\(\\text{Var}(\\eta_y)\\)"),
                      min = 0, max = 3, value = 1.0, step = 0.1),
          numericInput("t2_cov_BXBY", HTML("\\(\\text{Cov}(\\eta_x, \\eta_y)\\)"),
                       value = 0.3, min = -1, max = 2, step = 0.1),

          hr(),
          h5("State Variance"),
          fluidRow(
            column(6, numericInput("t2_var_p", HTML("\\(\\sigma^2_p\\)"),
                                   value = 1, min = 0.1, max = 3, step = 0.1)),
            column(6, numericInput("t2_var_q", HTML("\\(\\sigma^2_q\\)"),
                                   value = 1, min = 0.1, max = 3, step = 0.1))
          ),
          numericInput("t2_cov_pq", HTML("\\(\\text{Cov}(p,q)\\)"),
                       value = 0.1, min = -1, max = 1, step = 0.05),

          hr(),
          h5("Simulation Settings"),
          helpText("Select multiple wave counts to see how Nickell bias
                    attenuates with longer panels."),
          checkboxGroupInput("t2_waves_set", "Waves to compare:",
            choices = c("3" = 3, "4" = 4, "5" = 5, "7" = 7,
                        "10" = 10, "15" = 15, "20" = 20),
            selected = c(3, 5, 10, 20),
            inline = TRUE
          ),
          numericInput("t2_sample_size", "N per trial", value = 1000, min = 100, max = 10000, step = 100),
          sliderInput("t2_trials", "MC Trials per wave count", min = 10, max = 300, value = 50, step = 10),

          hr(),
          h5("Plot Aesthetics"),
          fluidRow(
            column(6, selectInput("t2_fe_color", "FE fill",
              choices = c("steelblue", "coral", "seagreen", "mediumpurple",
                          "goldenrod", "slategrey"),
              selected = "steelblue")),
            column(6, selectInput("t2_ols_color", "OLS fill",
              choices = c("grey60", "grey40", "khaki", "rosybrown",
                          "lightcoral", "thistle"),
              selected = "grey60"))
          ),

          hr(),
          actionButton("t2_run_btn", "Run Nickell Simulation", class = "btn-warning btn-lg btn-block"),
          br(), br(),
          textOutput("t2_progress_text"),
          hr(),
          h5("Download Results"),
          downloadButton("t2_download_csv", "Download CSV"),
          downloadButton("t2_download_rdata", "Download .RData"),
          br(), br(),
          h5("Download Plots"),
          fluidRow(
            column(4, downloadButton("t2_download_nickell_pdf", "Nickell (PDF)")),
            column(4, downloadButton("t2_download_nickell_png", "Nickell (PNG)")),
            column(4, downloadButton("t2_download_biasT_pdf", "Bias(T) (PDF)"))
          )
        ),

        mainPanel(
          width = 9,

          tabsetPanel(
            id = "nickell_subtabs",

            tabPanel("Nickell Bias Results",
              wellPanel(
                h4("Nickell Bias: \\(O(1/T)\\)"),
                p("The FE estimator removes trait bias (compare OLS vs FE), but
                  introduces a downward bias on autoregressive parameters that
                  shrinks as \\(T\\) increases. With short panels (\\(T=3,4\\)), the
                  \\(\\hat{\\beta}\\) estimates are substantially attenuated."),
                p(strong("Red dashed lines"), "= true parameter values.",
                  strong("Blue"), "= FE estimates.", strong("Grey"), "= Pooled OLS estimates.")
              ),
              plotOutput("t2_nickell_plot", height = "550px"),
              br(),
              plotOutput("t2_bias_by_T_plot", height = "350px"),
              br(),
              tableOutput("t2_bias_table")
            ),

            tabPanel("DGP Comparison: First Principles vs lavaan",
              wellPanel(
                h4("Do the two DGPs produce equivalent data?"),
                p("This tab generates one dataset from each method -- first
                  principles (direct MVN sampling + AR/CL recursion) and lavaan
                  (simulateData with RI-CLPM syntax) -- and compares the
                  implied covariance structures."),
                p("The comparison runs a single large-N draw from each and
                  displays the covariance matrices, marginal variances, and
                  person-mean centered lag correlations side by side.")
              ),
              actionButton("t2_compare_btn", "Generate Comparison Data", class = "btn-info btn-block"),
              br(), br(),
              fluidRow(
                column(6,
                  h5("First Principles DGP"),
                  verbatimTextOutput("t2_fp_summary")
                ),
                column(6,
                  h5("lavaan DGP"),
                  verbatimTextOutput("t2_lav_summary")
                )
              ),
              hr(),
              h5("Covariance Matrix Comparison (difference: First Principles - lavaan)"),
              verbatimTextOutput("t2_cov_diff"),
              hr(),
              h5("FE Regression Comparison (single large-N draw)"),
              fluidRow(
                column(6,
                  h6("First Principles -> FE"),
                  verbatimTextOutput("t2_fe_fp")
                ),
                column(6,
                  h6("lavaan -> FE"),
                  verbatimTextOutput("t2_fe_lav")
                )
              ),
              hr(),
              plotOutput("t2_compare_plot", height = "400px")
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

  # ---- Tab 1: Trait Bias ----
  t1_results <- reactiveVal(NULL)

  observeEvent(input$t1_run_btn, {
    n_trials <- input$t1_trials
    withProgress(message = "Running Trait Bias MC", value = 0, {
      results_list <- vector("list", n_trials)
      for (i in seq_len(n_trials)) {
        results_list[[i]] <- run_one_trial_clpm(
          waves = input$t1_waves, sample_size = input$t1_sample_size,
          beta_x = input$t1_beta_x, beta_y = input$t1_beta_y,
          omega_xy = input$t1_omega_xy, omega_yx = input$t1_omega_yx,
          var_p = input$t1_var_p, var_q = input$t1_var_q, cov_pq = input$t1_cov_pq,
          var_BX = input$t1_var_BX, var_BY = input$t1_var_BY, cov_BXBY = input$t1_cov_BXBY
        )
        incProgress(1 / n_trials, detail = paste("Trial", i, "of", n_trials))
      }
      res <- bind_rows(results_list)
      res$trial <- seq_len(nrow(res))
      t1_results(res)
    })
  })

  output$t1_progress_text <- renderText({
    res <- t1_results()
    if (is.null(res)) return("Press 'Run Simulation' to start.")
    paste0("Completed: ", nrow(res), " trials (", sum(res$converged, na.rm = TRUE), " converged)")
  })

  # Reactive plot builders (reused for display + download)
  t1_dist_plot <- reactive({
    res <- t1_results(); req(res)
    res_conv <- res %>% filter(converged); req(nrow(res_conv) > 0)

    true_vals <- data.frame(
      parameter = c(
        "hat(beta)[xx]", "hat(beta)[yy]",
        "hat(omega)[xy]", "hat(omega)[yx]"
      ),
      true_value = c(input$t1_beta_x, input$t1_beta_y,
                     input$t1_omega_xy, input$t1_omega_yx)
    )
    long_df <- res_conv %>%
      dplyr::select(trial, beta_x_est, beta_y_est, omega_xy_est, omega_yx_est) %>%
      pivot_longer(-trial, names_to = "parameter", values_to = "estimate") %>%
      mutate(parameter = recode(parameter,
        "beta_x_est"   = "hat(beta)[xx]",
        "beta_y_est"   = "hat(beta)[yy]",
        "omega_xy_est" = "hat(omega)[xy]",
        "omega_yx_est" = "hat(omega)[yx]"))

    scales_arg <- if (input$t1_free_y) "free" else "free_x"

    ggplot(long_df, aes(x = estimate)) +
      geom_histogram(aes(y = after_stat(density)), bins = input$t1_bins,
                     fill = input$t1_hist_color, alpha = 0.6, color = "white") +
      geom_density(color = "grey20", linewidth = 0.8) +
      geom_vline(data = true_vals, aes(xintercept = true_value),
                 linetype = "dashed", color = "red", linewidth = 1) +
      facet_wrap(~parameter, scales = scales_arg, ncol = 2,
                 labeller = label_parsed) +
      labs(
        title = bquote("CLPM Sampling Distributions  |  Trait Var:"
                        ~ eta[x] == .(input$t1_var_BX) ~ ","
                        ~ eta[y] == .(input$t1_var_BY)),
        subtitle = paste0("DGP: RI-CLPM | Estimator: CLPM | ",
                          nrow(res_conv), " trials | N = ", input$t1_sample_size),
        x = "CLPM Estimate", y = "Density"
      ) +
      theme_minimal(base_size = 14) +
      theme(strip.text = element_text(face = "bold", size = 13),
            plot.title = element_text(face = "bold"),
            panel.grid.minor = element_blank())
  })

  t1_bias_plot <- reactive({
    res <- t1_results(); req(res)
    res_conv <- res %>% filter(converged); req(nrow(res_conv) > 0)

    bias_df <- data.frame(
      Parameter = factor(
        c("hat(beta)[xx]", "hat(beta)[yy]",
          "hat(omega)[xy]", "hat(omega)[yx]"),
        levels = c("hat(beta)[xx]", "hat(beta)[yy]",
                   "hat(omega)[xy]", "hat(omega)[yx]")
      ),
      True = c(input$t1_beta_x, input$t1_beta_y,
               input$t1_omega_xy, input$t1_omega_yx),
      Mean_Est = c(mean(res_conv$beta_x_est, na.rm = TRUE),
                   mean(res_conv$beta_y_est, na.rm = TRUE),
                   mean(res_conv$omega_xy_est, na.rm = TRUE),
                   mean(res_conv$omega_yx_est, na.rm = TRUE))
    ) %>% mutate(Bias = Mean_Est - True)

    ggplot(bias_df, aes(x = Parameter, y = Bias, fill = Bias > 0)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_hline(yintercept = 0, linewidth = 0.8) +
      scale_x_discrete(labels = parse(text = levels(bias_df$Parameter))) +
      scale_fill_manual(values = c("TRUE" = input$t1_bias_pos_color,
                                   "FALSE" = input$t1_bias_neg_color)) +
      labs(
        title = bquote("Bias in CLPM Estimates  " ~ (bar(hat(theta)) - theta)),
        x = "", y = "Bias"
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$t1_sampling_dist_plot <- renderPlot({ t1_dist_plot() })
  output$t1_bias_summary_plot  <- renderPlot({ t1_bias_plot() })

  output$t1_bias_table <- renderTable({
    res <- t1_results(); req(res)
    res_conv <- res %>% filter(converged); req(nrow(res_conv) > 0)
    true <- c(input$t1_beta_x, input$t1_beta_y, input$t1_omega_xy, input$t1_omega_yx)
    ests <- c(mean(res_conv$beta_x_est, na.rm=T), mean(res_conv$beta_y_est, na.rm=T),
              mean(res_conv$omega_xy_est, na.rm=T), mean(res_conv$omega_yx_est, na.rm=T))
    sds <- c(sd(res_conv$beta_x_est,na.rm=T), sd(res_conv$beta_y_est,na.rm=T),
             sd(res_conv$omega_xy_est,na.rm=T), sd(res_conv$omega_yx_est,na.rm=T))
    data.frame(
      Parameter = c("beta_xx", "beta_yy", "omega_xy", "omega_yx"),
      `True` = true, `Mean Est.` = round(ests, 4),
      `Bias` = round(ests - true, 4), `SD` = round(sds, 4), check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  output$t1_convergence_table <- renderTable({
    res <- t1_results(); req(res)
    data.frame(
      Metric = c("Total Trials", "Converged", "ICC(X)", "ICC(Y)"),
      Value = c(nrow(res), sum(res$converged, na.rm=T),
                round(input$t1_var_BX / (input$t1_var_BX + input$t1_var_p), 3),
                round(input$t1_var_BY / (input$t1_var_BY + input$t1_var_q), 3))
    )
  }, striped = TRUE, bordered = TRUE)

  # Tab 1 downloads
  output$t1_download_csv <- downloadHandler(
    filename = function() {
      paste0("trait_bias_mc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      res <- t1_results(); req(res)
      write.csv(res, file, row.names = FALSE)
    }
  )

  output$t1_download_rdata <- downloadHandler(
    filename = function() {
      paste0("trait_bias_mc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
    },
    content = function(file) {
      res <- t1_results(); req(res)
      trait_bias_results <- res
      save(trait_bias_results, file = file)
    }
  )

  # Tab 1 plot downloads
  output$t1_download_dist_pdf <- downloadHandler(
    filename = function() paste0("clpm_distributions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) { ggsave(file, plot = t1_dist_plot(), width = 10, height = 6, device = "pdf") }
  )
  output$t1_download_dist_png <- downloadHandler(
    filename = function() paste0("clpm_distributions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) { ggsave(file, plot = t1_dist_plot(), width = 10, height = 6, dpi = 300, device = "png") }
  )
  output$t1_download_bias_pdf <- downloadHandler(
    filename = function() paste0("clpm_bias_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) { ggsave(file, plot = t1_bias_plot(), width = 9, height = 4, device = "pdf") }
  )


  # ---- Tab 2: Nickell Bias ----
  t2_results <- reactiveVal(NULL)

  observeEvent(input$t2_run_btn, {
    waves_set <- as.integer(input$t2_waves_set)
    req(length(waves_set) > 0)

    n_trials <- input$t2_trials
    total_runs <- length(waves_set) * n_trials * 2
    completed <- 0

    withProgress(message = "Running Nickell Bias MC", value = 0, {
      all_results <- list()

      for (w in waves_set) {
        for (method in c("first_principles", "lavaan")) {
          for (i in seq_len(n_trials)) {
            res_row <- run_one_nickell_trial(
              waves = w, sample_size = input$t2_sample_size,
              beta_x = input$t2_beta_x, beta_y = input$t2_beta_y,
              omega_xy = input$t2_omega_xy, omega_yx = input$t2_omega_yx,
              var_p = input$t2_var_p, var_q = input$t2_var_q, cov_pq = input$t2_cov_pq,
              var_BX = input$t2_var_BX, var_BY = input$t2_var_BY, cov_BXBY = input$t2_cov_BXBY,
              dgp_method = method
            )
            res_row$waves <- w
            res_row$trial <- i
            all_results[[length(all_results) + 1]] <- res_row

            completed <- completed + 1
            incProgress(1 / total_runs,
                        detail = paste0("T=", w, " ", method, " trial ", i, "/", n_trials))
          }
        }
      }

      t2_results(bind_rows(all_results))
    })
  })

  output$t2_progress_text <- renderText({
    res <- t2_results()
    if (is.null(res)) return("Press 'Run Nickell Simulation' to start.")
    n_conv <- sum(res$converged, na.rm = TRUE)
    paste0("Completed: ", nrow(res), " runs (", n_conv, " converged)")
  })

  # Reactive plot builders (reused for display + download)
  t2_nickell_plot_obj <- reactive({
    res <- t2_results(); req(res)
    res_conv <- res %>% filter(converged); req(nrow(res_conv) > 0)

    fe_long <- res_conv %>%
      dplyr::select(waves, dgp_method, trial, fe_ar_x, fe_ar_y, fe_cl_xy, fe_cl_yx) %>%
      pivot_longer(cols = starts_with("fe_"), names_to = "parameter", values_to = "estimate") %>%
      mutate(
        parameter = recode(parameter,
          "fe_ar_x" = "FE: AR(X)", "fe_ar_y" = "FE: AR(Y)",
          "fe_cl_xy" = "FE: CL(X->Y)", "fe_cl_yx" = "FE: CL(Y->X)"),
        estimator = "FE (demeaned)"
      )

    ols_long <- res_conv %>%
      dplyr::select(waves, dgp_method, trial, ols_ar_x, ols_ar_y, ols_cl_xy, ols_cl_yx) %>%
      pivot_longer(cols = starts_with("ols_"), names_to = "parameter", values_to = "estimate") %>%
      mutate(
        parameter = recode(parameter,
          "ols_ar_x" = "FE: AR(X)", "ols_ar_y" = "FE: AR(Y)",
          "ols_cl_xy" = "FE: CL(X->Y)", "ols_cl_yx" = "FE: CL(Y->X)"),
        estimator = "Pooled OLS"
      )

    plot_df <- bind_rows(fe_long, ols_long)

    true_vals <- data.frame(
      parameter = c("FE: AR(X)", "FE: AR(Y)", "FE: CL(X->Y)", "FE: CL(Y->X)"),
      true_value = c(input$t2_beta_x, input$t2_beta_y, input$t2_omega_xy, input$t2_omega_yx)
    )

    fe_col <- input$t2_fe_color
    ols_col <- input$t2_ols_color
    fe_line <- adjustcolor(fe_col, red.f = 0.7, green.f = 0.7, blue.f = 0.7)
    ols_line <- adjustcolor(ols_col, red.f = 0.7, green.f = 0.7, blue.f = 0.7)

    ggplot(plot_df, aes(x = estimate, fill = estimator, color = estimator)) +
      geom_density(alpha = 0.3, linewidth = 0.6) +
      geom_vline(data = true_vals, aes(xintercept = true_value),
                 linetype = "dashed", color = "red", linewidth = 0.9) +
      facet_grid(parameter ~ paste0("T = ", waves), scales = "free_y") +
      scale_fill_manual(values = c("FE (demeaned)" = fe_col, "Pooled OLS" = ols_col)) +
      scale_color_manual(values = c("FE (demeaned)" = fe_line, "Pooled OLS" = ols_line)) +
      labs(
        title = bquote("Nickell Bias: FE vs Pooled OLS by " ~ T),
        subtitle = paste0("Both DGP methods | N = ", input$t2_sample_size,
                          " | Var(", "\u03B7", "x) = ", input$t2_var_BX,
                          "  Var(", "\u03B7", "y) = ", input$t2_var_BY),
        x = "Estimate", y = "Density", fill = "Estimator", color = "Estimator"
      ) +
      theme_minimal(base_size = 13) +
      theme(strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold"),
            legend.position = "bottom")
  })

  t2_biasT_plot_obj <- reactive({
    res <- t2_results(); req(res)
    res_conv <- res %>% filter(converged); req(nrow(res_conv) > 0)

    bias_summary <- res_conv %>%
      group_by(waves, dgp_method) %>%
      summarise(
        `FE: AR(X) bias` = mean(fe_ar_x, na.rm = TRUE) - input$t2_beta_x,
        `FE: AR(Y) bias` = mean(fe_ar_y, na.rm = TRUE) - input$t2_beta_y,
        `FE: CL(X->Y) bias` = mean(fe_cl_xy, na.rm = TRUE) - input$t2_omega_xy,
        `FE: CL(Y->X) bias` = mean(fe_cl_yx, na.rm = TRUE) - input$t2_omega_yx,
        `OLS: AR(X) bias` = mean(ols_ar_x, na.rm = TRUE) - input$t2_beta_x,
        `OLS: AR(Y) bias` = mean(ols_ar_y, na.rm = TRUE) - input$t2_beta_y,
        .groups = "drop"
      ) %>%
      pivot_longer(cols = -c(waves, dgp_method), names_to = "parameter", values_to = "bias")

    ggplot(bias_summary, aes(x = waves, y = bias, color = parameter,
                              linetype = dgp_method, shape = dgp_method)) +
      geom_hline(yintercept = 0, linewidth = 0.6, color = "grey50") +
      geom_line(linewidth = 0.8) +
      geom_point(size = 3) +
      scale_linetype_manual(values = c("first_principles" = "solid", "lavaan" = "dashed")) +
      labs(
        title = bquote("Bias as a Function of " ~ T ~ " (Number of Waves)"),
        subtitle = "Nickell bias (FE) attenuates as T grows; OLS trait bias persists",
        x = "Number of Waves (T)", y = "Bias (Mean Est. - True)",
        color = "Parameter", linetype = "DGP Method", shape = "DGP Method"
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  })

  output$t2_nickell_plot  <- renderPlot({ t2_nickell_plot_obj() })
  output$t2_bias_by_T_plot <- renderPlot({ t2_biasT_plot_obj() })

  output$t2_bias_table <- renderTable({
    res <- t2_results(); req(res)
    res_conv <- res %>% filter(converged); req(nrow(res_conv) > 0)

    res_conv %>%
      group_by(Waves = waves, `DGP Method` = dgp_method) %>%
      summarise(
        `FE AR(X)` = round(mean(fe_ar_x, na.rm = TRUE), 4),
        `FE AR(Y)` = round(mean(fe_ar_y, na.rm = TRUE), 4),
        `FE CL(X->Y)` = round(mean(fe_cl_xy, na.rm = TRUE), 4),
        `FE CL(Y->X)` = round(mean(fe_cl_yx, na.rm = TRUE), 4),
        `OLS AR(X)` = round(mean(ols_ar_x, na.rm = TRUE), 4),
        `OLS AR(Y)` = round(mean(ols_ar_y, na.rm = TRUE), 4),
        N = n(),
        .groups = "drop"
      )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  # Tab 2 downloads
  output$t2_download_csv <- downloadHandler(
    filename = function() {
      paste0("nickell_bias_mc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      res <- t2_results(); req(res)
      write.csv(res, file, row.names = FALSE)
    }
  )

  output$t2_download_rdata <- downloadHandler(
    filename = function() {
      paste0("nickell_bias_mc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
    },
    content = function(file) {
      res <- t2_results(); req(res)
      nickell_bias_results <- res
      save(nickell_bias_results, file = file)
    }
  )

  # Tab 2 plot downloads
  output$t2_download_nickell_pdf <- downloadHandler(
    filename = function() paste0("nickell_bias_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) { ggsave(file, plot = t2_nickell_plot_obj(), width = 12, height = 8, device = "pdf") }
  )
  output$t2_download_nickell_png <- downloadHandler(
    filename = function() paste0("nickell_bias_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) { ggsave(file, plot = t2_nickell_plot_obj(), width = 12, height = 8, dpi = 300, device = "png") }
  )
  output$t2_download_biasT_pdf <- downloadHandler(
    filename = function() paste0("bias_by_T_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) { ggsave(file, plot = t2_biasT_plot_obj(), width = 10, height = 5, device = "pdf") }
  )


  # ---- Tab 2 sub-tab: DGP Comparison ----
  t2_compare_data <- reactiveVal(NULL)

  observeEvent(input$t2_compare_btn, {
    withProgress(message = "Generating comparison data (N=5000)...", value = 0.2, {

      big_n <- 5000
      w <- as.integer(input$t2_waves_set[1])
      if (is.na(w)) w <- 5

      fp_dat <- sim_riclpm_first_principles(
        waves = w, sample_size = big_n,
        beta_x = input$t2_beta_x, beta_y = input$t2_beta_y,
        omega_xy = input$t2_omega_xy, omega_yx = input$t2_omega_yx,
        var_p = input$t2_var_p, var_q = input$t2_var_q, cov_pq = input$t2_cov_pq,
        var_BX = input$t2_var_BX, var_BY = input$t2_var_BY, cov_BXBY = input$t2_cov_BXBY
      )

      incProgress(0.4, detail = "Generating lavaan data...")

      dgp_syntax <- build_riclpm_dgp(
        waves = w, beta_x = input$t2_beta_x, beta_y = input$t2_beta_y,
        omega_xy = input$t2_omega_xy, omega_yx = input$t2_omega_yx,
        var_p = input$t2_var_p, var_q = input$t2_var_q, cov_pq = input$t2_cov_pq,
        var_BX = input$t2_var_BX, var_BY = input$t2_var_BY, cov_BXBY = input$t2_cov_BXBY
      )
      lav_dat <- lavaan::simulateData(dgp_syntax, sample.nobs = big_n, int.ov.free = TRUE)

      incProgress(0.3, detail = "Computing comparisons...")

      fp_long <- demean_and_reshape(fp_dat, w)
      lav_long <- demean_and_reshape(lav_dat, w)

      fp_fe_y <- lm(within_y ~ within_x_lag + within_y_lag, data = fp_long)
      fp_fe_x <- lm(within_x ~ within_x_lag + within_y_lag, data = fp_long)
      lav_fe_y <- lm(within_y ~ within_x_lag + within_y_lag, data = lav_long)
      lav_fe_x <- lm(within_x ~ within_x_lag + within_y_lag, data = lav_long)

      t2_compare_data(list(
        fp_dat = fp_dat, lav_dat = lav_dat,
        fp_long = fp_long, lav_long = lav_long,
        fp_fe_y = fp_fe_y, fp_fe_x = fp_fe_x,
        lav_fe_y = lav_fe_y, lav_fe_x = lav_fe_x,
        waves = w
      ))
    })
  })

  output$t2_fp_summary <- renderPrint({
    d <- t2_compare_data(); req(d)
    cat("Observed variable covariance matrix (first 6 vars):\n")
    vars_show <- intersect(names(d$fp_dat), c("x1","x2","x3","y1","y2","y3"))
    print(round(cov(d$fp_dat[, vars_show]), 3))
    cat("\nColumn means:\n")
    print(round(colMeans(d$fp_dat[, vars_show]), 3))
    cat("\nColumn variances:\n")
    print(round(apply(d$fp_dat[, vars_show], 2, var), 3))
  })

  output$t2_lav_summary <- renderPrint({
    d <- t2_compare_data(); req(d)
    cat("Observed variable covariance matrix (first 6 vars):\n")
    vars_show <- intersect(names(d$lav_dat), c("x1","x2","x3","y1","y2","y3"))
    print(round(cov(d$lav_dat[, vars_show]), 3))
    cat("\nColumn means:\n")
    print(round(colMeans(d$lav_dat[, vars_show]), 3))
    cat("\nColumn variances:\n")
    print(round(apply(d$lav_dat[, vars_show], 2, var), 3))
  })

  output$t2_cov_diff <- renderPrint({
    d <- t2_compare_data(); req(d)
    vars_show <- intersect(names(d$fp_dat), c("x1","x2","x3","y1","y2","y3"))
    fp_cov <- cov(d$fp_dat[, vars_show])
    lav_cov <- cov(d$lav_dat[, vars_show])
    cat("Difference in covariance matrices (FP - lavaan):\n")
    cat("Values near 0 confirm equivalent DGPs (sampling noise only)\n\n")
    print(round(fp_cov - lav_cov, 3))
    cat("\nMax absolute difference:", round(max(abs(fp_cov - lav_cov)), 4), "\n")
  })

  output$t2_fe_fp <- renderPrint({
    d <- t2_compare_data(); req(d)
    cat("=== X equation: within_x ~ within_x_lag + within_y_lag ===\n")
    print(summary(d$fp_fe_x)$coefficients)
    cat("\n=== Y equation: within_y ~ within_x_lag + within_y_lag ===\n")
    print(summary(d$fp_fe_y)$coefficients)
  })

  output$t2_fe_lav <- renderPrint({
    d <- t2_compare_data(); req(d)
    cat("=== X equation: within_x ~ within_x_lag + within_y_lag ===\n")
    print(summary(d$lav_fe_x)$coefficients)
    cat("\n=== Y equation: within_y ~ within_x_lag + within_y_lag ===\n")
    print(summary(d$lav_fe_y)$coefficients)
  })

  output$t2_compare_plot <- renderPlot({
    d <- t2_compare_data(); req(d)

    fp_df <- d$fp_long %>% mutate(DGP = "First Principles")
    lav_df <- d$lav_long %>% mutate(DGP = "lavaan")
    both <- bind_rows(fp_df, lav_df)

    p1 <- ggplot(both, aes(x = within_x_lag, fill = DGP, color = DGP)) +
      geom_density(alpha = 0.3) +
      labs(title = "Within-person X lag", x = expression(tilde(x)[t-1])) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")

    p2 <- ggplot(both, aes(x = within_y_lag, fill = DGP, color = DGP)) +
      geom_density(alpha = 0.3) +
      labs(title = "Within-person Y lag", x = expression(tilde(y)[t-1])) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")

    p3 <- ggplot(both, aes(x = within_x, fill = DGP, color = DGP)) +
      geom_density(alpha = 0.3) +
      labs(title = "Within-person X (current)", x = expression(tilde(x)[t])) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")

    p4 <- ggplot(both, aes(x = within_y, fill = DGP, color = DGP)) +
      geom_density(alpha = 0.3) +
      labs(title = "Within-person Y (current)", x = expression(tilde(y)[t])) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")

    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2,
      top = grid::textGrob("Distribution Comparison: First Principles vs lavaan (demeaned variables)",
                            gp = grid::gpar(fontface = "bold", fontsize = 14)))
  })
}

shinyApp(ui = ui, server = server)
