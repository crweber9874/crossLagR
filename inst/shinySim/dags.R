## DAG DOT-source library for the crossLagR Shiny app.
## Each function returns a DiagrammeR-ready DOT string for a given waves count.
## The diagrams mirror those in the quarto book (quarto_book/01_Panel_Models.qmd,
## 02_Change.qmd, 04_Usami.qmd, 06_iscap.qmd) but were rewritten as small
## parameterized generators so the Shiny app can redraw at the user's chosen
## wave count.

## ---- shared helpers ----------------------------------------------------------

.dag_header <- function(title = "") {
  paste0(
    "digraph G {\n",
    "  graph [layout = neato, bgcolor = white, splines = true, ",
    "         labelloc=t, label = '", title, "', fontname='Helvetica', fontsize=14]\n",
    "  node [shape = box, style = filled, fillcolor = '#F8F9FA', ",
    "        fontname = 'Helvetica', fontsize = 11, color = '#343A40']\n",
    "  edge [color = '#343A40', fontname = 'Helvetica', fontsize = 10]\n"
  )
}

.dag_footer <- "}\n"

.obs_row <- function(prefix, waves, y, label_prefix = NULL) {
  ## prefix: "X" or "Y"; waves: integer; y: vertical coord; label_prefix: text
  if (is.null(label_prefix)) label_prefix <- tolower(prefix)
  paste(vapply(seq_len(waves), function(w) {
    sprintf("  %s%d [pos = '%d,%g!', label = '%s%d', width = 0.6, height = 0.3]",
            prefix, w, w - 1, y, label_prefix, w)
  }, character(1)), collapse = "\n")
}

.latent_row <- function(prefix, waves, y, label_prefix = NULL) {
  if (is.null(label_prefix)) label_prefix <- tolower(prefix)
  paste(vapply(seq_len(waves), function(w) {
    sprintf("  %s%d [pos = '%d,%g!', shape = circle, fillcolor = white, label = '%s%d', width = 0.45, height = 0.45]",
            prefix, w, w - 1, y, label_prefix, w)
  }, character(1)), collapse = "\n")
}

.ar_paths <- function(prefix, waves, label = "", color = "#1F77B4") {
  if (waves < 2) return("")
  paste(vapply(2:waves, function(w) {
    sprintf("  %s%d -> %s%d [color='%s', penwidth = 1.5, arrowsize = 0.6, label='%s']",
            prefix, w - 1, prefix, w, color, label)
  }, character(1)), collapse = "\n")
}

.cl_paths <- function(from, to, waves, label = "", color = "#D62728") {
  if (waves < 2) return("")
  paste(vapply(2:waves, function(w) {
    sprintf("  %s%d -> %s%d [color='%s', penwidth = 1.2, arrowsize = 0.6, label='%s']",
            from, w - 1, to, w, color, label)
  }, character(1)), collapse = "\n")
}

.measurement_paths <- function(latent_prefix, obs_prefix, waves) {
  paste(vapply(seq_len(waves), function(w) {
    sprintf("  %s%d -> %s%d [label = '1', arrowsize = 0.45, color='#6C757D']",
            latent_prefix, w, obs_prefix, w)
  }, character(1)), collapse = "\n")
}

.ri_loadings <- function(ri_node, latent_prefix, waves) {
  paste(vapply(seq_len(waves), function(w) {
    sprintf("  %s -> %s%d [label = '1', arrowsize = 0.4, color='#6C757D']",
            ri_node, latent_prefix, w)
  }, character(1)), collapse = "\n")
}

.within_time_cov <- function(a, b, waves, color = "#6C757D") {
  paste(vapply(seq_len(waves), function(w) {
    sprintf("  %s%d -> %s%d [dir=both, style=dashed, splines=curved, constraint=false, arrowsize=0.4, color='%s']",
            a, w, b, w, color)
  }, character(1)), collapse = "\n")
}

## ---- DAG generators (one per estimator) --------------------------------------

dag_CLPM <- function(waves = 4) {
  paste0(
    .dag_header("Cross-Lagged Panel Model (CLPM)"),
    .obs_row("X", waves,  1), "\n",
    .obs_row("Y", waves, -1), "\n",
    .ar_paths("X", waves, label = "ar_x", color = "#1F77B4"), "\n",
    .ar_paths("Y", waves, label = "ar_y", color = "#1F77B4"), "\n",
    .cl_paths("X", "Y", waves, label = "cl_xy", color = "#D62728"), "\n",
    .cl_paths("Y", "X", waves, label = "cl_yx", color = "#D62728"), "\n",
    .within_time_cov("X", "Y", waves), "\n",
    .dag_footer
  )
}

dag_RICLPM <- function(waves = 4) {
  paste0(
    .dag_header("Random-Intercept CLPM (RI-CLPM)"),
    sprintf("  Ix [pos = '%g,3!', shape = circle, fillcolor='#E8F1FA', label = 'I_x', width = 0.5, height = 0.5]\n", -1),
    sprintf("  Iy [pos = '%g,-3!', shape = circle, fillcolor='#FDECEC', label = 'I_y', width = 0.5, height = 0.5]\n", -1),
    .obs_row("X", waves,  2), "\n",
    .latent_row("q", waves,  1, label_prefix = "q"), "\n",
    .latent_row("p", waves, -1, label_prefix = "p"), "\n",
    .obs_row("Y", waves, -2), "\n",
    .ri_loadings("Ix", "X", waves), "\n",
    .ri_loadings("Iy", "Y", waves), "\n",
    .measurement_paths("q", "X", waves), "\n",
    .measurement_paths("p", "Y", waves), "\n",
    .ar_paths("q", waves, label = "ar_x", color = "#1F77B4"), "\n",
    .ar_paths("p", waves, label = "ar_y", color = "#1F77B4"), "\n",
    .cl_paths("q", "p", waves, label = "cl_xy", color = "#D62728"), "\n",
    .cl_paths("p", "q", waves, label = "cl_yx", color = "#D62728"), "\n",
    .within_time_cov("q", "p", waves), "\n",
    "  Ix -> Iy [dir=both, style=dashed, color='#6C757D', constraint=false, arrowsize=0.5, ",
    "tailport=sw, headport=nw, minlen=8]\n",
    .dag_footer
  )
}

dag_ALT <- function(waves = 4) {
  paste0(
    .dag_header("Autoregressive Latent Trajectory (ALT)"),
    "  Ix [pos = '-2,3!', shape = circle, fillcolor='#E8F1FA', label = 'I_x', width = 0.5, height = 0.5]\n",
    "  Sx [pos = '-2,2!', shape = circle, fillcolor='#E8F1FA', label = 'S_x', width = 0.5, height = 0.5]\n",
    "  Iy [pos = '-2,-2!', shape = circle, fillcolor='#FDECEC', label = 'I_y', width = 0.5, height = 0.5]\n",
    "  Sy [pos = '-2,-3!', shape = circle, fillcolor='#FDECEC', label = 'S_y', width = 0.5, height = 0.5]\n",
    .obs_row("X", waves,  1), "\n",
    .obs_row("Y", waves, -1), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  Ix -> X%d [label='1', arrowsize=0.4, color='#6C757D']\n  Sx -> X%d [label='%d', arrowsize=0.4, color='#6C757D']",
              w, w, w - 1)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  Iy -> Y%d [label='1', arrowsize=0.4, color='#6C757D']\n  Sy -> Y%d [label='%d', arrowsize=0.4, color='#6C757D']",
              w, w, w - 1)
    }, character(1)), collapse = "\n"), "\n",
    .ar_paths("X", waves, label = "ar_x", color = "#1F77B4"), "\n",
    .ar_paths("Y", waves, label = "ar_y", color = "#1F77B4"), "\n",
    .cl_paths("X", "Y", waves, label = "cl_xy", color = "#D62728"), "\n",
    .cl_paths("Y", "X", waves, label = "cl_yx", color = "#D62728"), "\n",
    "  Ix -> Iy [dir=both, style=dashed, color='#6C757D', constraint=false]\n",
    "  Sx -> Sy [dir=both, style=dashed, color='#6C757D', constraint=false]\n",
    .dag_footer
  )
}

dag_LGM <- function(waves = 4) {
  paste0(
    .dag_header("Latent Growth Model (LGM)"),
    "  Ix [pos = '-2,2.6!', shape=circle, fillcolor='#E8F1FA', label='I_x', width=0.5, height=0.5]\n",
    "  Sx [pos = '-2,1.4!', shape=circle, fillcolor='#E8F1FA', label='S_x', width=0.5, height=0.5]\n",
    "  Iy [pos = '-2,-1.4!', shape=circle, fillcolor='#FDECEC', label='I_y', width=0.5, height=0.5]\n",
    "  Sy [pos = '-2,-2.6!', shape=circle, fillcolor='#FDECEC', label='S_y', width=0.5, height=0.5]\n",
    .obs_row("X", waves,  1), "\n",
    .obs_row("Y", waves, -1), "\n",
    paste(vapply(seq_len(waves), function(w) {
      sprintf("  Ix -> X%d [label='1', arrowsize=0.4, color='#6C757D']\n  Sx -> X%d [label='%d', arrowsize=0.4, color='#6C757D']",
              w, w, w - 1)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(seq_len(waves), function(w) {
      sprintf("  Iy -> Y%d [label='1', arrowsize=0.4, color='#6C757D']\n  Sy -> Y%d [label='%d', arrowsize=0.4, color='#6C757D']",
              w, w, w - 1)
    }, character(1)), collapse = "\n"), "\n",
    "  Ix -> Iy [dir=both, style=dashed, color='#6C757D']\n",
    "  Sx -> Sy [dir=both, style=dashed, color='#6C757D']\n",
    "  Ix -> Sx [dir=both, style=dashed, color='#6C757D']\n",
    "  Iy -> Sy [dir=both, style=dashed, color='#6C757D']\n",
    .dag_footer
  )
}

dag_LCMSR <- function(waves = 4) {
  paste0(
    .dag_header("LCM-SR / GCLM"),
    "  Ix [pos = '-2,3!', shape=circle, fillcolor='#E8F1FA', label='I_x', width=0.5, height=0.5]\n",
    "  Sx [pos = '-2,2!', shape=circle, fillcolor='#E8F1FA', label='S_x', width=0.5, height=0.5]\n",
    "  Iy [pos = '-2,-2!', shape=circle, fillcolor='#FDECEC', label='I_y', width=0.5, height=0.5]\n",
    "  Sy [pos = '-2,-3!', shape=circle, fillcolor='#FDECEC', label='S_y', width=0.5, height=0.5]\n",
    .obs_row("X", waves,  2.6), "\n",
    .latent_row("q", waves,  1), "\n",
    .latent_row("p", waves, -1), "\n",
    .obs_row("Y", waves, -2.6), "\n",
    paste(vapply(seq_len(waves), function(w) {
      sprintf("  Ix -> X%d [label='1', arrowsize=0.4, color='#6C757D']\n  Sx -> X%d [label='%d', arrowsize=0.4, color='#6C757D']",
              w, w, w - 1)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(seq_len(waves), function(w) {
      sprintf("  Iy -> Y%d [label='1', arrowsize=0.4, color='#6C757D']\n  Sy -> Y%d [label='%d', arrowsize=0.4, color='#6C757D']",
              w, w, w - 1)
    }, character(1)), collapse = "\n"), "\n",
    .measurement_paths("q", "X", waves), "\n",
    .measurement_paths("p", "Y", waves), "\n",
    .ar_paths("q", waves, label = "ar_x", color = "#1F77B4"), "\n",
    .ar_paths("p", waves, label = "ar_y", color = "#1F77B4"), "\n",
    .cl_paths("q", "p", waves, label = "cl_xy", color = "#D62728"), "\n",
    .cl_paths("p", "q", waves, label = "cl_yx", color = "#D62728"), "\n",
    "  Ix -> Iy [dir=both, style=dashed, color='#6C757D']\n",
    "  Sx -> Sy [dir=both, style=dashed, color='#6C757D']\n",
    .dag_footer
  )
}

dag_LCHANGE <- function(waves = 4) {
  ## Bivariate latent change with constant-change accumulating factors
  paste0(
    .dag_header("Latent Change Score Model (LCS)"),
    .obs_row("X", waves,  2), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  dX%d [pos='%d,0.6!', shape=circle, fillcolor='#E8F1FA', label='dx%d', width=0.5, height=0.5]",
              w, w - 1, w)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  dY%d [pos='%d,-0.6!', shape=circle, fillcolor='#FDECEC', label='dy%d', width=0.5, height=0.5]",
              w, w - 1, w)
    }, character(1)), collapse = "\n"), "\n",
    .obs_row("Y", waves, -2), "\n",
    "  Ax [pos='-1.5,1.4!', shape=circle, fillcolor='#E8F1FA', label='A_x', width=0.5, height=0.5]\n",
    "  Ay [pos='-1.5,-1.4!', shape=circle, fillcolor='#FDECEC', label='A_y', width=0.5, height=0.5]\n",
    .ar_paths("X", waves, label = "1", color = "#6C757D"), "\n",
    .ar_paths("Y", waves, label = "1", color = "#6C757D"), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  dX%d -> X%d [label='1', arrowsize=0.4, color='#6C757D']\n  dY%d -> Y%d [label='1', arrowsize=0.4, color='#6C757D']",
              w, w, w, w)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  X%d -> dX%d [label='ar_x', arrowsize=0.5, color='#1F77B4']\n  Y%d -> dY%d [label='ar_y', arrowsize=0.5, color='#1F77B4']",
              w - 1, w, w - 1, w)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  Y%d -> dX%d [label='cl_yx', arrowsize=0.5, color='#D62728']\n  X%d -> dY%d [label='cl_xy', arrowsize=0.5, color='#D62728']",
              w - 1, w, w - 1, w)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  Ax -> dX%d [label='1', arrowsize=0.4, color='#6C757D']\n  Ay -> dY%d [label='1', arrowsize=0.4, color='#6C757D']",
              w, w)
    }, character(1)), collapse = "\n"), "\n",
    .dag_footer
  )
}

dag_BB <- function(waves = 4) {
  paste0(
    .dag_header("Bollen-Brand Dynamic Panel"),
    "  eta_x [pos='-2,2.8!', shape=circle, fillcolor='#FFF3CD', label='η_x', width=0.6, height=0.6]\n",
    "  eta_y [pos='-2,-2.8!', shape=circle, fillcolor='#FFF3CD', label='η_y', width=0.6, height=0.6]\n",
    .obs_row("X", waves,  1.5), "\n",
    .obs_row("Y", waves, -1.5), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  eta_x -> X%d [label='1', arrowsize=0.4, color='#856404']\n  eta_y -> Y%d [label='1', arrowsize=0.4, color='#856404']",
              w, w)
    }, character(1)), collapse = "\n"), "\n",
    .ar_paths("X", waves, label = "ar_x", color = "#1F77B4"), "\n",
    .ar_paths("Y", waves, label = "ar_y", color = "#1F77B4"), "\n",
    .cl_paths("X", "Y", waves, label = "cl_xy", color = "#D62728"), "\n",
    .cl_paths("Y", "X", waves, label = "cl_yx", color = "#D62728"), "\n",
    "  eta_x -> eta_y [dir=both, style=dashed, color='#856404']\n",
    .within_time_cov("X", "Y", 1), "\n",  ## allow x1 ↔ y1 covariance
    "  X1 -> eta_x [dir=both, style=dotted, color='#6C757D']\n",
    "  X1 -> eta_y [dir=both, style=dotted, color='#6C757D']\n",
    "  Y1 -> eta_x [dir=both, style=dotted, color='#6C757D']\n",
    "  Y1 -> eta_y [dir=both, style=dotted, color='#6C757D']\n",
    .dag_footer
  )
}

dag_TSO <- function(waves = 4) {
  paste0(
    .dag_header("Trait-State-Occasion (TSO)"),
    "  Tx [pos='-2,3!', shape=circle, fillcolor='#E8F1FA', label='T_x', width=0.55, height=0.55]\n",
    "  Ty [pos='-2,-3!', shape=circle, fillcolor='#FDECEC', label='T_y', width=0.55, height=0.55]\n",
    .obs_row("X", waves,  1), "\n",
    .obs_row("Y", waves, -1), "\n",
    .ri_loadings("Tx", "X", waves), "\n",
    .ri_loadings("Ty", "Y", waves), "\n",
    paste(vapply(seq_len(waves), function(w) {
      sprintf("  s_x%d [pos='%d,2.0!', shape=circle, fillcolor='white', label='s%d', width=0.35, height=0.35]\n  s_y%d [pos='%d,-2.0!', shape=circle, fillcolor='white', label='s%d', width=0.35, height=0.35]",
              w, w - 1, w, w, w - 1, w)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(seq_len(waves), function(w) {
      sprintf("  s_x%d -> X%d [arrowsize=0.4, color='#6C757D']\n  s_y%d -> Y%d [arrowsize=0.4, color='#6C757D']",
              w, w, w, w)
    }, character(1)), collapse = "\n"), "\n",
    .within_time_cov("s_x", "s_y", waves), "\n",
    "  Tx -> Ty [dir=both, style=dashed, color='#6C757D']\n",
    .dag_footer
  )
}

dag_CTSEM <- function(waves = 4) {
  ## CTSEM is a continuous-time SEM. Sketch the structure.
  paste0(
    .dag_header("Continuous-Time SEM (CTSEM)"),
    .obs_row("X", waves,  1), "\n",
    .obs_row("Y", waves, -1), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  X%d -> X%d [color='#1F77B4', label='exp(A·Δt)', arrowsize=0.6]\n  Y%d -> Y%d [color='#1F77B4', arrowsize=0.6]",
              w - 1, w, w - 1, w)
    }, character(1)), collapse = "\n"), "\n",
    paste(vapply(2:waves, function(w) {
      sprintf("  X%d -> Y%d [color='#D62728', arrowsize=0.6]\n  Y%d -> X%d [color='#D62728', arrowsize=0.6]",
              w - 1, w, w - 1, w)
    }, character(1)), collapse = "\n"), "\n",
    .within_time_cov("X", "Y", waves), "\n",
    .dag_footer
  )
}

dag_RI <- function(waves = 4) {
  paste0(
    .dag_header("Mixed model with random intercept (RI)"),
    "  uX [pos='-2,1.5!', shape=circle, fillcolor='#E8F1FA', label='u_x', width=0.5, height=0.5]\n",
    "  uY [pos='-2,-1.5!', shape=circle, fillcolor='#FDECEC', label='u_y', width=0.5, height=0.5]\n",
    .obs_row("X", waves,  1), "\n",
    .obs_row("Y", waves, -1), "\n",
    .ri_loadings("uX", "X", waves), "\n",
    .ri_loadings("uY", "Y", waves), "\n",
    .ar_paths("X", waves, label = "ar_x", color = "#1F77B4"), "\n",
    .ar_paths("Y", waves, label = "ar_y", color = "#1F77B4"), "\n",
    .cl_paths("X", "Y", waves, label = "cl_xy", color = "#D62728"), "\n",
    .cl_paths("Y", "X", waves, label = "cl_yx", color = "#D62728"), "\n",
    .dag_footer
  )
}

dag_OLS <- function(waves = 4) {
  paste0(
    .dag_header("OLS wide-format regression"),
    .obs_row("X", waves,  1), "\n",
    .obs_row("Y", waves, -1), "\n",
    .ar_paths("X", waves, label = "ar_x", color = "#1F77B4"), "\n",
    .ar_paths("Y", waves, label = "ar_y", color = "#1F77B4"), "\n",
    .cl_paths("X", "Y", waves, label = "cl_xy", color = "#D62728"), "\n",
    .cl_paths("Y", "X", waves, label = "cl_yx", color = "#D62728"), "\n",
    .dag_footer
  )
}

## ---- master dispatch ---------------------------------------------------------

get_dag_dot <- function(estimator, waves = 4) {
  fn <- switch(
    estimator,
    "CLPM"    = dag_CLPM,
    "RICLPM"  = dag_RICLPM,
    "ALT"     = dag_ALT,
    "LGM"     = dag_LGM,
    "LCMSR"   = dag_LCMSR,
    "LCHANGE" = dag_LCHANGE,
    "BB"      = dag_BB,
    "TSO"     = dag_TSO,
    "CTSEM"   = dag_CTSEM,
    "RI"      = dag_RI,
    "OLS"     = dag_OLS,
    NULL
  )
  if (is.null(fn)) return(NULL)
  ## Clamp waves to a reasonable range to keep diagrams legible
  fn(min(max(as.integer(waves), 2L), 8L))
}
