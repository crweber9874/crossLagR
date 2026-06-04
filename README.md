# README


## The structure of this repository is as follows:

- `README.md`: This file, which provides an overview of the repository.
- `R/`: Contains all the R scripts for the project, packaged into `crossLagR`
There are five simulation files. These files simulate datasets using `lavaan.sim` based on particular data-generating processes (DGPs).
For instance, `simCLPM` includes code to run simulations under a cross-lagged regression DGP. On the other hand, `simRICLPM3`
simulates a random-intercept cross-lagged panel with three variables, and `simRICLPM` with two. Finally, `simRICLPMl` specifies a latent-variable
simulation when more than one indicator per time point is included. All files with `sim` as a prefix are simulation files.

This directory also includes estimation functions, which are useful for comparing estimators across DGPs. These are prefixed with `estimate`.
Thus far, only `estimateCLPM`, `estimateRICLPM`, `estimateStTr` (a trait-state model), and `estimateHLM` (a hierarchical model in `rstan`) 
have been thoroughly tested.





- `src/`: Contains the source code for the project.
- `tests/`: Contains unit tests for the project.
- `examples/`: Contains example scripts demonstrating how to use the project.
- `LICENSE`: The license file for the project.


## Project Overview
This repository contains the `crossLagR` project, which provides simulation and estimation tools for longitudinal panel models (for example, CLPM, RI-CLPM, and related trait-state and hierarchical variants). It includes source code, documentation, tests, and examples to support applied use and replication.

# Description
# Installation
# Usage
# Contributing

## Technologies
- R
- Shiny

## Use of AI

AI has very much changed my typical workflow. While I would have typically considered
myself fairly well versed in R, *Shiny*, and a bit of Python, AI has made it far easier
to transition from hours of tedious coding to orchestrating and monitoring these tasks,
and producing accurate and interpretable output that might be used for scientific publication.
The *Claude Code CLI* is used to quickly structure code, which must necessarily be monitored --
and has been in this circumstance -- as errors, some of which can cascade across one's analytic
workflow, do occur. AI assistance was used for:

- **Boilerplate code and scaffolding framework** (roxygen2 skeletons, testthat templates,
  package structure). Packages require functions built upon functions,
  uniform labeling, and documentation placeholders, much of which is
  mechanical and well-suited to AI drafting.
- **Expanding base code in functions.** Many of the functions in this package
  are `lavaan` models. `lavaan` model code consists of a long string of
  `lavaan` commands, in the tradition of a structural equation model
  as specified in packages like *Mplus* and *LISREL*. We built the
  base models and used AI to expand upon these models to include more
  complex parameterizations and simulations (e.g., extending the bivariate
  CLPM to the RI-CLPM, latent-variable variants, and three-variable
  generalizations).
- **Generating helper functions.** This package builds upon other packages,
  such as `lavaan`, `brms`, `blavaan`, `lme4`, and the `tidyverse`.
  Helper functions are used to link input/output across the
  different estimators in the package, as well as repetitive
  data wrangling exercises, like when one applies the same `dplyr::recode()`
  logic across many variables, and pivoting data.
- **Generating Shiny UI/server boilerplate.** The Shiny app draws in the
  code developed in the book chapter *Simulation*, linking these elements
  in a common interface, using reactive programming, to allow the user
  to specify different simulation protocols. AI was used to develop and
  refine this architecture.

### How AI was *not* used

- AI was not used to interpret or describe statistical results reported
  in the accompanying user guide. 
- AI was not used to make modeling decisions (priors, identification
  constraints, parameter labeling conventions).
- AI was not used to validate correctness against published benchmarks (e.g., Hamaker et al. 2015; Usami et al. 2019).
- AI was not used to write the text of the manuscript,
  vignettes, or scientific interpretation.
- AI was not used to alter any of the original `.rda` data objects shipped with
  the package.

All analytical decisions, model specifications, parameterizations, and scientific conclusions are those
of the authors. All AI-generated code was reviewed by the authors, tested in the functions
`testthat`, and validated against known results before inclusion.
Responsibility for correctness rests with the authors.
