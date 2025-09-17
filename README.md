# README


## The structure of this repository is as follows:

- `README.md`: This file, which provides an overview of the repository.
- `R/`: Contains all the R scripts for the project, packaged into `crossLagR`
There are five simulation files. These files simply simulate a dataset using `lavaan.sim` based on a particular data generating process (DGP).
When `simCLPM` for instance includes the code to run a simulation assuming a cross lagged regression DGP. On the other hand, `simRICLPM3` 
simulates a random intercept cross lagged panel with three variables, `simRICLPM` with two. Finally , `simRICLPMl` specifices a latent variable
simulation, when more than one indicator per time point is included. All files with `sim` as a prefix are simulation files.

This directory also includes estimation functions, which is useful to compare estimators in varying DGPs. These are prefixed with `estimate`.
Thus far, only `estimateCLPM`, `estimateRICLPM`, `estimateStTr` (a trait-state model), and `estimateHLM` (a hierarchical model in `rstan`) 
have been thoroughly tested.





- `src/`: Contains the source code for the project.
- `tests/`: Contains unit tests for the project.
- `examples/`: Contains example scripts demonstrating how to use the project.
- `LICENSE`: The license file for the project.


## Project Overview
This repository contains a project that aims to [briefly describe the purpose of the project, e.g., "provide a set of tools for data analysis"]. It includes various components such as source code, documentation, tests, and examples to help users understand and utilize the project effectively.

# Description
# Installation
# Usage
# Contributing

## Technologies
- R
- Shiny
- Claude.AI - This project's code was developed with the assistance of an AI code assistant to assist with with base code, error catching, and debugging.
