# Package index

## Estimators (lavaan syntax generators)

Functions that emit lavaan model syntax for panel-data models.

- [`estimateCLPM()`](https://crweber9874.github.io/crossLagR/reference/estimateCLPM.md)
  : estimateCLPM. This follows the Usami, Murayama, & Hamaker (2019)
  unified framework for longitudinal models.
- [`estimateRICLPM()`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPM.md)
  : estimateRICLPM
- [`estimateRICLPM3()`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPM3.md)
  : estimateRICLPM3
- [`estimateRICLPMl()`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPMl.md)
  : estimateRICLPMl
- [`estimateRICLPM_nolag()`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPM_nolag.md)
  : estimateRICLPM_nolag
- [`estimateRI()`](https://crweber9874.github.io/crossLagR/reference/estimateRI.md)
  : estimateRI
- [`estimateLChange()`](https://crweber9874.github.io/crossLagR/reference/estimateLChange.md)
  : estimateLChange
- [`estimateLGM()`](https://crweber9874.github.io/crossLagR/reference/estimateLGM.md)
  : estimateLGM
- [`estimateALT()`](https://crweber9874.github.io/crossLagR/reference/estimateALT.md)
  : estimateALT, following Usami, Murayama, & Hamaker (2019) unified
  frameowk.
- [`estimateLCMSR()`](https://crweber9874.github.io/crossLagR/reference/estimateLCMSR.md)
  : estimateLCMSR. This is just the GCLPM. Structure the means here.
- [`estimateBollen_and_Brand()`](https://crweber9874.github.io/crossLagR/reference/estimateBollen_and_Brand.md)
  : estimateBollen_and_Brand
- [`estimateTSO()`](https://crweber9874.github.io/crossLagR/reference/estimateTSO.md)
  : estimateTSO
- [`estimateStTr()`](https://crweber9874.github.io/crossLagR/reference/estimateStTr.md)
  : estimateStTr
- [`estimateAllisonChamberlainFI()`](https://crweber9874.github.io/crossLagR/reference/estimateAllisonChamberlainFI.md)
  : estimateAllisonChamberlainFI
- [`estimateCLPM_MI()`](https://crweber9874.github.io/crossLagR/reference/estimateCLPM_MI.md)
  : estimateCLPM_MI
- [`estimateRICLPM_MI()`](https://crweber9874.github.io/crossLagR/reference/estimateRICLPM_MI.md)
  : estimateRICLPM_MI
- [`estimateALT_MI()`](https://crweber9874.github.io/crossLagR/reference/estimateALT_MI.md)
  : estimateALT_MI
- [`estimateLCMSR_MI()`](https://crweber9874.github.io/crossLagR/reference/estimateLCMSR_MI.md)
  : estimateLCMSR_MI
- [`estimateBollen_and_Brand_MI()`](https://crweber9874.github.io/crossLagR/reference/estimateBollen_and_Brand_MI.md)
  : estimateBollen_and_Brand_MI

## Data simulators

Generate panel data under known parameters.

- [`simCLPM()`](https://crweber9874.github.io/crossLagR/reference/simCLPM.md)
  : simCLPM
- [`simCLPMu()`](https://crweber9874.github.io/crossLagR/reference/simCLPMu.md)
  : simCLPMu
- [`simCLPM_timeInvariantU()`](https://crweber9874.github.io/crossLagR/reference/simCLPM_timeInvariantU.md)
  : simCLPM_timeInvariantU
- [`simRICLPM()`](https://crweber9874.github.io/crossLagR/reference/simRICLPM.md)
  : simRICLPM
- [`simRICLPM3()`](https://crweber9874.github.io/crossLagR/reference/simRICLPM3.md)
  : simRICLPM
- [`simRICLPMl()`](https://crweber9874.github.io/crossLagR/reference/simRICLPMl.md)
  : simulateRICLPMl
- [`simLChange()`](https://crweber9874.github.io/crossLagR/reference/simLChange.md)
  : simLChange
- [`simLGM()`](https://crweber9874.github.io/crossLagR/reference/simLGM.md)
  : simLGM
- [`simWithinBetween()`](https://crweber9874.github.io/crossLagR/reference/simWithinBetween.md)
  : simWithinBetween
- [`sim_trait_state()`](https://crweber9874.github.io/crossLagR/reference/sim_trait_state.md)
  : sim_trait_state
- [`baseRICLPM_sim()`](https://crweber9874.github.io/crossLagR/reference/baseRICLPM_sim.md)
  : Simulate Data from a Random Intercept Cross-Lagged Panel Model
  (RI-CLPM)

## Monte Carlo runners

Replication wrappers that fit estimators across many simulated datasets.

- [`monteCarloAllisonChamberlainFI()`](https://crweber9874.github.io/crossLagR/reference/monteCarloAllisonChamberlainFI.md)
  : monteCarloAllisonChamberlainFI
- [`monteCarloCLPM()`](https://crweber9874.github.io/crossLagR/reference/monteCarloCLPM.md)
  : monteCarloCLPM
- [`monteCarloCTSEM()`](https://crweber9874.github.io/crossLagR/reference/monteCarloCTSEM.md)
  : monteCarloCTSEM
- [`monteCarloConfounder()`](https://crweber9874.github.io/crossLagR/reference/monteCarloConfounder.md)
  : monteCarloConfounder
- [`monteCarloFixed()`](https://crweber9874.github.io/crossLagR/reference/monteCarloFixed.md)
  : monteCarloFixed
- [`monteCarloLChange()`](https://crweber9874.github.io/crossLagR/reference/monteCarloLChange.md)
  : monteCarloLChange
- [`monteCarloOLS()`](https://crweber9874.github.io/crossLagR/reference/monteCarloOLS.md)
  : monteCarloOLS
- [`monteCarloRI()`](https://crweber9874.github.io/crossLagR/reference/monteCarloRI.md)
  : monteCarloRI
- [`monteCarloRICLPM()`](https://crweber9874.github.io/crossLagR/reference/monteCarloRICLPM.md)
  : monteCarloRICLPM
- [`monteCarloTraitState()`](https://crweber9874.github.io/crossLagR/reference/monteCarloTraitState.md)
  : monteCarloTraitState
- [`monteCarlo_OLS()`](https://crweber9874.github.io/crossLagR/reference/monteCarlo_OLS.md)
  : simpleSim
- [`run_mc_sims()`](https://crweber9874.github.io/crossLagR/reference/run_mc_sims.md)
  : run_mc_sims

## Helpers

Utilities for reshaping, plotting, and continuous-time SEM.

- [`reshape_long_sim_cr()`](https://crweber9874.github.io/crossLagR/reference/reshape_long_sim_cr.md)
  : reshape_long_sim_cr
- [`compCTSEM()`](https://crweber9874.github.io/crossLagR/reference/compCTSEM.md)
  : compCTSEM
- [`riclpm_graph()`](https://crweber9874.github.io/crossLagR/reference/riclpm_graph.md)
  : riclpm_graph
- [`withinBetween()`](https://crweber9874.github.io/crossLagR/reference/withinBetween.md)
  : withinBetween
- [`plotWithinBetween()`](https://crweber9874.github.io/crossLagR/reference/plotWithinBetween.md)
  : plotWithinBetween
- [`zero.one()`](https://crweber9874.github.io/crossLagR/reference/zero.one.md)
  : Create a zero one normalization
- [`generate_riclpm_indicator_lists_for_estimation()`](https://crweber9874.github.io/crossLagR/reference/generate_riclpm_indicator_lists_for_estimation.md)
  : generate_riclpm_indicator_lists_for_estimation
