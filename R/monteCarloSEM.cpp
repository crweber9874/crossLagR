#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Function to simulate data based on RICLPM
arma::mat simulateRICLPM(int trials, int waves, double within_person_stability_y, double within_person_stability_x,
                         double cross_lag_x, double cross_lag_y, double variance_q, double variance_p,
                         double variance_between_x, double variance_between_y, double cov_pq, int sample_size) {

  arma::mat results(trials * sample_size, 2 * waves); // Initialize results matrix

  for (int trial = 0; trial < trials; ++trial) {
    arma::vec between_x = arma::randn(sample_size) * std::sqrt(variance_between_x);
    arma::vec between_y = arma::randn(sample_size) * std::sqrt(variance_between_y);

    for (int i = 0; i < sample_size; ++i) {
      arma::vec x(waves);
      arma::vec y(waves);

      // Initial values
      x(0) = between_x(i) + arma::randn() * std::sqrt(variance_p);
      y(0) = between_y(i) + arma::randn() * std::sqrt(variance_q);

      for (int t = 1; t < waves; ++t) {
        x(t) = within_person_stability_x * x(t - 1) + cross_lag_y * y(t - 1) + arma::randn() * std::sqrt(variance_p);
        y(t) = within_person_stability_y * y(t - 1) + cross_lag_x * x(t - 1) + arma::randn() * std::sqrt(variance_q);
      }

      // Store results for this individual and trial
      for (int t = 0; t < waves; ++t) {
        results(trial * sample_size + i, t) = x(t);
        results(trial * sample_size + i, t + waves) = y(t);
      }
    }
  }

  return results;
}

// [[Rcpp::export]]
arma::mat monteCarloCTSEM(int trials = 1, int waves = 5, std::string dgp = "riclpm",
                          double within_person_stability_y = 0.2, double within_person_stability_x = 0.2,
                          double cross_lag_x = 0.0, double cross_lag_y = 0.0,
                          double variance_q = 0.5, double variance_p = 0.5,
                          double variance_between_x = 0.5, double variance_between_y = 0.5,
                          double cov_pq = 0, int sample_size = 5000) {

  arma::mat results;

  if (dgp == "riclpm") {
    results = simulateRICLPM(trials, waves, within_person_stability_y, within_person_stability_x,
                             cross_lag_x, cross_lag_y, variance_q, variance_p,
                             variance_between_x, variance_between_y, cov_pq, sample_size);
  } else {
    Rcerr << "Unsupported DGP: " << dgp << std::endl;
    return results; // Return empty matrix if dgp is not supported
  }

  return results;
}
