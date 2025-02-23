% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphRI.R
\name{riclpm_graph}
\alias{riclpm_graph}
\title{riclpm_graph}
\usage{
riclpm_graph(model = fit, waves = c(3, 4, 5, 6))
}
\arguments{
\item{waves}{Number of waves in the fitted model.}

\item{fit}{A fitted lavaan model object.}
}
\value{
RI-CLPM SEM plot
}
\description{
This function draws a tidySEM graph of a RI-CLPM model

Create a simulated dataset of a cross-lagged model with a specified number of waves and structural parameters.
}
