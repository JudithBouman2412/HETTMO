% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_spline_functions.R
\name{find_poly}
\alias{find_poly}
\title{find_poly}
\usage{
find_poly(x, y)
}
\arguments{
\item{x}{x values}

\item{y}{y values}
}
\value{
coefficients of the fitted polynomial
}
\description{
Use the Lagrange algorithm to fit a polynomial of degree (number of points -1)
through a set of points.
Code inspired by: https://en.wikibooks.org/wiki/Algorithm_Implementation/Mathematics/Polynomial_interpolation#:~:text=Lagrange\%20interpolation\%20is\%20an\%20algorithm,xi\%2C\%20yi).
}
