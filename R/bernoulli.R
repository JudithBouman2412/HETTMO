#' Example function to run the example Bernoulli model
#'
#' @param y The binary outcome of the trials.
#' @return The draws for the Bernoulli chance-of-success `theta` parameter.
#'
#' @export
bernoulli <- function(y) {

  stan_file <- system.file("bernoulli.stan", package = "democmdstanr")

  data <- list(
    N = length(y),
    y = y
  )
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(data = data, chains = 1)
  return(fit$draws("theta"))
}

