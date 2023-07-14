#' Save our model results without cmdstan dependency
#'
#' @param stan_result model result that needs to be saved
#' @param base_name base of naming the result files
#' @param dir directory where files should be saved
#'
#' @return
#' @export
#'
#' @examples
save_stan_result <- function(stan_result, base_name, dir){
  temp_file_posterior <- paste0(dir, base_name, "_posterior.RDS" )
  posterior <- stan_result$samples_posterior
  posterior$save_object(file = temp_file_posterior)
  temp_file_prior <- paste0(dir, base_name, "_prior.RDS" )
  prior <- stan_result$samples_prior
  prior$save_object(file = temp_file_prior)
  saveRDS(stan_result$data_list, file = paste0(dir, base_name, "_datalist.RDS"))
  saveRDS(stan_result$diagnostics, file = paste0(dir, base_name, "_diagnostics.RDS"))
}
#' Read in model results
#'
#' @param base_name
#' @param dir data directory where the object is saved
#'
#' @return list of posterior samples, prior samples, model parameters and diagnostics
#' @export
#'
#' @examples
read_stan_result <- function(base_name, dir ){
  samples_posterior <- readRDS( file =  paste0(dir, base_name, "_posterior.RDS" ) )
  samples_prior <- readRDS( file =  paste0(dir, base_name, "_prior.RDS" ) )
  data_list <- readRDS( file =  paste0(dir, base_name, "_datalist.RDS" ) )
  diagnostics <- readRDS( file =  paste0(dir, base_name, "_diagnostics.RDS" ) )
  return(list(samples_posterior=samples_posterior, samples_prior=samples_prior, data_list=data_list, diagnostics=diagnostics))
}
