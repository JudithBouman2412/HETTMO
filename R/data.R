#' Synthetic contact matrices per country by Prem et al. (2021)
#'
#' @format ## `contact_all`
#' A list of 177 arrays
#' \describe{
#'   \item{country}{synthetic contact matrix per country}
#' }
#' @source <https://doi.org/10.1371/journal.pcbi.1009098>
"contact_all"


#' SARS-CoV-2 related data for the Canton of Geneva, Switzerland
#' in 2020.
#'
#' @format ## `GE_data`
#' List of 5 items
#' \describe{
#'   \item{seroprevalence study 1}{number of tests performed, number of positive tests per age-group}
#'   \item{seroprevalence study 2}{number of tests performed, number of positive tests per age-group}
#'   \item{seroprevalence study 3}{number of tests performed, number of positive tests per age-group}
#'   \item{weekly stratified cases}{Number of positive PCR tests per week per age-class in Geneva during 2020}
#'   \item{population distribution Geneva}{Number of inhabitants per age-class}
#' }
#' @source <https://doi.org/10.1371/journal.pcbi.1009098>
"GE_data"


#' Simulated non-stratified data, using the parameters as specified in the manuscript (function "set_parameters")
#'
#' @format ## `simulated_nonstratified`
#' List of 4 items
#' \describe{
#'   \item{sim_data}{simulated positive PCR tests per week}
#'   \item{re_sim}{Effective reproduction number per week}
#'   \item{n_infected_survey}{number of positive serological tests}
#'   \item{prob_trans_sim}{Rho(t) per week}
#' }
#' @source
"simulated_nonstratified"

#' Simulated stratified data, using the parameters as specified in the manuscript (function "set_parameters")
#'
#' @format ## `simulated_stratified`
#' List of 3 items
#' \describe{
#'   \item{sim_data}{array simulated positive PCR tests per week and age-class}
#'   \item{beta_sim}{Rho(t) per week and age-class}
#'   \item{n_infected_survey}{number of positive serological tests per age-class}
#' }
#' @source
"simulated_stratified"
