
#' standata
#'
#' Function to generate a list with all input variables needed to run the Stan
#' model as specified by the input variables
#'
#' @param data (vector, array) Vector (unstratified model) or array (stratified model) with the
#' (weekly) case data that the model should be fitted to.
#' @param homogeneous TRUE if an unstratified model should be fitted and FALSE if
#' the fitted model should be stratified
#' @param simulated TRUE if simulated data is used and FALSE if
#' @param type type of time-dependence that is included: "BM" (Brownian motion),
#'  "spline" (splines) or "GP" (Gaussian process)
#' @param sampling the sampling distribution that should be used: "p" (Poisson),
#' "qp" (quasi-Poisson) or "negbin" (negative binomial)
#' @param seroprev_dat dataframe of seroprevalence data if homogeneous model is used,
#' should include one column called "num_tested" with the number of tests performed and
#' one column called "num_pos_tests" with the number of positive tests
#' @param additional_data a list with additional data elements, only relevant when
#' the stratified model is used. For the simulated data, additional
#' data only contains the data from Prem et al. (2021) to create the contact matrix.
#' For the real case data it needs to include all seroprevalence datasets and a
#' vector called "pop_GE" with the population distribution over the age-classes for the
#' Canton of Geneva.
#'
#' @return List with all input variables needed to run the Stan model
#' @export
#'
#' @examples
standata <- function( data,
                      homogeneous = TRUE,
                      simulated = TRUE,
                      type = "spline",
                      sampling = "qp",
                      seroprev_dat,
                      additional_data){

  # priors general for all models
  priors <- list(p_I0 = c(1,2),
                 p_R0 = c(2,1))

  # fixed parameters general for all models
  fixed_pars <- list(generation_time=5.2/7,
                     fraction_pre=0.5) # set monthly knots for the Bern data

  # parameters specific for sampling distribution
  if (sampling=="qp"){
    priors$p_theta = c(1/10)
  } else if ( sampling =="negbin"){
    priors$p_phi = 1
  }

  if (simulated){
    # parameters specific for simulated data
    fixed_pars$sens = c(1,1)
    fixed_pars$spec = c(1,1)
    fixed_pars$ts = 1:45
    fixed_pars$num_t = 45
    fixed_pars$beta_fixed = 0.085
    fixed_pars$t_survey_start = c(20, 45) #begin_week
    fixed_pars$t_survey_end = c(20, 45) #end_week
    fixed_pars$num_serosurvey = 2

    if (homogeneous){
      # homogeneous specific parameters
      fixed_pars$contact = 77
      fixed_pars$n_tested_survey = c(5000, 5000)
      fixed_pars$n_infected_survey = seroprev_dat$n_infected_survey
      fixed_pars$popsize = 100000
      fixed_pars$data_pre=data
    } else {
      # heterogeneous specific parameters
      contact_matrix <- create_contactmatrix_GE(additional_data)
      fixed_pars$data_pre=t(data)
      fixed_pars$popdist=contact_matrix[[2]]
      fixed_pars$contact=contact_matrix[[1]]*7
      fixed_pars$n_tested_survey = cbind(c(3000,6000,6000),
                                         c(3000,6000,6000))# 5 percent of total population was tested
      fixed_pars$n_infected_survey = seroprev_dat
      fixed_pars$beta_fixed = 0.1
      fixed_pars$num_class = 3
      fixed_pars$num_age = 3
      fixed_pars$num_sex = 1
      fixed_pars$p_beta = c(1,0.1)
    }

  } else {
    # Parameters for COVID-19 data from Canton of Geneva
    fixed_pars$sens = c(0.93,0.93)
    fixed_pars$spec = c(1.0, 1.0)
    fixed_pars$ts = 1:45
    fixed_pars$num_t = 45
    fixed_pars$num_serosurvey = 2
    fixed_pars$beta_fixed = 0.085

    if (homogeneous){
      # homogeneous specific parameters
      fixed_pars$contact = 50
      fixed_pars$n_tested_survey = c(sum(seroprev_dat[[1]]$num_tested),
                                     sum(seroprev_dat[[2]]$num_tested))
      fixed_pars$n_infected_survey = c(sum(seroprev_dat[[1]]$num_pos_tests),
                                       sum(seroprev_dat[[2]]$num_pos_tests))
      fixed_pars$popsize = 499480
      fixed_pars$data_pre=data
      fixed_pars$t_survey_start = c(6, 39) #begin_week
      fixed_pars$t_survey_end = c(10, 45) #end_week
    } else {
      data(contact_all)
      contact_matrix <- create_contactmatrix_GE(contact_all)
      # heterogeneous specific parameters
      fixed_pars$data_pre=t(data)
      fixed_pars$popdist=additional_data$pop_GE
      fixed_pars$contact=contact_matrix[[1]]*7
      fixed_pars$t_survey_start = c(6, 39) #begin_week
      fixed_pars$t_survey_end = c(10, 45) #end_week

      fixed_pars$n_infected_survey = cbind(additional_data$sero_1[[1]]$num_pos_tests,
                                           additional_data$sero_2[[1]]$num_pos_tests)
      fixed_pars$n_tested_survey = cbind(additional_data$sero_1[[1]]$num_tested,
                                         additional_data$sero_2[[1]]$num_tested)
      fixed_pars$num_class = 3
      fixed_pars$num_age = 3
      fixed_pars$num_sex = 1
      fixed_pars$p_beta = c(1, 0.1)
    }
  }

  if (type=="spline"){

    fixed_pars$num_knots=13
    fixed_pars$spline_degree=3
    fixed_pars$knots=seq(0,48,4)

  } else if (type=="GP"){

    fixed_pars$M_f1 = 15
    fixed_pars$c_f1 = 1.5

  } else if (type=="BM"){
    fixed_pars$p_sigma_BM=0.1
  }


  return(c(priors, fixed_pars))
}

#' create_function_initial_values
#'
#' @param data_list (list) List of variables needed to run the Stan model as generated by
#' "standata".
#' @param homogeneous (logical) TRUE if an unstratified model should be fitted and FALSE if
#' the fitted model should be stratified.
#' @param type (string) Type of time-dependence that is included: "BM" (Brownian motion),
#'  "spline" (splines) or "GP" (Gaussian process).
#' @param sampling (string) The sampling distribution that should be used: "p" (Poisson),
#' "qp" (quasi-Poisson) or "negbin" (negative binomial).
#'
#' @return List of initial values for running the specified Stan model.
#' @export
#'
#' @examples
create_function_initial_values <- function(data_list,
                                           homogeneous = TRUE,
                                           type = "spline",
                                           sampling = "qp" ){

  initial_value_list <- list()

  if(homogeneous){
    # General parameters
    initial_value_list$R0=stats::rgamma(1,data_list$p_R0[1],data_list$p_R0[2])
    initial_value_list$I0_raw=stats::rgamma(1,data_list$p_I0[1])
    initial_value_list$pi_ = stats::rbeta(data_list$num_serosurvey, 2, 2)

    # parameters based on method of time-dependence
    if (type=="spline"){
      initial_value_list$alpha_init=stats::rgamma(data_list$num_knots+1, 2.5, 5)
    } else if (type=="GP"){
      initial_value_list$beta_f1 = rnorm(data_list$M_f1)
      initial_value_list$lambda_f1 = rexp(1,5)
    } else if (type=="BM"){
      initial_value_list$sigmaBM=truncnorm::rtruncnorm(1, 0, a=0, b=Inf, data_list$p_sigma_BM)
      initial_value_list$eta_noise=stats::rnorm(data_list$num_t-1,0,1)
    }
  } else {
    initial_value_list$beta=rep(truncnorm::rtruncnorm(1,0,1,data_list$p_beta[1], data_list$p_beta[2]),3)
    initial_value_list$I0_raw=stats::rgamma(1,data_list$p_I0[1])
    initial_value_list$alpha=rbind( stats::rgamma(data_list$num_knots+2, 2.5, 5),
                                    stats::rgamma(data_list$num_knots+2, 2.5, 5),
                                    stats::rgamma(data_list$num_knots+2, 2.5, 5))
    initial_value_list$pi_ = rbind(stats::rbeta(data_list$num_serosurvey, 1, 1),
                                   stats::rbeta(data_list$num_serosurvey, 1, 1),
                                   stats::rbeta(data_list$num_serosurvey, 1, 1))
  }

  # parameters specific for sampling distribution
  if (sampling=="qp"){
    initial_value_list$theta=rexp(1, data_list$p_theta ) + 2
  } else if ( sampling =="negbin"){
    initial_value_list$phi_inv = rexp(1, data_list$p_phi)
  }

  return(initial_value_list)
}

#' sampling_model
#'
#' Function that fits the selected stan model to the supplied case data.
#'
#' @param data_list (list) List of variables for the selected Stan model as generated by "standata".
#' @param model (string) The selected model for fitting the data: default option is "qpoisson_spline_unstrat.stan".
#' Other options: "poisson_spline_unstrat", "negbin_spline_unstrat", "qpoisson_BM_unstrat",
#' "qpoisson_GP_unstrat", "qpoisson_spline_strat" and "qpoisson_spline_strat_GE".
#' @param save_result (string) The path to the file where the result of the model fit should be saved
#' @param type (string) type of time-dependence that is included: "BM" (Brownian motion),
#'  "spline" (splines) or "GP" (Gaussian process).
#' @param sampling (string) The sampling distribution that should be used: "p" (Poisson),
#' "qp" (quasi-Poisson) or "negbin" (negative binomial).
#' @param homogeneous (logical) TRUE if an unstratified model should be fitted and FALSE if
#' the fitted model should be stratified.
#' @param sampler (integer) Choice of sampler used to solve the ODE-system within the Stan model,
#' options are: (0) rk45, (1) adams, (2) bdf, (3) ckrk and (4) Trapeziodal solver.
#' @param rtol (real) The relative tolerance for solvers 0-3.
#' @param atol (real) The absoluted tolerance for solvers 0-3.
#' @param max_num_steps (positive interger) Maximum number of steps in mcmc.
#'
#' @return (list) List with prior samples (samples_prior), posterior samples (samples_posterior),
#' data list used for analysis (data_list) and diagnostics of the model fit (diagnostics).
#' @export
#'
#' @examples
sampling_model <- function( data_list,
                            model = "qpoisson_spline_unstrat.stan",
                            save_result,
                            type = "spline",
                            sampling = "qp",
                            homogeneous = TRUE,
                            sampler = 4,
                            rtol=1e-6,
                            atol = 1e-6,
                            max_num_steps=1000) {

  # Add MCMC specific parameters to data_list
  data_list_within = within(data_list,{
    sampler = sampler;
    rtol = rtol;
    atol = rtol;
    max_num_steps = max_num_steps})

  # Load stan model with cmdstanr
  stan_file <- system.file( model, package = "HETTMO")
  model_compiled = cmdstanr::cmdstan_model(stan_file = stan_file ,
                                        #include_paths = "inst",
                                        force_recompile = TRUE)

  # prior predictive check
  data_list_within$inference = 0

  samples_prior = model_compiled$sample(
    data = data_list_within,
    seed = 123,
    chains = 1, parallel_chains = 1,
    iter_warmup = 500, iter_sampling = 500,
    max_treedepth = 1e2)

  #samples_prior$save_output_files(dir="various_tests/compareLikihood/")

  # create function to sample from initial conditions
  init_fun <- function() { create_function_initial_values(data_list_within, homogeneous = homogeneous,
                                                          type = type, sampling = sampling )  }

  # sample from posterior
  data_list_within$inference = 1

  if ( homogeneous){
    samples_posterior = model_compiled$sample(
      data = data_list_within,
      chains = 4, parallel_chains = 4,
      iter_warmup = 500, iter_sampling = 500,
      init=init_fun)
  } else {
    samples_posterior = model_compiled$sample(
      data = data_list_within,
      chains = 4, parallel_chains = 4,
      iter_warmup = 500, iter_sampling = 500,
      init=init_fun, adapt_delta = 0.95 )
  }

  #samples_posterior$save_output_files("various_tests/compareLikihood/")

  result <- list(samples_posterior=samples_posterior,
                                samples_prior=samples_prior,
                                data_list=data_list_within,
                                diagnostics=samples_posterior$cmdstan_diagnose())

  saveRDS(result, file = save_result)

  return(result)
}
