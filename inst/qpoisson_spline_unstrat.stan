//
// Basic transmission model with stratification and time dependence in the form of splines
// in the variables.

functions {
  #include "/functions/generic_functions.stan"
  #include "/functions/spline_functions.stan"
  #include "/functions/SEIR_ode_spline.stan"
  #include "/functions/prevalence_noS.stan"
}

// load data objects
data {
  // load basic data objects that are needed independent of the exact specifications of the model
  #include "/data/generic_data.stan"

  // load data objects related to the choosen time dependency
  #include "/data/splines_data.stan"

  // use pre-defined input data to fit the model
  array[num_t] int data_pre;
  array[num_serosurvey] int t_survey_start;
  array[num_serosurvey] int t_survey_end;
  array[num_serosurvey] int n_infected_survey;
  array[num_serosurvey] int n_tested_survey;

  real beta_fixed;
  array[num_serosurvey] real sens;
  array[num_serosurvey] real spec;

  real fraction_pre;

  // load data objects for the type of stratification choosen
  real contact;

  // load data element that defines what sampler should be used
  int sampler;

}

transformed data {

  #include "/data/generic_transformed_data.stan"
  #include "/data/splines_transformed_data.stan"

  array[8] int DIM = {num_comp,num_t,num_prev, popsize, num_basis, spline_degree, 20, num_serosurvey};
}

parameters {
  // generic parameters
  real<lower=0> R0;           // transmission probability per contact
  real<lower=0> I0_raw;                 // initial seed (in number of individuals)
  vector<lower=0, upper=1>[num_serosurvey] pi_;

  // paremeters dependent on type of time dependence
  vector<lower=0>[num_basis-1]  alpha_init;              // coefficients for spline

  //for quasi poisson model
  real<lower=2> theta;
}

transformed parameters {

  #include "/data/generic_transformed_para.stan"

  // generic steps
  vector[num_t] asc_incidence;

  vector[num_basis] alpha;
  alpha[1] = (beta/beta_fixed); // can also be negative, as we use inv.logit insdide ODE
  for (i in 2:num_basis)
    alpha[i] =  alpha_init[i-1];      // a[i-1] + alpha_init[i-1]*kappa;

  array[num_t] vector[num_comp] y;

  // run ODE solver --> type of solver is selected based on variable "sampler"
  if (sampler == 0 ){
    y = ode_rk45_tol(
      seir_0d_spline,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, alpha, b_hat, order,       // parameters
      tau, gamma, contact, beta_fixed,              // data
      DIM                               // metadata
      );
      // calculate asc_incidence of the infection
      asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 1 ){
    y = ode_adams_tol(
      seir_0d_spline,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, alpha, b_hat, order, // parameters
      tau, gamma, contact, beta_fixed,              // data
      DIM                               // metadata
      );
    // calculate incidence of the infection
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 2){
    y = ode_bdf_tol(
      seir_0d_spline,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, alpha, b_hat, order, // parameters
      tau, gamma, contact, beta_fixed,              // data
      DIM                               // metadata
      );
    // calculate incidence of the infection
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 3){
    y = ode_ckrk_tol(
      seir_0d_spline,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, alpha, b_hat, order, // parameters
      tau, gamma, contact, beta_fixed,              // data
      DIM                               // metadata
      );
    // calculate incidence of the infection
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 4){
    y = solve_ode_system_trapezoidal( rep_vector(0.0,num_comp),
                                        ts,
                                        I0, knots, alpha, b_hat, order,
                                        tau, gamma, contact, beta_fixed,
                                        DIM
                                        );
    // calculate prevalence of the infection
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  }

  // given the SEIR dynamics, calculate the probability to observe a seropositive individual in the population
  array[num_serosurvey] real cum_inf_frac;
  for (q in 1:num_serosurvey){
    cum_inf_frac[q] = mean(to_vector(y[t_survey_start[q]:t_survey_end[q], 4])) / popsize;
  }

}

model {
  // Priors
  R0 ~ gamma(p_R0[1],p_R0[2]);

  I0_raw ~ gamma(p_I0[1]^2/p_I0[2]^2,p_I0[1]/p_I0[2]^2);
  alpha_init ~ gamma(2.5, 5); // C95 for rho(t): 0.08312116 1.28325020

  theta ~ exponential( p_theta );

  pi_ ~ beta(2,2);

  // quasi poisson model
  if (inference==1) {
    for (i in 1:num_serosurvey){
      n_infected_survey[i] ~ binomial(n_tested_survey[i], cum_inf_frac[i]*sens[i] + (1-cum_inf_frac[i])*(1-spec[i]));
    }
    target += neg_binomial_2_log_lpmf( data_pre | log(asc_incidence), to_array_1d( to_row_vector(asc_incidence)/(theta-1)) );
  }

}

generated quantities {

  // log-likelihood for LOO
  vector[num_t] log_lik;

  //  poisson
  if (inference==1) {
    for (t in 1:num_t)
        log_lik[t] = neg_binomial_2_log_lpmf( data_pre[t] | log(asc_incidence)[t], asc_incidence[t]/(theta-1));
        array[num_t] int confirmed_cases_simulated = data_pre;
  }

 // real log_sero = binomial_lpmf(n_infected_survey | n_tested_survey, p_infected_survey*sens + (1-p_infected_survey)*(1-spec));

  array[num_t] int confirmed_cases_predicted = qpoisson_rng(num_t, asc_incidence, theta);

  // calculate probability of infection
  vector[num_t] rho;
  for (i in 1:num_t){
    rho[i] = analytical_bspline( b_hat, alpha, ts[i], knots, order );
  }

  // calculate effective reproduction number
  vector[num_t] R_eff;
  for (i in 1:num_t){
    R_eff[i] = (((y[i,1]+popsize)/popsize) * rho[i] * beta_fixed ) * (contact/gamma);
  }

}


