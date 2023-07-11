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
  int t_survey_start;
  int t_survey_end;
  int n_infected_survey;
  int n_tested_survey;

  int t_detectionSwitch;
  real beta_fixed;
  real sens;
  real spec;
  real p_detect1;

  real fraction_pre;

  // load data objects for the type of stratification choosen
  real contact;

  // load data element that defines what sampler should be used
  int sampler;

}

transformed data {

  #include "/data/generic_transformed_data.stan"
  #include "/data/splines_transformed_data.stan"

  array[7] int DIM = {num_comp,num_t,num_prev, popsize, num_basis, spline_degree, 20};
}

parameters {
  // generic parameters
  real<lower=0> R0;           // transmission probability per contact
  real<lower=0> I0_raw;                 // initial seed (in number of individuals)
  //real<lower=0.4, upper=0.6> fraction_pre; // fraction of the generation time spend in compartment E
  //real<lower=0, upper=1> p_detect1; // ascertainmentrate first wave
  real<lower=0, upper=1> p_detect2; // ascertainmentrate second wave

  // paremeters dependent on type of time dependence
  vector<lower=0>[num_basis-1]  a_raw;              // coefficients for spline

  //for quasi poisson model
  real<lower=1.5> theta;
}

transformed parameters {

  #include "/data/generic_transformed_para.stan"

  // generic steps
  row_vector[num_t] incidence;

  vector[num_basis] a;
  a[1] = (beta/beta_fixed); // can also be negative, as we use inv.logit insdide ODE
  for (i in 2:num_basis)
    a[i] =  a_raw[i-1];      // a[i-1] + a_raw[i-1]*kappa;

  array[num_t] vector[num_comp] y;

  // run ODE solver --> type of solver is selected based on variable "sampler"
  if (sampler == 0 ){
    y = ode_rk45_tol(
      seir_0d_spline,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, a, b_hat, order,       // parameters
      tau, gamma, contact, beta_fixed,              // data
      DIM                               // metadata
      );
      // calculate incidence of the infection
      incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
  } else if (sampler == 1 ){
    y = ode_adams_tol(
      seir_0d_spline,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, a, b_hat, order, // parameters
      tau, gamma, contact, beta_fixed,              // data
      DIM                               // metadata
      );
    // calculate incidence of the infection
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
  } else if (sampler == 2){
    y = ode_bdf_tol(
      seir_0d_spline,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, a, b_hat, order, // parameters
      tau, gamma, contact, beta_fixed,              // data
      DIM                               // metadata
      );
    // calculate incidence of the infection
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
  } else if (sampler == 3){
    y = ode_ckrk_tol(
      seir_0d_spline,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, a, b_hat, order, // parameters
      tau, gamma, contact, beta_fixed,              // data
      DIM                               // metadata
      );
    // calculate incidence of the infection
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
  } else if (sampler == 4){
    y = solve_ode_system_trapezoidal( rep_vector(0.0,num_comp),
                                        ts,
                                        I0, knots, a, b_hat, order,
                                        tau, gamma, contact, beta_fixed,
                                        DIM
                                        );
    // calculate prevalence of the infection
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
  }

  real p_infected_survey = mean(to_vector(y[t_survey_start:t_survey_end, 4])) / popsize;

}

model {
  // Priors
  R0 ~ gamma(p_R0[1],p_R0[2]);

  I0_raw ~ gamma(p_I0[1]^2/p_I0[2]^2,p_I0[1]/p_I0[2]^2);
  a_raw ~ normal((p_R0[1]*gamma/contact)/beta_fixed, 0.1);

  theta ~ normal( p_theta[1], p_theta[2]);
  //fraction_pre ~ uniform(0.4, 0.6);
  //p_detect1 ~ normal(0.1,0.1);
  p_detect2 ~ normal(0.5,0.1);

  // quasi poisson model
  if (inference==1) {
    n_infected_survey ~ binomial(n_tested_survey,p_infected_survey*sens + (1-p_infected_survey)*(1-spec)); // do I need to combine this in the target += statement directly?
    target += neg_binomial_2_log_lpmf( data_pre | log(incidence), to_array_1d( to_row_vector(incidence)/(theta-1)) );
  }

}

generated quantities {

  // log-likelihood for LOO
  vector[num_t] log_lik;

  //  poisson
  if (inference==1) {
    for (t in 1:num_t)
        log_lik[t] = neg_binomial_2_log_lpmf( data_pre[t] | log(incidence)[t], incidence[t]/(theta-1));
        array[num_t] int I_t_simulated = data_pre;
  }

 // real log_sero = binomial_lpmf(n_infected_survey | n_tested_survey, p_infected_survey*sens + (1-p_infected_survey)*(1-spec));

  array[num_t] int I_t_predicted = qpoisson_rng(num_t, incidence, theta);

  // calculate probability of infection
  vector[num_t] prob_infection;
  for (i in 1:num_t){
    prob_infection[i] = analytical_bspline( b_hat, a, ts[i], knots, order );
  }

  // calculate effective reproduction number
  vector[num_t] R_eff;
  for (i in 1:num_t){
    R_eff[i] = (y[i,1]/popsize) * prob_infection[i] * beta_fixed  * (contact/gamma);
  }

}


