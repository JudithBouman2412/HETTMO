//
// Basic transmission model with stratification and time dependence in the form of splines
// in the variables.

functions {
  #include "/functions/generic_functions.stan"
  #include "/functions/GP_functions.stan"
  #include "/functions/SEIR_ode_GP.stan"
  #include "/functions/prevalence_noS.stan" // prevalence function without stratification
}

// load data objects
data {
  // load basic data objects that are needed independent of the exact specifications of the model
  // data objects that are generic independent of the model specifications
  int num_t;
  int num_serosurvey;
  array[num_t] real ts;
  int popsize;

  // priors
  vector[2] p_I0; // expected initial seed (mean, sd)
  vector[2] p_R0; // expected beta (alpha, beta)
  vector[2] p_theta;

  // fixed quantities
  real generation_time;
  real fraction_pre;

  // control parameters
  real rtol;
  real atol;
  int max_num_steps;
  int inference;

  // load data objects related to the choosen time dependency
  real<lower=0> c_f1;  // factor c to determine the boundary value L
  int<lower=1> M_f1;   // number of basis functions for f1

  // use pre-defined input data to fit the model
  array[num_t] int data_pre;
  array[num_serosurvey] int t_survey_start;
  array[num_serosurvey] int t_survey_end;
  array[num_serosurvey] int n_infected_survey;
  array[num_serosurvey] int n_tested_survey;

  real beta_fixed;
  array[num_serosurvey] real sens;
  array[num_serosurvey] real spec;

  // load data objects for the type of stratification choosen
  real contact;

  // load data element that defines what sampler should be used
  int sampler;

}

transformed data {

  // generic transformed data block
  int num_comp = 4; // total compartments
  int num_prev = 3; // index of the I compartment
  real t0 = 0.;

  real alpha_f1 = 0.5;

  // normalize time series
  real tsmean = mean(ts);
  //real tsds = sd(ts);
  //real ysd = sd(ts);
  array[num_t] real ts_norm = to_array_1d( (to_row_vector(ts) - tsmean)/(ts[num_t]-tsmean) ); // make sure x is centered around 0

  // compute boundary value
  real L_f1 = c_f1*ts_norm[num_t];

  // compute basis functions for f1
  matrix[num_t, M_f1] PHI_f1 = PHI_EQ(num_t, M_f1, L_f1, ts_norm); // use transformed ts data instead of the real ts

  array[7] int DIM = {num_comp, num_t, num_prev, M_f1, popsize, 20, num_serosurvey}; // Do we need any information here?

}

parameters {
  // generic parameters
  real<lower=0> R0;           // transmission probability per contact
  real<lower=0> I0_raw;                 // initial seed (in number of individuals)
  vector<lower=0, upper=1>[num_serosurvey] pi_;

  // paremeters dependent on type of time dependence
  vector[M_f1] beta_f1;         // basis function coefficients for f1
  real<lower=0> lambda_f1;      // lengthscale of f1

  //for quasi poisson model
  real<lower=1.5> theta;
}

transformed parameters {
  // generic steps
  real tau = 1./(generation_time*fraction_pre);
  real gamma = 1./(generation_time*(1-fraction_pre));

  real I0 = (1+I0_raw);
  real beta = ((R0*gamma)/contact)/beta_fixed;

  // prepare GP based on sampled variables
  vector[M_f1] diagSPD_f1 = diagSPD_EQ(alpha_f1, lambda_f1, L_f1, M_f1);

  vector[num_t] asc_incidence;

  array[num_t] vector[num_comp] y;

  // run ODE solver --> type of solver is selected based on variable "sampler"
  if (sampler == 0 ){
    y = ode_rk45_tol( seir_0d_GP,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)
      t0,                       // initial time from centered time scale
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, beta, beta_f1, diagSPD_f1, PHI_f1, // parameters
      tau, gamma, contact, beta_fixed,              // data
      ts, DIM, M_f1               // metadata
      );
    // calculate incidence of the infection
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 1 ){
    y = ode_adams_tol( seir_0d_GP,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)
      t0,                       // initial time from centered time scale
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, beta, beta_f1,  diagSPD_f1, PHI_f1, // parameters
      tau, gamma, contact, beta_fixed,              // data
      ts, DIM,  M_f1                                // metadata
      );
    // calculate incidence of the infection
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 2){
    y = ode_bdf_tol( seir_0d_GP,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)
      t0,                       // initial time from centered time scale
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, beta, beta_f1,  diagSPD_f1, PHI_f1, // parameters
      tau, gamma, contact, beta_fixed,              // data
      ts, DIM, M_f1                                // metadata
      );
    // calculate incidence of the infection
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 3){
    y = ode_ckrk_tol( seir_0d_GP,
      rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)
      t0,                       // initial time from centered time scale
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, beta, beta_f1,  diagSPD_f1, PHI_f1, // parameters
      tau, gamma, contact, beta_fixed,              // data
      ts, DIM, M_f1                            // metadata
      );
    // calculate incidence of the infection
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 4){
    y = solve_ode_system_trapezoidal(   rep_vector(0.0,num_comp),
                                        ts,
                                        I0, L_f1, M_f1, PHI_f1, beta,
                                        beta_f1,  diagSPD_f1,
                                        tau, gamma, contact, beta_fixed,
                                        DIM  );
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
  // SEIR parameters
  R0 ~ gamma(p_R0[1],p_R0[2]);
  I0_raw ~ gamma(p_I0[1]^2/p_I0[2]^2,p_I0[1]/p_I0[2]^2);
  pi_ ~ beta(1,1);

  // GP parameters
  beta_f1 ~ normal(0, 1);
  lambda_f1 ~ exponential(5);

  theta ~ normal( p_theta[1], p_theta[2]);

  // likelihood
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

  //  quasi poisson
  if (inference==1) {
    for (t in 1:num_t)
        log_lik[t] = neg_binomial_2_log_lpmf( data_pre[t] | log(asc_incidence)[t], asc_incidence[t]/(theta-1));
        array[num_t] int confirmed_cases_simulated = data_pre;
  }

  array[num_t] int confirmed_cases_predicted = qpoisson_rng(num_t, asc_incidence, theta);

  // calculate probability of infection
  array[num_t] real rho;
  // vector[M_f1] diagSPD_f1_sim = diagSPD_EQ(alpha_f1, lambda_f1, L_f1, M_f1);
  for (i in 1:num_t){
    rho[i] = ((beta) + inv_logit(to_row_vector(ana_phi(PHI_f1, M_f1, i, num_t, ts)) * (diagSPD_f1 .* beta_f1))) ;
  }

    // calculate effective reproduction number
  vector[num_t] R_eff;
  for (i in 1:num_t){
    R_eff[i] = (((y[i,1]+popsize)/popsize) * rho[i] * beta_fixed ) * (contact/gamma);
  }


}
