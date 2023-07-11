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
  array[num_t] real ts;
  int popsize;

  // priors
  vector[2] p_I0; // expected initial seed (mean, sd)
  vector[2] p_R0; // expected beta (alpha, beta)
  vector[2] p_theta;

  // fixed quantities
  real generation_time;

  real beta_fixed;
  real sens;
  real spec;
  real p_detect1;
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
  int t_survey_start;
  int t_survey_end;
  int n_infected_survey;
  int n_tested_survey;

  int t_detectionSwitch;

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

  array[6] int DIM = {num_comp, num_t, num_prev, M_f1, popsize, 20}; // Do we need any information here?

}

parameters {
  // generic parameters
  real<lower=0> R0;
  real<lower=0> I0_raw;                 // initial seed (in number of individuals)
  //real<lower=0.4, upper=0.6> fraction_pre; // fraction of the generation time spend in compartment E
  //real<lower=0, upper=1> p_detect1; // ascertainmentrate first wave
  real<lower=0, upper=1> p_detect2; // ascertainmentrate second wave

  // paremeters dependent on type of time dependence
  vector[M_f1] beta_f1;         // basis function coefficients for f1
  real<lower=0> lambda_f1;      // lengthscale of f1
  //real<lower=0> alpha_f1;       // scale of f1

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

  row_vector[num_t] incidence;

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
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
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
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
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
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
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
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
  } else if (sampler == 4){
    y = solve_ode_system_trapezoidal(   rep_vector(0.0,num_comp),
                                        ts,
                                        I0, L_f1, M_f1, PHI_f1, beta,
                                        beta_f1,  diagSPD_f1,
                                        tau, gamma, contact, beta_fixed,
                                        DIM  );
    // calculate prevalence of the infection
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
  }

  real p_infected_survey = mean(to_vector(y[t_survey_start:t_survey_end, 4])) / popsize;

}

model {
  // Priors
  // SEIR parameters
  R0 ~ gamma(p_R0[1],p_R0[2]);
  //R0 ~ uniform(0.3,0.7);
  I0_raw ~ gamma(p_I0[1]^2/p_I0[2]^2,p_I0[1]/p_I0[2]^2);
  //fraction_pre ~ uniform(0.4, 0.6);
  //p_detect1 ~ normal(0.1,0.1);
  p_detect2 ~ beta(2,2);

  // GP parameters
  beta_f1 ~ normal(0, 1);
  lambda_f1 ~ exponential(5); // scale the data?
  //alpha_f1 ~ exponential(1);

  theta ~ normal( p_theta[1], p_theta[2]);

  // likelihood
  if (inference==1) {
    n_infected_survey ~ binomial(n_tested_survey, p_infected_survey*sens + (1-p_infected_survey)*(1-spec)); // do I need to combine this in the target += statement directly?
    target += neg_binomial_2_log_lpmf( data_pre | log(incidence), to_array_1d( to_row_vector(incidence)/(theta-1)) );
  }
}

generated quantities {

  // log-likelihood for LOO
  vector[num_t] log_lik;

  //  quasi poisson
  if (inference==1) {
    for (t in 1:num_t)
        log_lik[t] = neg_binomial_2_log_lpmf( data_pre[t] | log(incidence)[t], incidence[t]/(theta-1));
        array[num_t] int I_t_simulated = data_pre;
  }

  real log_sero = binomial_lpmf(n_infected_survey | n_tested_survey, p_infected_survey*sens + (1-p_infected_survey)*(1-spec));

  array[num_t] int I_t_predicted = qpoisson_rng(num_t, incidence, theta);

  // calculate probability of infection
  array[num_t] real prob_infection;
  // vector[M_f1] diagSPD_f1_sim = diagSPD_EQ(alpha_f1, lambda_f1, L_f1, M_f1);
  for (i in 1:num_t){
    prob_infection[i] = ((beta) + inv_logit(to_row_vector(ana_phi(PHI_f1, M_f1, i, num_t, ts)) * (diagSPD_f1 .* beta_f1))) ;
  }

    // calculate effective reproduction number
  vector[num_t] R_eff;
  for (i in 1:num_t){
    R_eff[i] = (y[i,1]/popsize) * prob_infection[i] * beta_fixed * (contact/gamma);
  }


}
