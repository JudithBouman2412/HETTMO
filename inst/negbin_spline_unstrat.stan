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
  // data objects that are generic independent of the model specifications
  int num_t;
  array[num_t] int ts;
  int popsize;

  // priors
  vector[2] p_I0; // expected initial seed (mean, sd)
  vector[2] p_R0; // expected beta (alpha, beta)
  vector[2] p_theta;

  // fixed quantities
  real generation_time;
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

  // control parameters
  real rtol;
  real atol;
  int max_num_steps;
  int inference;

  // load data objects related to the choosen time dependency

  // define periods
  int num_knots;
  int spline_degree;
  vector[num_knots] knots;

  // use pre-defined input data to fit the model
  array[num_t] int data_pre;

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

  #include "/data/splines_transformed_data.stan"

  array[7] int DIM = {num_comp, num_t, num_prev, popsize, num_basis, spline_degree,20};

}

parameters {
  // generic parameters
  // real<lower=0, upper=1> beta;           // transmission probability per contact
  real<lower=0> R0;
  real<lower=0> I0_raw;                 // initial seed (in number of individuals)
  //real<lower=0.4, upper=0.6> fraction_pre; // fraction of the generation time spend in compartment E
  //real<lower=0, upper=1> p_detect1; // ascertainmentrate first wave
  real<lower=0, upper=1> p_detect2; // ascertainmentrate second wave

  // paremeters dependent on type of time dependence
  vector<lower=0>[num_basis-1] a_raw;              // coefficients for spline
  real<lower=0> phi_inv;

}

transformed parameters {
  // generic steps
  real tau = 1./(generation_time*fraction_pre);
  real gamma = 1./(generation_time*(1-fraction_pre));

  real I0 = (1+I0_raw);

  // calculate beta from R0
  real beta = ((R0*gamma)/contact)/beta_fixed;

  real phi = 1./phi_inv;
  row_vector[num_t] incidence;

  vector[num_basis] a;
  a[1] = (beta); // can also be negative, as we use inv.logit insdide ODE
  for (i in 2:num_basis)
    a[i] =  a_raw[i-1];      // a[i-1] + a_raw[i-1]*kappa;

  // run ODE solver --> type of solver is selected based on variable "sampler"
  array[num_t] vector[num_comp] y = ode_rk45_tol(
    seir_0d_spline,
    rep_vector(0.0,num_comp),         // initial values = 0 (handled within the ODE)  to_vector(append_array( {log(1-I0), -50, log(I0)}, rep_array(-50, num_comp-3) ) )
    t0,                               // initial time = 0
    ts,                               // evaluation times
    rtol, atol, max_num_steps,        // tolerances
    I0, knots, a, b_hat, order,       // parameters
    tau, gamma, contact, beta_fixed,             // data
    DIM                               // metadata
    );
    incidence = get_incidence(y, DIM, popsize, atol, I0, p_detect1, p_detect2, t_detectionSwitch);
    real p_infected_survey = mean(to_vector(y[t_survey_start:t_survey_end, 4])) / popsize;

}

model {
  // Priors
  I0_raw ~ gamma(p_I0[1]^2/p_I0[2]^2,p_I0[1]/p_I0[2]^2);
  //beta ~ beta(p_beta[1],p_beta[2]);
  R0 ~ gamma(p_R0[1],p_R0[2]);
  a_raw ~ normal((p_R0[1]*gamma/contact)/beta_fixed, 0.1);
  //fraction_pre ~ uniform(0.4, 0.6);
  //p_detect1 ~ normal(0.1,0.1);
  p_detect2 ~ normal(0.5,0.1);
  phi_inv ~ exponential(p_phi);

  // likelihood
  if(inference==1){
      n_infected_survey ~ binomial(n_tested_survey,p_infected_survey*sens + (1-p_infected_survey)*(1-spec));
     target += neg_binomial_2_lpmf( data_pre | incidence, phi );
  }


}

generated quantities {


  // log-likelihood for LOO
  vector[num_t] log_lik;

  for (t in 1:num_t)
    log_lik[t] = neg_binomial_2_lpmf(data_pre[t] | incidence[t], phi);
    array[num_t] int I_t_simulated = data_pre;

  // posterior predictive check
  array[num_t] int I_t_predicted = neg_binomial_2_rng(incidence, phi); //  poisson_rng(incidence); Also adjust, because likelihood is calculated differently?

  // calculate probability of infection
  array[num_t] real<lower=0, upper=1> prob_infection;
  for (i in 1:num_t){
    prob_infection[i] = inv_logit( analytical_bspline( b_hat, a, ts[i], knots, order ));
  }

  // calculate effective reproduction number
  vector[num_t] R_eff;
  for (i in 1:num_t){
    R_eff[i] = prob_infection[i]*(contact/(1/(generation_time*(1-fraction_pre))));
  }



}


