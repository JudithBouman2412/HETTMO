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
  int num_serosurvey;
  array[num_t] int ts;
  int popsize;

  // priors
  vector[2] p_I0; // expected initial seed (mean, sd)
  vector[2] p_R0; // expected beta (alpha, beta)
  real p_phi;

  // fixed quantities
  real generation_time;
  array[num_serosurvey] int t_survey_start;
  array[num_serosurvey] int t_survey_end;
  array[num_serosurvey] int n_infected_survey;
  array[num_serosurvey] int n_tested_survey;

  real beta_fixed;
  array[num_serosurvey] real sens;
  array[num_serosurvey] real spec;

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

  array[8] int DIM = {num_comp, num_t, num_prev, popsize, num_basis, spline_degree, 20, num_serosurvey};

}

parameters {
  // generic parameters
  real<lower=0> R0;
  real<lower=0> I0_raw;                 // initial seed (in number of individuals)
  vector<lower=0, upper=1>[num_serosurvey] pi_;

  // paremeters dependent on type of time dependence
  vector<lower=0>[num_basis-1] alpha_init;              // coefficients for spline
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
  vector[num_t] asc_incidence;

  vector[num_basis] a;
  a[1] = (beta); // can also be negative, as we use inv.logit insdide ODE
  for (i in 2:num_basis)
    a[i] =  alpha_init[i-1];      // a[i-1] + alpha_init[i-1]*kappa;

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
    asc_incidence = get_incidence(y, DIM, popsize, atol, I0, pi_, t_survey_start, t_survey_end);

  // given the SEIR dynamics, calculate the probability to observe a seropositive individual in the population
  array[num_serosurvey] real cum_inf_frac;
  for (q in 1:num_serosurvey){
    cum_inf_frac[q] = mean(to_vector(y[t_survey_start[q]:t_survey_end[q], 4])) / popsize;
  }

}

model {
  // Priors
  I0_raw ~ gamma(p_I0[1]^2/p_I0[2]^2,p_I0[1]/p_I0[2]^2);
  R0 ~ gamma(p_R0[1],p_R0[2]);
  alpha_init ~ normal((p_R0[1]*gamma/contact)/beta_fixed, 0.1);
  pi_[i,] ~ beta(2,2);
  phi_inv ~ exponential(p_phi);

  // likelihood
  if(inference==1){
      for (i in 1:num_serosurvey){
        n_infected_survey[i] ~ binomial(n_tested_survey[i], cum_inf_frac[i]*sens[i] + (1-cum_inf_frac[i])*(1-spec[i]));
      }
     target += neg_binomial_2_lpmf( data_pre | asc_incidence, phi );
  }


}

generated quantities {


  // log-likelihood for LOO
  vector[num_t] log_lik;

  for (t in 1:num_t)
    log_lik[t] = neg_binomial_2_lpmf(data_pre[t] | asc_incidence[t], phi);
    array[num_t] int confirmed_cases_simulated = data_pre;

  // posterior predictive check
  array[num_t] int confirmed_cases_predicted = neg_binomial_2_rng(asc_incidence, phi); //  poisson_rng(incidence); Also adjust, because likelihood is calculated differently?

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


