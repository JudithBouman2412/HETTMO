//
// Basic transmission model with stratification and time dependence in the form of splines
// in the variables.

functions {
  #include "/functions/generic_functions.stan"
  #include "/functions/spline_functions.stan"
  #include "/functions/SEIR_ode_spline.stan"
  #include "/functions/prevalence_S_simple.stan" // prevalence function without stratification
}

// load data objects
data {
  // load basic data objects that are needed independent of the exact specifications of the model
  // data objects that are generic independent of the model specifications
  int num_t;
  int num_class;
  int num_serosurvey;
  array[num_t] real ts;

  // priors
  vector[2] p_I0; // expected initial seed (mean, sd)
  vector[2] p_beta; // expected beta (alpha, beta)
  vector[2] p_theta;

  // fixed quantities
  real generation_time;
  real fraction_pre;
  array[num_serosurvey] int t_survey_start;
  array[num_serosurvey] int t_survey_end;
  array[num_class, num_serosurvey] int n_infected_survey;
  array[num_class, num_serosurvey] int n_tested_survey;

  real beta_fixed;
  array[num_serosurvey] real sens;
  array[num_serosurvey] real spec;

  // control parameters
  real rtol;
  real atol;
  int max_num_steps;
  int inference;

  #include "/data/stratified_data.stan"

  // load data objects related to the choosen time dependency
  #include "/data/splines_data.stan"

  // use pre-defined input data to fit the model
  array[num_class, num_t] int data_pre;

   // load data element that defines what sampler should be used
  int sampler;

}

transformed data {
  #include "/data/generic_transformed_data.stan"
  #include "/data/splines_transformed_data.stan"

  // specific parameters for stratified version
  int num_eq = num_comp*num_class;

  array[9] int DIM = {num_comp, num_t, num_prev, num_age, num_sex, num_class, 20, num_eq, num_serosurvey};

}

parameters {
  // generic parameters
  real<lower=0> I0_raw;                 // initial seed (in number of individuals)

  // paremeters dependent on type of time dependence
  array[num_class, num_basis] real<lower=0> alpha;              // coefficients for spline
  array[num_class, num_serosurvey] real<lower=0, upper=1> pi_;

  //for quasi poisson model
  real<lower=1.5> theta;
}

transformed parameters {

  // fixed parameters
  real tau = 1./(generation_time*fraction_pre);
  real gamma = 1./(generation_time*(1-fraction_pre));

  real I0 = (1+I0_raw);

  // generic steps
  array[num_class, num_t] real asc_incidence;

  array[num_t] vector[num_eq] y;
  if (sampler == 0 ){
    y = ode_rk45_tol( seir_2d_spline,
       rep_vector(0.0,num_eq),           // initial values = 0 (handled within the ODE)
        t0,                               // initial time = 0
        ts,                               // evaluation times
        rtol, atol, max_num_steps,        // tolerances
        I0, knots, alpha, b_hat, order,  // parameters
        tau, gamma, contact, beta_fixed, popdist,     // data
        DIM                               // metadata
        );

    // extract, rescale and format prevalence, simulate data
    asc_incidence = get_incidence(y, DIM, pi_, t_survey_start, t_survey_end);
  } else if (sampler == 4 ){
    y = solve_ode_system_trapezoidal_2d( rep_vector(0.0,num_eq),
                                    ts,
                                    I0, knots, alpha, b_hat, order,
                                    tau, gamma, contact, beta_fixed,
                                    popdist,
                                    DIM
                                    );
    // extract, rescale and format prevalence, simulate data
    asc_incidence = get_incidence(y, DIM, pi_, t_survey_start, t_survey_end);
  }

  // given the SEIR dynamics, calculate the probability to observe a seropositive individual in the population
  array[num_class, num_serosurvey] real cum_inf_frac;
  for (i in 1:num_class){
    for (q in 1:num_serosurvey){
      cum_inf_frac[i, q] = mean(to_vector(y[t_survey_start[q]:t_survey_end[q], ind(4,i, num_class) ])) / popdist[i];
    }
  }

}

model {
  // Priors
  I0_raw ~ gamma(p_I0[1]^2/p_I0[2]^2,p_I0[1]/p_I0[2]^2);
  //beta ~ normal(p_beta[1],p_beta[2]);

  for (i in 1:num_class){
    alpha[i][] ~ normal( (p_beta[1]), 0.1 ); // better informative prior --> update to more stable method for overfitting see ref. manual
  }

  theta ~ normal( p_theta[1], p_theta[2]);
  //fraction_pre ~ uniform(0.4, 0.6);
  pi_ ~ normal(0.5,0.1);

  // likelihood
  if(inference==1){
    for (j in 1:num_class){
      for (i in 1:num_serosurvey){
        n_infected_survey[j,i] ~ binomial(n_tested_survey[j, i], cum_inf_frac[j,i]*sens[i] + (1-cum_inf_frac[j,i])*(1-spec[i]));
      }
      for (i in 1:num_t){
        target += neg_binomial_2_log_lpmf( data_pre[j,i] | log(asc_incidence[j,i]), asc_incidence[j,i]/(theta-1) );
      }
    }
  }

}

generated quantities {

  array[num_class] vector[num_t] confirmed_cases_predicted;

  for (i in 1:num_class){
    confirmed_cases_predicted[i][] = to_vector(qpoisson_rng(num_t, to_vector(asc_incidence[i,]), theta));
  }

  // calculate transmission probability per age group
  array[num_class] vector[num_t] rho;
  for ( j in 1:num_class){
   for (i in 1:num_t){
      rho[j][i] = ( analytical_bspline( b_hat, to_vector(alpha[j,]), ts[i], knots, order ));
   }
  }

}


