//
// Basic transmission model with stratification and time dependence in the form of splines
// in the variables.

functions {
  #include "/functions/generic_functions.stan"
  #include "/functions/spline_functions.stan"
  #include "/functions/SEIR_ode_spline.stan"
  #include "/functions/prevalence_GE.stan" // prevalence function without stratification
}

// load data objects
data {
  // load basic data objects that are needed independent of the exact specifications of the model
  // data objects that are generic independent of the model specifications
  int num_t;
  int num_class;
  array[num_t] real ts;

  // priors
  vector[2] p_I0; // expected initial seed (mean, sd)
  vector[2] p_beta; // expected beta (alpha, beta)
  vector[2] p_theta;

  // fixed quantities
  real generation_time;
  real beta_fixed;
  real sens;
  real spec;
  array[2] int t_detectionSwitch;
  real fraction_pre;

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

  // First serosurvey
  int t_survey_start_1;
  int t_survey_end_1;
  array[num_class] int n_infected_survey_1;
  array[num_class] int n_tested_survey_1;

  // Second serosurvey
  int t_survey_start_2;
  int t_survey_end_2;
  array[num_class] int n_infected_survey_2;
  array[num_class] int n_tested_survey_2;

  // Third serosurvey
  int t_survey_start_3;
  int t_survey_end_3;
  array[num_class] int n_infected_survey_3;
  array[num_class] int n_tested_survey_3;

  // load data element that defines what sampler should be used
  int sampler;

}

transformed data {
  #include "/data/generic_transformed_data.stan"
  #include "/data/splines_transformed_data.stan"

  // specific parameters for stratified version
  int num_eq = num_comp*num_class;

  array[7] int DIM = {num_comp, num_t, num_prev, num_age, num_sex, num_class, num_eq};

}

parameters {
  // generic parameters
  array[num_class] real<lower=0, upper=1> beta;               // transmission probability per contact
  real<lower=0> I0_raw;                 // initial seed (in number of individuals)

  // paremeters dependent on type of time dependence
  array[num_class, num_basis-1] real<lower=0> a_raw;              // coefficients for spline
  //real<lower=0.4, upper=0.6> fraction_pre; // fraction of the generation time spend in compartment E
  // real<lower=0, upper=1> p_detect1; // ascertainmentrate first wave
  array[num_class] real<lower=0, upper=1> p_detect1; // ascertainmentrate second wave
  array[num_class] real<lower=0, upper=1> p_detect2; // ascertainmentrate second wave
  array[num_class] real<lower=0, upper=1> p_detect3; // ascertainmentrate second wave

  //for quasi poisson model
  real<lower=1.5> theta;
}

transformed parameters {

  // fixed parameters
  real tau = 1./(generation_time*fraction_pre);
  real gamma = 1./(generation_time*(1-fraction_pre));

  real I0 = (1+I0_raw);

  // generic steps
  array[num_class, num_t] real incidence;

  array[num_class, num_basis] real a;

  for (q in 1:num_class){
    a[q,1] = (beta[q])/beta_fixed; // can also be negative, as we use inv.logit insdide ODE
    for (i in 2:num_basis)
      a[q,i] =  a_raw[q, i-1];      // a[i-1] + a_raw[i-1]*kappa;
  }

  array[num_t] vector[num_eq] y;

  y = ode_rk45_tol( seir_2d_spline,
      rep_vector(0.0,num_eq),           // initial values = 0 (handled within the ODE)
      t0,                               // initial time = 0
      ts,                               // evaluation times
      rtol, atol, max_num_steps,        // tolerances
      I0, knots, a, b_hat, order,  // parameters
      tau, gamma, contact, beta_fixed, popdist,     // data
      DIM                               // metadata
      );

  // extract, rescale and format prevalence, simulate data
  incidence = get_incidence(y, DIM, p_detect1, p_detect2, p_detect3, t_detectionSwitch );

  array[num_class] real p_infected_survey_1;
  array[num_class] real p_infected_survey_2;
  array[num_class] real p_infected_survey_3;

  for (i in 1:num_class){
    p_infected_survey_1[i] = mean(to_vector(y[t_survey_start_1:t_survey_end_1, ind(4,i, num_class) ])) / popdist[i];
    p_infected_survey_2[i] = mean(to_vector(y[t_survey_start_2:t_survey_end_2, ind(4,i, num_class) ])) / popdist[i];
    p_infected_survey_3[i] = mean(to_vector(y[t_survey_start_3:t_survey_end_3, ind(4,i, num_class) ])) / popdist[i];
  }

}

model {
  // Priors
  I0_raw ~ gamma(p_I0[1]^2/p_I0[2]^2,p_I0[1]/p_I0[2]^2);
  beta ~ normal(p_beta[1],p_beta[2]);

  for (i in 1:num_class){
    a_raw[i][] ~ normal( (p_beta[1]), 0.1 ); // better informative prior --> update to more stable method for overfitting see ref. manual
  }

  theta ~ normal( p_theta[1], p_theta[2]);
  //fraction_pre ~ uniform(0.4, 0.6);
  p_detect1 ~ normal(0.1,0.1);
  p_detect2 ~ normal(0.5,0.1);
  p_detect3 ~ normal(0.5,0.1);

  // likelihood
  if(inference==1){
    for (j in 1:num_class){
      n_infected_survey_1[j] ~ binomial(n_tested_survey_1[j], p_infected_survey_1[j]*sens + (1-p_infected_survey_1[j])*(1-spec));
      n_infected_survey_2[j] ~ binomial(n_tested_survey_2[j], p_infected_survey_2[j]*sens + (1-p_infected_survey_2[j])*(1-spec));
      n_infected_survey_3[j] ~ binomial(n_tested_survey_3[j], p_infected_survey_3[j]*sens + (1-p_infected_survey_3[j])*(1-spec));
      for (i in 1:num_t){
        target += neg_binomial_2_log_lpmf( data_pre[j,i] | log(incidence[j,i]), incidence[j,i]/(theta-1) );
      }
    }
  }

}

generated quantities {

  array[num_class] vector[num_t] I_t_predicted;

  for (i in 1:num_class){
    I_t_predicted[i][] = to_vector(qpoisson_rng(num_t, to_row_vector(incidence[i,]), theta));
  }

  // calculate transmission probability per age group
  array[num_class] vector[num_t] prob_infection;
  for ( j in 1:num_class){
   for (i in 1:num_t){
      prob_infection[j][i] = ( analytical_bspline( b_hat, to_vector(a[j,]), ts[i], knots, order ));
   }
  }

  // cumulative incidence
  vector[num_class] cum_inc;

  for (j in 1:num_class){
    cum_inc[j] = sum( incidence[j,] );
  }

}


