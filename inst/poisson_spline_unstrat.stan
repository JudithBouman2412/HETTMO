//
// Basic transmission model with stratification and time dependence in the form of splines
// in the variables.

functions {
    // Generic functions
    int[] sequence(int from, int to){
      int s = to - from + 1;
      int v[s];
      for(i in 1:s){
        v[i] = from + i - 1;
      }
      return v;
  }

  // get the vector index for a given compartment index and age_sex index
  int ind(int comp, int age_sex, int num_class){
    return(age_sex+(comp-1)*num_class);
  }

  // get the vector index for several compartment indexes and age_sex indexes
  int[] mind(int[] comp, int[] age_sex, int num_class){
    int s1=size(comp);
    int s2=size(age_sex);
    int v[s1*s2];
    int index=1;
    for(i in comp){
      for(j in age_sex){
        v[index]=j+(i-1)*num_class;
        index=index+1;
      }
    }
    return v;
  }

  // quasi poisson, random number generator
  int[] qpoisson_rng(int n, row_vector mu, real theta  ){
    int out[n] = neg_binomial_2_rng( mu, to_array_1d( mu/(theta-1) ) );
    return(out);
  }


  // infer value at time t between two time points using logistic switch
  real between_points( vector beta_weekly, real t, data int num_t, real[] ts ){

    int i = 1;
    real tdif = fabs(t - i);
    real k = 100; // slope
    real lowerbnd;
    real upperbnd;
    real approx_t;

    while (tdif>1 && i<(size(ts)-1)){
      i = i + 1;
      tdif = fabs(t-i);
    }

    lowerbnd = beta_weekly[i];
    upperbnd = beta_weekly[i+1];

    approx_t = lowerbnd + (upperbnd-lowerbnd)/(1+exp(-k*(t-(ts[i+1]+ts[i])/2)));

    return approx_t;
  }

  // infer value at time t between two time points using logistic switch
  real[] between_points_class( real[,] beta_weekly, real t, data int num_t, real[] ts ){

    int i = 1;
    real tdif = fabs(t - i);
    real k = 100; // slope
    int num_class = dims(beta_weekly)[1];
    vector[num_class] lowerbnd;
    vector[num_class] upperbnd;
    real approx_t[num_class];

    while (tdif>1 && i<(size(ts)-1)){
      i = i + 1;
      tdif = fabs(t-i);
    }

    lowerbnd = to_vector( beta_weekly[,i]);
    upperbnd = to_vector( beta_weekly[,i+1]);

    approx_t = to_array_1d(lowerbnd + (upperbnd-lowerbnd)/(1+exp(-k*(t-(ts[i+1]+ts[i])/2))));

    return approx_t;
  }

    // spline functions

  vector build_b_spline(real[] t, data vector ext_knots, data int index, data int order);
  vector build_b_spline(real[] t, data vector ext_knots, data int index, data int order){
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    index:        the indexex of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[index] <= t[i]) && (t[i] < ext_knots[index+1]);
    else {
      if (ext_knots[index] != ext_knots[index+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[index], size(t))) /
             (ext_knots[index+order-1] - ext_knots[index]);
      if (ext_knots[index+1] != ext_knots[index+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[index+1], size(t))) /
                 (ext_knots[index+order] - ext_knots[index+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, index, order-1) +
                 w2 .* build_b_spline(t, ext_knots, index+1, order-1);
    }
    return b_spline;
  }

  matrix find_support(vector ext_knots, data int order, vector knots, data int num_basis) {
    int num_interval = dims(knots)[1] - 1;
    // define t for which polynomial is well defined for each B-spline
    matrix[num_basis, num_interval] support_bs;

    for (i in 1:num_basis) {
      vector[num_interval] support = rep_vector(0, num_interval);
      for (j in 1:num_interval) {
        if (knots[j] >= ext_knots[i]) {
          if (knots[j + 1] <= ext_knots[i + order]) {
            support[j] = 1;
          }
        }
      }
      support_bs[i,] = to_row_vector(support);
    }

    return support_bs;
  }

  vector find_poly( data real[] x, data row_vector y){
    // given m points x,y find polynomial of degree m-1 through these points
    int m = size(x);

    // vector initializing the polynomial
    vector[m] thepoly = rep_vector(0, m);

    for (i in 1:m){
      vector[m] theterm = rep_vector(0, m);

      real r = 1;
      for (j in 1:m){
        if (i!=j){
          r = r * (x[i]-x[j]);
        }
      }

      theterm[m] = y[i]/r;

      for (j in 1:m){
        if (i!=j){
          for (k in 2:m){
            theterm[k-1] += theterm[k];
            theterm[k] *= (-x[j]);
          }
        }
      }

      for (j in 1:m){
        thepoly[j] = thepoly[j]+theterm[j];
      }
    }

    thepoly = rev(thepoly);

    return thepoly;
  }

  matrix ana_B( data matrix B, data vector ext_knots, data array[] real int_knots, data int order, data vector knots, data int num_basis, data matrix support_bs){
      // empty matrix to save estimated coefficients for each b_spline and each area at which it is defined
      int num_interval = size(knots)-1;
      matrix[num_basis, num_interval*order] b_hat = rep_matrix(0.0, num_basis, num_interval*order);

      // fit polynomal for each b_spline at the interval at which it is defined
      for (i in 1:num_basis){
        for (j in 1:num_interval){
          if (support_bs[i,j]){
            b_hat[i,((j-1)*order+1):(j*order)] = to_row_vector( find_poly(int_knots[((j-1)*(order-1)+1):((order-1)*j+1)], B[i,((j-1)*(order-1)+1):((order-1)*j+1)]) ); //extend to use any degree/order
          }
        }
      }

      return b_hat;
    }

  // function to calculate full spline for new values
  real analytical_bspline( data matrix b_hat, vector a, real t_new, data vector knots, data int order ) {

    int num_interval = size(knots)-1;
    int n_b_spline = dims(b_hat)[1];

    matrix[n_b_spline, order] b_hat_now;

    for (i in 1:num_interval){
      if (t_new>=knots[i]){
        if (t_new<knots[i+1]){
          b_hat_now = b_hat[,((i-1)*order+1):(i*order)];
        }
      }
    }

    vector[n_b_spline] b_spline;

    for (i in 1:n_b_spline){
      b_spline[i] = a[i] * (b_hat_now[i,1] +b_hat_now[i,2]*t_new +b_hat_now[i,3]*t_new^2 + b_hat_now[i,4]*t_new^3);
    }

    real b_spline_1 = sum(b_spline);

    return(b_spline_1);
  }

    // ode spline functions
  // ODE for basic SEIR model with spline function for time dependence in beta and no stratification
  vector seir_0d_spline(real t,
                        vector y,
                        real I0,
                        data vector knots,
                        vector a,
                        data matrix b_hat,
                        data int order,
                        real tau,
                        real gamma,
                        data real contact,
                        data real beta_fixed,
                        data array[] int DIM) {

    // array with all necessary dimensions
    int num_comp = DIM[1];
    int popsize = DIM[4];
    // define compartments
    vector[num_comp] dydt;
    // define initial values
    vector[num_comp] init = to_vector(append_array( {(popsize-I0), 0, (I0)}, rep_array(0, num_comp-3))); // define time-varying force of infection, by creating splines
    real beta_spline = analytical_bspline(b_hat, a, t, knots, order); // logit transformation for beta_spline?
    // define ODEs
    real finf_t = beta_fixed * (beta_spline) * contact * (y[3]+init[3]);
    // define ODEs
    dydt[1] = -finf_t * (y[1]+init[1])/popsize;
    dydt[2] =  finf_t * (y[1]+init[1])/popsize - tau * (y[2]+init[2]);
    dydt[3] =  tau * (y[2]+init[2]) - gamma * (y[3]+init[3]);
    dydt[4] =  gamma * (y[3]+init[3]);

    return dydt;
  }

  // ODE for basic SEIR model with spline function for time dependence in beta and no stratification
  vector seir_0d_spline_logtrans(real t, vector y, real I0, data vector knots, vector a,
                        data matrix b_hat, data int order, real tau, real gamma,
                        data real contact, data real beta_fixed, data array[] int DIM) {

    // array with all necessary dimensions
    int num_comp = DIM[1];
    int popsize = DIM[4];
    // define compartments
    vector[num_comp] dydt;

    // define initial values
    vector[num_comp] init = to_vector(append_array( {log(1-I0), log(1e-8), log(I0)}, rep_array(log(1e-8), num_comp-3))); // define time-varying force of infection, by creating splines
    //print(init);
    real beta_spline = analytical_bspline(b_hat, a, t, knots, order); // logit transformation for beta_spline?
    // define ODEs
    dydt[1] = - beta_fixed * (beta_spline)*contact*(exp(y[3]+init[3]));
    dydt[2] =  beta_fixed * ((beta_spline)*contact*(exp(y[3]+init[3]))*(exp(y[1]+init[1])))/(exp(y[2]+init[2])) - tau;
    dydt[3] =  tau*(exp(y[2]+init[2]))/(exp(y[3]+init[3])) - gamma;
    dydt[4] =  gamma*(exp(y[3]+init[3]))/(exp(y[4]+init[4]));

    return dydt;
  }


  // ODE for stratified SEIR model without time dependence in beta
  vector seir_2d_spline(real t,
                        vector y,
                        real I0,
                        data vector knots,
                        array[,] real a,
                        data matrix b_hat,
                        data int order,
                        real tau,
                        real gamma,
                        data matrix contact,
                        data real beta_fixed,
                        data vector popdist,
                        data array[] int DIM) {
    // dimensions
    int num_comp = DIM[1];
    int num_class = DIM[6];
    int num_eq = DIM[7];

    // define equations
    vector[num_eq] dydt;

    // define initial values
    vector[num_class] initI = to_vector([ I0 , I0 , I0]);
    vector[num_class] initS = popdist - I0;

    // define time-varying force of infection
    vector[num_class] beta_spline;
    for (i in 1:num_class){
      beta_spline[i] = analytical_bspline( b_hat, to_vector(a[i,]), t, knots, order ); // logit transformation for beta_spline?
    }

    vector[num_class] f_inf_t = beta_fixed * contact * ((to_vector(y[mind({3},seq(1,num_class), num_class )]) + to_vector(initI)) ./ to_vector(popdist));

     // define ODEs
    for(i in 1:num_class){
      dydt[ind(1,i, num_class)] = - (beta_spline[i]) * f_inf_t[i] * (y[ind(1,i, num_class)]+initS[i]);
      dydt[ind(2,i, num_class)] =  (beta_spline[i]) * f_inf_t[i] * (y[ind(1,i, num_class)]+initS[i]) - tau * y[ind(2,i, num_class)];
      dydt[ind(3,i, num_class)] =  tau * y[ind(2,i, num_class)] - gamma * (y[ind(3,i, num_class)]+initI[i]);
      dydt[ind(4,i, num_class)] =  gamma * (y[ind(3,i, num_class)]+initI[i]);
    }
    return dydt;
  }

  // function to approximate the solution of an ODE function using the Trapeziodal rule
  array[] vector solve_ode_system_trapezoidal( vector initial_state,
                                             array[] real ts,
                                             real I0,
                                             data vector knots,
                                             vector a,
                                             data matrix b_hat,
                                             data int order,
                                             real tau,
                                             real gamma,
                                             data real contact,
                                             data real beta_fixed,
                                             data array[] int DIM
                                             ) {

    int num_comp = DIM[1];
    int popsize = DIM[4];
    int N_hsub = DIM[7];

    // Initialize the output vector to store the solution
    array[size(ts)] vector[num_comp] ys;
    ys[1] = initial_state;
    vector[num_comp] ytemp;

    // Loop through each time point and solve the ODE system using the trapezoidal rule
    for (i in 2:size(ts)) {

      real h = ts[i] - ts[i-1];

      real stepsize_sub = h/N_hsub;
      array[N_hsub+1] real ts_sub; // divide the time steps between ts1 and ts2 in N_hsub steps
      ts_sub[1] = ts[i-1];
      array[N_hsub+1] vector[num_comp] ys_sub; // create array of vectors to save the sub steps
      ys_sub[1] = ys[i-1];

      for (q in 1:N_hsub) {
        ytemp = ys_sub[q];

        ts_sub[q+1] = ts_sub[q] + stepsize_sub;

        vector[num_comp] k1 = seir_0d_spline(ts_sub[q], ytemp, I0, knots, a, b_hat, order, tau, gamma, contact, beta_fixed, DIM);

        for (j in 1:num_comp) ytemp[j] = ytemp[j] + stepsize_sub * k1[j];

        vector[num_comp] k2 = seir_0d_spline(ts_sub[q+1], ytemp, I0, knots, a, b_hat, order, tau, gamma, contact, beta_fixed, DIM);

        for (j in 1:num_comp) ys_sub[q+1][j] = ys_sub[q][j] + stepsize_sub * (k1[j] + k2[j]) / 2;
      }

      vector[num_comp] k3 = seir_0d_spline(ts_sub[N_hsub], ytemp, I0, knots, a, b_hat, order, tau, gamma, contact, beta_fixed, DIM);

      for (j in 1:num_comp) ytemp[j] = ytemp[j] + stepsize_sub * k3[j];

      vector[num_comp] k4 = seir_0d_spline(ts[i], ytemp, I0, knots, a, b_hat, order, tau, gamma, contact, beta_fixed, DIM);

      for (j in 1:num_comp) ys[i][j] = ys_sub[N_hsub][j] + stepsize_sub * (k3[j] + k4[j]) / 2;

    }

    return ys;
  }

    // prevalence functions without stratification

  // prevalence function without stratification

  // Get prevalence from result of ODE of SEIR model
  array[] real get_prevalence(array[] vector y, data array[] int DIM, data real popsize, data real atol, real I0) {
    // dimensions
    int num_comp = DIM[1];
    int num_t = DIM[2];
    int num_prev = DIM[3];
    // extract, rescale and format prevalence
    array[num_t] real prev = to_array_1d( ( to_vector(y[,num_prev]) + I0 + 2*atol ) * popsize );
    return prev;
  }


  row_vector get_incidence(array[] vector y, data array[] int DIM, data real popsize, data real atol, real I0,
                            real p_detect1, real p_detect2, data real t_detectionSwitch ) {
    // incidence defined as people moving from I to R
    // dimensions
    int num_comp = DIM[1];
    int num_t = DIM[2];
    int num_prev = DIM[3];

    // extract, rescale and format prevalence
    row_vector[num_t] inc;
    inc[1] = p_detect1 *(y[1, num_prev+1] + 1);

    for (i in 2:num_t){ // maybe do a vector multiplication instead to make it faster?
      if (i < t_detectionSwitch){
        inc[i] = p_detect1 * (y[i,num_prev+1] - y[i-1,num_prev+1] + 1);
      } else {
        inc[i] = p_detect2 * (y[i,num_prev+1] - y[i-1,num_prev+1] + 1);
      }
    }
    return inc;
  }

  array[] real get_LOGprevalence(array[] vector y, data array[] int DIM, data real popsize, data real atol, real I0) {

    // dimensions
    int num_t = DIM[2];
    int num_prev = DIM[3];

    // extract, rescale and format prevalence
    array[num_t] real prev;
    prev = to_array_1d( exp( to_vector( y[,num_prev]) + log(I0) ) * popsize );

    return prev;
  }

  array[] real get_LOGincidence(array[] vector y, data array[] int DIM, data real popsize, data real atol, real I0) {
    // dimensions
    int num_comp = DIM[1];
    int num_t = DIM[2];
    int num_prev = DIM[3];

    // extract, rescale and format prevalence
    array[num_t] real inc;
    inc[1] =  ( exp(y[1,num_prev] + log(I0)) + exp(y[1,num_prev+1] + log(I0)) + 2*atol ) * popsize ;
    for (i in 2:num_t){
      inc[i] = ( exp(y[i,num_prev]+ log(I0)) - exp(y[i-1,num_prev]+ log(I0)) + exp(y[i,num_prev+1]-10) - exp(y[i-1,num_prev+1]-10) + 2*atol ) * popsize;
    }

    return inc;
  }

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

  // specific for splines
  // transformed data block for splines: create b_hat
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  int order = spline_degree + 1;

  // create extended knots for fitting the splineC
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));

    // Create internal knots to have exactly enough points for each interval to fit a polynomial of degree 3
  int num_int = (num_knots-1)*(order-1)+1;
  array[num_int] real int_knots;
  array[order-1] real int_now;
  for (i in 2:num_knots){
    real dif = floor( (knots[i]-knots[i-1])/(order-1));
    for (j in 1:(order-1)){
      int_now[j] = knots[i-1] + (j-1)*dif;
    }
    int_knots[((i-2)*(order-1)+1):((i-1)*(order-1))] = int_now;
  }
  int_knots[num_int] = knots[num_knots];

  // Build the spline --> this spline stays the same now
  matrix[num_basis, num_int] B;

  for (indicator in 1:num_basis){
    B[indicator,:] = to_row_vector(build_b_spline( int_knots, ext_knots, indicator, spline_degree + 1));
  }
  B[num_knots + spline_degree - 1, num_int] = 1; // --> this is in the manual, but I do not understand why...

  // find support for each spline
  matrix[num_basis, (num_knots-1) ] support_bs = find_support( ext_knots, order, knots, num_basis );

  // find coefficients for each support of each spline
  matrix[num_basis, (num_knots-1)*order] b_hat = ana_B( B, ext_knots, int_knots, order, knots, num_basis, support_bs);


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

}

transformed parameters {
  // generic steps
  real tau = 1./(generation_time*fraction_pre);
  real gamma = 1./(generation_time*(1-fraction_pre));

  real I0 = (1+I0_raw);

  // calculate beta from R0
  real beta = ((R0*gamma)/contact)/beta_fixed;

  //real phi = 1./phi_inv;
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

  // poisson model
  if (inference==1) {
    n_infected_survey ~ binomial(n_tested_survey,p_infected_survey*sens + (1-p_infected_survey)*(1-spec));
    target += poisson_lpmf( data_pre | incidence );
  }

}

generated quantities {

  // log-likelihood for LOO
  vector[num_t] log_lik;

  //  poisson
  for (t in 1:num_t)
      log_lik[t] = poisson_lpmf( data_pre[t] | incidence[t] );
      array[num_t] int I_t_simulated = data_pre;

  array[num_t] int I_t_predicted = poisson_rng( incidence );

  // calculate probability of infection
  vector[num_t] prob_infection;
  for (i in 1:num_t){
    prob_infection[i] = inv_logit( analytical_bspline( b_hat, a, ts[i], knots, order ));
  }

}


