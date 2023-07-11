// ODE for basic SEIR model with brownian motion for time dependence in beta and no stratification
vector seir_0d_BM(real t,
                      vector y,
                      real I0,
                      vector beta_weekly,
                      real tau,
                      real gamma,
                      data real contact,
                      data real beta_fixed,
                      data array[] real ts,
                      data array[] int DIM) {

  // array with all necessary dimensions
  int num_comp = DIM[1];
  int num_t = DIM[2];
  int popsize = DIM[4];
  // define compartments
  vector[num_comp] dydt;
  // define initial values
  vector[num_comp] init = to_vector(append_array( {(popsize-I0), 0, (I0)}, rep_array(0, num_comp-3))); // define time-varying force of infection, by creating splines
  real beta_BM = between_points( beta_weekly, t, num_t, ts );
  // define ODEs
  real finf_t = beta_fixed * beta_BM * contact * (y[3]+init[3]);
  // define ODEs
  dydt[1] = -finf_t * (y[1]+init[1])/popsize;
  dydt[2] =  finf_t * (y[1]+init[1])/popsize - tau * (y[2]+init[2]);
  dydt[3] =  tau * (y[2]+init[2]) - gamma * (y[3]+init[3]);
  dydt[4] =  gamma * (y[3]+init[3]);

  return dydt;
}

// ODE for stratified SEIR model without time dependence in beta
vector seir_2d_BM(real t,
                      vector y,
                      real I0,
                      real tau,
                      array[,] real beta_weekly,
                      real gamma,
                      data matrix contact,
                      data real beta_fixed,
                      data vector popdist,
                      data array[] real ts,
                      data array[] int DIM) {
  // dimensions
  int num_comp = DIM[1];
  int num_class = DIM[6];
  int num_eq = DIM[7];
  int num_t = DIM[2];

  // define equations
  vector[num_eq] dydt;

  // define initial values
  vector[num_class] initI = to_vector([ I0 , I0 , I0]);
  vector[num_class] initS = popdist - I0;

  // beta given by the BM process
  array[num_class] real beta_BM = between_points_class( beta_weekly, t, num_t, ts );

  // define time-varying force of infection
  vector[num_class] f_inf_t = beta_fixed * contact * ((to_vector(y[mind({3},seq(1,num_class), num_class )]) + to_vector(initI)) ./ to_vector(popdist));

   // define ODEs
  for(i in 1:num_class){
    dydt[ind(1,i, num_class)] = - (beta_BM[i]) * f_inf_t[i] * (y[ind(1,i, num_class)]+initS[i]);
    dydt[ind(2,i, num_class)] =  (beta_BM[i]) * f_inf_t[i] * (y[ind(1,i, num_class)]+initS[i]) - tau * y[ind(2,i, num_class)];
    dydt[ind(3,i, num_class)] =  tau * y[ind(2,i, num_class)] - gamma * (y[ind(3,i, num_class)]+initI[i]);
    dydt[ind(4,i, num_class)] =  gamma * (y[ind(3,i, num_class)]+initI[i]);
  }
  return dydt;
}

// ODE for basic SEIR model with Brownian motion for time dependence in beta and no stratification
// special function with stepwise difference in transmission, only to evaluate with trapeziodal solver
vector seir_0d_BM_trap(real t,
                      vector y,
                      real I0,
                      real beta_now,
                      real tau,
                      real gamma,
                      data real contact,
                      data real beta_fixed,
                      data array[] int DIM) {

  // array with all necessary dimensions
  int num_comp = DIM[1];
  int num_t = DIM[2];
  int popsize = DIM[4];
  // define compartments
  vector[num_comp] dydt;
  // define initial values
  vector[num_comp] init = to_vector(append_array( {(popsize-I0), 0, (I0)}, rep_array(0, num_comp-3))); // define time-varying force of infection, by creating splines
  // define ODEs
  real finf_t = beta_fixed * beta_now * contact * (y[3]+init[3]);
  // define ODEs
  dydt[1] = -finf_t * (y[1]+init[1])/popsize;
  dydt[2] =  finf_t * (y[1]+init[1])/popsize - tau * (y[2]+init[2]);
  dydt[3] =  tau * (y[2]+init[2]) - gamma * (y[3]+init[3]);
  dydt[4] =  gamma * (y[3]+init[3]);

  return dydt;
}

vector seir_2d_BM_trap(real t,
                      vector y,
                      real I0,
                      array[,] real beta_weekly,
                      real tau,
                      real gamma,
                      matrix contact,
                      real beta_fixed,
                      vector popdist,
                      array[] real ts,
                      array[] int DIM) {

  // dimensions
  int num_comp = DIM[1];
  int num_class = DIM[6];
  int num_eq = DIM[7];
  int num_t = DIM[2];

  // define equations
  vector[num_eq] dydt;

  // define initial values
  vector[num_class] initI = to_vector([ I0 , I0 , I0]);
  vector[num_class] initS = popdist - I0;

  // beta given by the BM process
  array[num_class] real beta_BM = between_points_class( beta_weekly, t, num_t, ts );

  // define time-varying force of infection
  vector[num_class] f_inf_t = beta_fixed * contact * ((to_vector(y[mind({3},seq(1,num_class), num_class )]) + to_vector(initI)) ./ to_vector(popdist));

   // define ODEs
  for(i in 1:num_class){
    dydt[ind(1,i, num_class)] = - (beta_BM[i]) * f_inf_t[i] * (y[ind(1,i, num_class)]+initS[i]);
    dydt[ind(2,i, num_class)] =  (beta_BM[i]) * f_inf_t[i] * (y[ind(1,i, num_class)]+initS[i]) - tau * y[ind(2,i, num_class)];
    dydt[ind(3,i, num_class)] =  tau * y[ind(2,i, num_class)] - gamma * (y[ind(3,i, num_class)]+initI[i]);
    dydt[ind(4,i, num_class)] =  gamma * (y[ind(3,i, num_class)]+initI[i]);
  }
  return dydt;
}

// function to approximate the solution of an ODE function using the Trapeziodal rule
array[] vector solve_ode_system_trapezoidal( vector initial_state,
                                           array[] real ts,
                                           real I0,
                                           vector beta_weekly,
                                           real tau,
                                           real gamma,
                                           data real contact,
                                           data real beta_fixed,
                                           data array[] int DIM
                                           ) {

  int num_comp = DIM[1];
  int popsize = DIM[4];
  int N_hsub = DIM[5];

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

      vector[num_comp] k1 = seir_0d_BM_trap(ts_sub[q], ytemp, I0, beta_weekly[i-1], tau, gamma, contact, beta_fixed, DIM);

      for (j in 1:num_comp) ytemp[j] = ytemp[j] + stepsize_sub * k1[j];

      vector[num_comp] k2 = seir_0d_BM_trap(ts_sub[q+1], ytemp, I0, beta_weekly[i-1], tau, gamma, contact, beta_fixed, DIM);

      for (j in 1:num_comp) ys_sub[q+1][j] = ys_sub[q][j] + stepsize_sub * (k1[j] + k2[j]) / 2;
    }

    vector[num_comp] k3 = seir_0d_BM_trap(ts_sub[N_hsub], ytemp, I0, beta_weekly[i-1], tau, gamma, contact, beta_fixed, DIM);

    for (j in 1:num_comp) ytemp[j] = ytemp[j] + stepsize_sub * k3[j];

    vector[num_comp] k4 = seir_0d_BM_trap(ts[i], ytemp, I0, beta_weekly[i-1], tau, gamma, contact, beta_fixed, DIM);

    for (j in 1:num_comp) ys[i][j] = ys_sub[N_hsub][j] + stepsize_sub * (k3[j] + k4[j]) / 2;

  }

  return ys;
}


// function to approximate the solution of an ODE function using the Trapeziodal rule
array[] vector solve_ode_system_trapezoidal_2D( vector initial_state,
                                           array[] real ts,
                                           real I0,
                                           real tau,
                                           array[,] real beta_weekly,
                                           real gamma,
                                           matrix contact,
                                           real beta_fixed,
                                           vector popdist,
                                           array[] int DIM
                                           ) {

  int num_comp = DIM[1];
  int popsize = DIM[4];
  int N_hsub = DIM[5];
  int num_class = size(popdist);
  int num_eq = DIM[7];

  // Initialize the output vector to store the solution
  array[size(ts)] vector[num_eq] ys;
  ys[1][] = initial_state;
  vector[num_eq] ytemp;

  // Loop through each time point and solve the ODE system using the trapezoidal rule
  for (i in 2:size(ts)) {

    real h = ts[i] - ts[i-1];

    real stepsize_sub = h/N_hsub;
    array[N_hsub+1] real ts_sub; // divide the time steps between ts1 and ts2 in N_hsub steps
    ts_sub[1] = ts[i-1];
    array[N_hsub+1] vector[num_eq] ys_sub; // create array of vectors to save the sub steps
    ys_sub[1][] = ys[i-1][];

    for (q in 1:N_hsub) {
      ytemp = ys_sub[q][];

      ts_sub[q+1] = ts_sub[q] + stepsize_sub;

      vector[num_eq] k1 = seir_2d_BM_trap(ts_sub[q], ytemp, I0, beta_weekly, tau, gamma, contact, beta_fixed, popdist,ts,  DIM);

      for (j in 1:num_eq) ytemp[j] = ytemp[j] + stepsize_sub * k1[j];

      vector[num_eq] k2 = seir_2d_BM_trap(ts_sub[q+1], ytemp, I0, beta_weekly, tau, gamma, contact, beta_fixed, popdist,ts,  DIM);

      for (j in 1:num_eq) ys_sub[q+1][j] = ys_sub[q][j] + stepsize_sub * (k1[j] + k2[j]) / 2;
    }

    vector[num_eq] k3 = seir_2d_BM_trap(ts_sub[N_hsub], ytemp, I0, beta_weekly, tau, gamma, contact, beta_fixed, popdist, ts, DIM);

    for (j in 1:num_eq) ytemp[j] = ytemp[j] + stepsize_sub * k3[j];

    vector[num_eq] k4 = seir_2d_BM_trap(ts[i], ytemp, I0, beta_weekly, tau, gamma, contact, beta_fixed, popdist, ts, DIM);

    for (j in 1:num_eq) ys[i][j] = ys_sub[N_hsub][j] + stepsize_sub * (k3[j] + k4[j]) / 2;

  }

  return ys;
}


