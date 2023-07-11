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
