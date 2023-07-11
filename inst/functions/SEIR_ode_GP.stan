vector seir_0d_GP(real t, vector y, real I0, real beta, vector beta_f1, vector diagSPD_f1,  matrix PHI, real tau, real gamma,
                  data real contact, data real beta_fixed, array[] real ts, array[] int DIM, data int M_f1) {
  // array with all necessary dimensions
  int num_comp = DIM[1];
  int num_t = DIM[2];
  int popsize = DIM[5];

  // define compartments
  vector[num_comp] dydt;
  // define initial values
  vector[num_comp] init = to_vector(append_array({popsize-I0,0.0, I0}, rep_array(0.0, num_comp-3)));
  // define time-varying force of infection, by creating building block of GP
  row_vector[M_f1] PHI_t = ana_phi(PHI, M_f1, t, num_t, ts ); // scale ts such that l is correct? and then scale back...? Does that slow down when I do it every step?
  // compute f1
  real f1_t = PHI_t * (diagSPD_f1 .* beta_f1); // f1_t is scaled to have mean 0, multiply with intercept?=beta? take logit_inv as with spline?
  real finf_t = beta_fixed * ((beta) + inv_logit(f1_t)) * contact * (y[3]+init[3]); // times two to have mean at beta
  // define ODEs
  dydt[1] = -finf_t * (y[1]+init[1])/popsize;
  dydt[2] =  finf_t * (y[1]+init[1])/popsize - tau * (y[2]+init[2]);
  dydt[3] =  tau * (y[2]+init[2]) - gamma * (y[3]+init[3]);
  dydt[4] =  gamma * (y[3]+init[3]);
  return dydt;
}


// function to approximate the solution of an ODE function using the Trapeziodal rule
array[] vector solve_ode_system_trapezoidal( vector initial_state,
                                           array[] real ts,
                                           real I0,
                                           data real L_f1,
                                           data int M_f1,
                                           matrix PHI,
                                           real beta,
                                           vector beta_f1,
                                           vector diagSPD_f1,
                                           real tau,
                                           real gamma,
                                           data real contact,
                                           data real beta_fixed,
                                           data array[] int DIM
                                           ) {

  int num_comp = DIM[1];
  int popsize = DIM[4];
  int N_hsub = DIM[6];

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

      vector[num_comp] k1 = seir_0d_GP(ts_sub[q], ytemp, I0, beta, beta_f1, diagSPD_f1,
                              PHI, tau, gamma, contact, beta_fixed, ts, DIM, M_f1);

      for (j in 1:num_comp) ytemp[j] = ytemp[j] + stepsize_sub * k1[j];

      vector[num_comp] k2 =  seir_0d_GP(ts_sub[q+1], ytemp, I0, beta, beta_f1, diagSPD_f1, PHI, tau, gamma, contact, beta_fixed, ts, DIM, M_f1);

      for (j in 1:num_comp) ys_sub[q+1][j] = ys_sub[q][j] + stepsize_sub * (k1[j] + k2[j]) / 2;
    }

    vector[num_comp] k3 =  seir_0d_GP(ts_sub[N_hsub], ytemp, I0, beta, beta_f1, diagSPD_f1, PHI, tau, gamma, contact, beta_fixed, ts, DIM, M_f1);

    for (j in 1:num_comp) ytemp[j] = ytemp[j] + stepsize_sub * k3[j];

    vector[num_comp] k4 = seir_0d_GP(ts[i], ytemp, I0, beta, beta_f1, diagSPD_f1, PHI, tau, gamma, contact, beta_fixed, ts, DIM, M_f1);

    for (j in 1:num_comp) ys[i][j] = ys_sub[N_hsub][j] + stepsize_sub * (k3[j] + k4[j]) / 2;

  }

  return ys;
}

