// basis function (exponentiated quadratic kernel)
matrix PHI_EQ(data int N, data int M, data real L, data array[] real x) {
  matrix[N,M] A = rep_matrix(pi()/(2*L) * (to_vector(x)+L), M);
  vector[M] B = linspaced_vector(M, 1, M);
  matrix[N,M] PHI = sin(diag_post_multiply(A, B))/sqrt(L);
  //for (m in 1:M) PHI[,m] = PHI[,m] - mean(PHI[,m]); // remove demeaning, see https://discourse.mc-stan.org/t/hilbert-space-approximate-gp-prior-why-use-de-meaned-basis-functions-phi/25852
  return PHI;
}

// spectral density (exponentiated quadratic kernel)
vector diagSPD_EQ(real alpha, real lambda, data real L, data int M) {
  vector[M] B = linspaced_vector(M, 1, M);
  return sqrt( alpha^2 * sqrt(2*pi()) * lambda * exp(-0.5*(lambda*pi()/(2*L))^2*B^2) );
}


// Given PHI, create a function that returns a value any value of t
// for stepwise version of phi
// also make a version that creates a linear approximation between the point?
row_vector ana_phi( matrix PHI, data int M, real t, data int num_t, array[] real ts ){

  int i = 1;
  real tdif = abs(t - i);
  real k = 10; // slope
  int ts_length = size(ts);

  while (tdif>1 && i < (ts_length-1)){
    i = i + 1;
    tdif = abs(t-i);
  }

  row_vector[M] lowerbnd = PHI[i,];
  row_vector[M] upperbnd = PHI[i+1,];

  row_vector[M] PHI_approx_t = lowerbnd + (upperbnd-lowerbnd)/(1+exp(-k*(t-(ts[i+1]+ts[i])/2)));

  return PHI_approx_t;
}

// for the original GP method: find values of the function and find intermediate value
vector GP_original( real lambda, real alpha, vector beta_f1, data array[] real ts,
                    real sigma, data int num_t ){

  matrix[num_t, num_t] cov = cov_exp_quad( ts, alpha, lambda )  + diag_matrix( rep_vector( sigma, num_t )); // to stabilize calculations, do we need that?
  matrix[num_t, num_t] L_cov = cholesky_decompose(cov);

  // this:  vector[num_t] random_function = multi_normal_cholesky_rng( rep_vector(0, num_t), L_cov);
  // or this: ???
  vector[num_t] random_function = L_cov * beta_f1;

  return random_function;
}

real ana_random( vector random_function, real t, data array[] real ts ){

  int i = 1;
  real tdif = abs(t - i);
  real k = 100; // slope

  while (tdif>1){
    i = i + 1;
    tdif = abs(t-i);
  }

  real lowerbnd = random_function[i];
  real upperbnd = random_function[i+1];

  real approx_value = lowerbnd + (upperbnd-lowerbnd)/(1+exp(-k*(t-(ts[i+1]+ts[i])/2)));

  return approx_value;
}
