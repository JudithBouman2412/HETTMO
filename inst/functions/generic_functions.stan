array[] int seq(int from, int to){
  int s = to - from + 1;
  array[s] int v;
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
array[] int mind(array[] int comp, array[] int age_sex, int num_class){
  int s1=size(comp);
  int s2=size(age_sex);
  array[s1*s2] int v;
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
array[] int qpoisson_rng(int n, row_vector mu, real theta  ){
  array[n] int out = neg_binomial_2_rng( mu, to_array_1d( mu/(theta-1) ) );
  return(out);
}

// infer value at time t between two time points using logistic switch
real between_points( vector beta_weekly, real t, data int num_t, array[] real ts ){

  int i = 1;
  real tdif = abs(t - i);
  real k = 100; // slope

  while (tdif>1 && i<(size(ts)-1)){
    i = i + 1;
    tdif = abs(t-i);
  }

  real lowerbnd = beta_weekly[i];
  real upperbnd = beta_weekly[i+1];

  real approx_t = lowerbnd + (upperbnd-lowerbnd)/(1+exp(-k*(t-(ts[i+1]+ts[i])/2)));

  return approx_t;
}


// infer value at time t between two time points using logistic switch
array[] real between_points_class( array[,] real beta_weekly, real t, data int num_t, array[] real ts ){

  int i = 1;
  real tdif = abs(t - i);
  real k = 100; // slope
  int num_class = dims(beta_weekly)[1];

  while (tdif>1 && i<(size(ts)-1)){
    i = i + 1;
    tdif = abs(t-i);
  }

  vector[num_class] lowerbnd = to_vector( beta_weekly[,i]);
  vector[num_class] upperbnd = to_vector( beta_weekly[,i+1]);

  array[num_class] real approx_t = to_array_1d(lowerbnd + (upperbnd-lowerbnd)/(1+exp(-k*(t-(ts[i+1]+ts[i])/2))));

  return approx_t;
}
