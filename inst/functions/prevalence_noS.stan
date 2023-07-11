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
