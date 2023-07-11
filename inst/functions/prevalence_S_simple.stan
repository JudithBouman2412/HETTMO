// Incidence function for stratified ODE function
array[,] real get_incidence( array[] vector y , data array[] int DIM,
                         real p_detect1, real p_detect2, data real t_detectionSwitch ) {
  // extract, rescale and format prevalence
  int num_t = DIM[2];
  int num_class = DIM[6];
  array[num_class, num_t] real inc;

  for (j in 1:num_class){
    inc[j,1] = p_detect1 * (y[1][ind(4, j, num_class)] + 1);
    for (i in 2:num_t){ // maybe do a vector multiplication instead to make it faster?
      if (i < t_detectionSwitch){
        inc[j,i] = p_detect1 * (y[i][ind(4, j, num_class)] - y[i-1][ind(4, j, num_class)] + 1);
      } else {
        inc[j,i] = p_detect2 * (y[i][ind(4, j, num_class)] - y[i-1][ind(4, j, num_class)] + 1);
      }
    }
  }

  return inc;
}
