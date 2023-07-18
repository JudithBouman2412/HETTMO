// Incidence function for stratified ODE function
array[,] real get_incidence( array[] vector y , data array[] int DIM,
                         array[,] real pi_, data array[] real t_survey_start, data array[] int t_survey_end) {
  // extract, rescale and format prevalence
  int num_t = DIM[2];
  int num_class = DIM[6];
  int num_serosurvey = DIM[9];
  array[num_class, num_t] real inc;
  array[num_class] vector[num_t] ascertainment_perweek;

  for (j in 1:num_class){
    ascertainment_perweek[j] = rep_each( to_vector(pi_[j,]) , t_survey_end );
  }

  for (j in 1:num_class){
    inc[j,1] = pi_[1] * (y[1][ind(4, j, num_class)] + 1);
    for (i in 2:num_t){ // maybe do a vector multiplication instead to make it faster?
        inc[j,i] = ascertainment_perweek[[j]][i] * (y[i][ind(4, j, num_class)] - y[i-1][ind(4, j, num_class)] + 1);
    }
  }

  return inc;
}
