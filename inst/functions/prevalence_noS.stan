// Get incidence
vector get_incidence(array[] vector y, data array[] int DIM, data real popsize, data real atol, real I0,
                          vector pi_, data array[] real t_survey_start, data array[] int t_survey_end ) {
  // incidence defined as people moving from I to R
  // function outputs the ascertained incidence rather than the total incidence defined by the model
  // dimensions
  int num_comp = DIM[1];
  int num_t = DIM[2];
  int num_prev = DIM[3];
  int num_seroprev = DIM[size(DIM)];

  // extract, rescale and format prevalence
  vector[num_t] inc;
  inc[1] = pi_[1] *(y[1, num_prev+1] + 1);

  vector[num_t] ascertainment_perweek = rep_each( pi_ , t_survey_end );

  for (i in 2:num_t){ // maybe do a vector multiplication instead to make it faster?
      inc[i] = ascertainment_perweek[i] * (y[i,num_prev+1] - y[i-1,num_prev+1] + 1);
  }

  return inc;
}
