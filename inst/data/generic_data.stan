// data objects that are generic independent of the model specifications
int num_t;
array[num_t] real ts;
int popsize;

// priors
vector[2] p_I0; // expected initial seed (mean, sd)
vector[2] p_R0; // expected beta (alpha, beta)
vector[2] p_theta;

// fixed quantities
real generation_time;

// control parameters
real rtol;
real atol;
int max_num_steps;
int inference;
