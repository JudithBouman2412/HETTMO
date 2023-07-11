// create eta_weekly
// create Brownian Motion vector for transmission rate
vector[num_t] eta_weekly;
eta_weekly[1] = log(beta);

for (i in 2:num_t) eta_weekly[i] = eta_daily[i-1] + sigmaBM * eta_noise_mat[i-1]; 

// exponential transformation
vector[num_t] beta_weekly = exp(eta_weekly);

