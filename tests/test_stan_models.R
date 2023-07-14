# Testing if updated pacakge works

library(dplyr)

data("GE_data")

weekly_GE <- GE_data[[4]] %>% select(date, cases_total)

first_wave <- weekly_GE$cases_total

# Testing different sampling models
data_list_p = standata( data = first_wave, homogeneous = TRUE, simulated = FALSE,
                         type = "spline", sampling = "p", seroprev_dat = list(GE_data[[1]][[1]], GE_data[[2]][[1]], GE_data[[3]][[1]] ) )

result_p <- sampling_model( data_list = data_list_qp,
                                          model = "poisson_spline_unstrat.stan",
                                          save_result = "data/test_poisson.RDS",
                                          sampler  = 4, rtol=1e-6, atol = 1e-6, max_num_steps=1000,
                                          sampling = "p" )

data_list_qp = standata( data = first_wave, homogeneous = TRUE, simulated = FALSE,
                        type = "spline", sampling = "qp", seroprev_dat = list(GE_data[[1]][[1]], GE_data[[2]][[1]], GE_data[[3]][[1]] ) )

result_qp <- sampling_model( data_list = data_list_qp,
                                         model = "qpoisson_spline_unstrat.stan",
                                         save_result = "data/test_qpoisson.RDS",
                                         sampler  = 4, rtol=1e-6, atol = 1e-6, max_num_steps=1000,
                                         sampling = "qp" )


data_list_nb = standata( data = first_wave, homogeneous = TRUE, simulated = FALSE,
                         type = "spline", sampling = "negbin", seroprev_dat = list(GE_data[[1]][[1]], GE_data[[2]][[1]], GE_data[[3]][[1]] ) )

result_nb <- sampling_model( data_list = data_list_nb,
                                          model = "negbin_spline_unstrat.stan",
                                          save_result = "data/test_negbin.RDS",
                                          sampler  = 4, rtol=1e-6, atol = 1e-6, max_num_steps=1000,
                                          sampling = "negbin" )
# Testing different time-varying methods
data_list_BM = standata( data = first_wave, homogeneous = TRUE, simulated = FALSE,
                         type = "BM", sampling = "qp", seroprev_dat = list(GE_data[[1]][[1]], GE_data[[2]][[1]], GE_data[[3]][[1]] ) )

result_BM <- sampling_model( data_list = data_list_BM,
                             type = "BM",
                                          model = "qpoisson_BM_unstrat.stan",
                                          save_result = "data/test_qpoisson_BM.RDS",
                                          sampler  = 4, rtol=1e-6, atol = 1e-6, max_num_steps=1000,
                                          sampling = "qp" )

data_list_GP = standata( data = first_wave, homogeneous = TRUE, simulated = FALSE,
                         type = "GP", sampling = "qp", seroprev_dat = list(GE_data[[1]][[1]], GE_data[[2]][[1]], GE_data[[3]][[1]] ) )

result_GP <- sampling_model( data_list = data_list_GP,
                             type = "GP",
                             model = "qpoisson_GP_unstrat.stan",
                             save_result = "data/test_qpoisson_GP.RDS",
                             sampler  = 4, rtol=1e-6, atol = 1e-6, max_num_steps=1000,
                             sampling = "qp" )

# Testing stratified models
data(simulated_stratified)
sim_data_3 = simulated_stratified[[1]]
beta_sim_3 = simulated_stratified[[2]]
n_infected_survey_strat = simulated_stratified[[3]]

data_list_strat <-  standata(  sim_data_3,  homogeneous = FALSE,
                                  seroprev_dat = n_infected_survey_strat, additional_data = contact_all,
                                  simulated = TRUE)

result_strat <- sampling_model( data_list = data_list_strat,
                             type = "spline",
                             model = "qpoisson_spline_strat.stan",
                             save_result = "data/test_strat.RDS",
                             sampler  = 1, rtol=1e-6, atol = 1e-6, max_num_steps=1000,
                             sampling = "qp" )

