# fit1 = test_poisson_result_p
# fit2 = test_poisson_result_qp
# fit3 = test_poisson_result_nb
# data = weekly_GE

#' Title
#'
#' @param fit1
#' @param fit2
#' @param fit3
#' @param data
#' @param GE_data
#'
#' @return
#' @export
#'
#' @examples
plot_compare_sampling_GE = function (fit1, fit2, fit3, data, GE_data) {
  ggplot2::theme_set(ggplot2::theme_bw())
  cust_cols = c("chartreuse1","chartreuse3","darkgreen")

  # lab-confirmed cases
  dat_ = fit1$samples_posterior$summary(c("confirmed_cases_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  dat_$date = data$date
  post_1 = fit1$samples_posterior$summary(c("confirmed_cases_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  post_1$date = data$date
  post_2 = fit2$samples_posterior$summary(c("confirmed_cases_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  post_2$date = data$date
  post_3 = fit3$samples_posterior$summary(c("confirmed_cases_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  post_3$date = data$date

  post_1$type = "Poisson"
  post_2$type = "Quasi-Poisson"
  post_3$type = "Negative-binomial"

  post_all = bind_rows(post_1,post_2,post_3)

  # cumulative infected
  popsize = sum(GE_data[[5]])
  serop_date = tibble(pos=NA,tested=NA,date=NA)
  for(i in 1:3) {
    serop_date[i,"pos"] = sum(GE_data[[i]][[1]]$num_pos_tests)
    serop_date[i,"tested"] = sum(GE_data[[i]][[1]]$num_tested)
    serop_date[i,"date"] = mean(c(lubridate::ymd(GE_data[[i]][[2]]),lubridate::ymd(GE_data[[i]][[3]])))
  }
  serop_date$prop = serop_date$pos/serop_date$tested

  cum1 = fit1$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::filter(d2=="4]") %>%
    dplyr::mutate(date=data$date,
                  type="Poisson",
                  cum_inf_prop = median/popsize,
                  cum_inf_prop_lwb = q5/popsize,
                  cum_inf_prop_upb = q95/popsize)
  cum2 = fit2$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::filter(d2=="4]") %>%
    dplyr::mutate(date=data$date,
                  type="Quasi-Poisson",
                  cum_inf_prop = median/popsize,
                  cum_inf_prop_lwb = q5/popsize,
                  cum_inf_prop_upb = q95/popsize)
  cum3 = fit3$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::filter(d2=="4]") %>%
    dplyr::mutate(date=data$date,
                  type="Negative-binomial",
                  cum_inf_prop = median/popsize,
                  cum_inf_prop_lwb = q5/popsize,
                  cum_inf_prop_upb = q95/popsize)
  cum_all = bind_rows(cum1,cum2,cum3)  %>%
    dplyr::mutate(type=factor(type,levels=c("Poisson","Quasi-Poisson","Negative-binomial")))

  # assemble
  scaling = 50000
  g1 = post_all %>%
    dplyr::mutate(type=factor(type,levels=c("Poisson","Quasi-Poisson","Negative-binomial"))) %>%
    ggplot2::ggplot( aes(x=date)) +
    geom_ribbon(aes(ymin=q5,ymax=q95,fill=type), alpha=.5) +
    geom_line(aes(y=median,colour=type)) +
    geom_point(data=dat_,aes(y=median,shape="Confirmed cases (Geneva)")) +

    geom_point(data=serop_date,aes(y=prop*scaling,shape="Seroprevalence (Geneva)")) +
    geom_ribbon(data=cum_all,aes(ymin=cum_inf_prop_lwb*scaling,ymax=cum_inf_prop_upb*scaling),fill="grey50",alpha=.3) +
    geom_line(data=cum_all,aes(y=cum_inf_prop*scaling),linetype=2,colour="grey50") +

    facet_wrap(~type,ncol=1) +
    scale_shape_manual(values=c("Confirmed cases (Geneva)"=1,"Seroprevalence (Geneva)"=3)) +
    scale_x_date(date_breaks="2 month", date_labels="%b. %Y") +
    scale_colour_manual(values=cust_cols) +
    scale_fill_manual(values=cust_cols) +
    guides(fill = "none", colour="none")  +
    labs(x=element_blank(),y="Laboratory-confirmed cases",shape=NULL,fill=NULL,colour=NULL) +
    theme(legend.position=c(.3,.95),
          legend.background = element_blank()) +
    scale_y_continuous(
      "Laboratory-confirmed cases",
      sec.axis = sec_axis(~ . /scaling, name = "Cumulative infected",labels=scales::percent))

  return(g1)


  #
  # if(FALSE) {
  #   g1 = ggplot2::ggplot() + ggplot2::geom_ribbon(data = post_1,
  #                                                 aes(x = date, ymin = q5, ymax = q95, fill = "Posterior predictive check"),
  #                                                 alpha = 0.3) + ggplot2::geom_line(data = post_1, aes(x = date,
  #                                                                                                      y = median), colour = col_post) + ggplot2::geom_point(data = dat_,
  #                                                                                                                                                            aes(x = date, y = median, colour = "Case data Bern")) +
  #     ggplot2::scale_fill_manual(values = c(`Posterior predictive check` = col_post)) +
  #     ggplot2::scale_colour_manual(values = c(`Case data Bern` = col_simul)) +
  #     ggplot2::labs(x = "Time", y = "N",
  #                   colour = NULL, fill = NULL) +
  #     ggplot2::theme(legend.position = "bottom", text = ggplot2::element_text(size = 16))
  #   g2 = ggplot2::ggplot() + ggplot2::geom_ribbon(data = post_2,
  #                                                 aes(x = date, ymin = q5, ymax = q95, fill = "Posterior predictive check"),
  #                                                 alpha = 0.3) + ggplot2::geom_line(data = post_2, aes(x = date,
  #                                                                                                      y = median), colour = col_post) + ggplot2::geom_point(data = dat_,
  #                                                                                                                                                            aes(x = date, y = median, colour = "Case data Bern")) +
  #     ggplot2::scale_fill_manual(values = c(`Posterior predictive check` = col_post)) +
  #     ggplot2::scale_colour_manual(values = c(`Case data Bern` = col_simul)) +
  #     ggplot2::labs(x = "Time (weeks)", y = "# pos. tests",
  #                   colour = NULL, fill = NULL, title = " Quasi Poisson model") +
  #     ggplot2::theme(legend.position = "bottom", text = ggplot2::element_text(size = 16))
  #   g3 = ggplot2::ggplot() + ggplot2::geom_ribbon(data = post_3,
  #                                                 aes(x = date, ymin = q5, ymax = q95, fill = "Posterior predictive check"),
  #                                                 alpha = 0.3) + ggplot2::geom_line(data = post_3, aes(x = date,
  #                                                                                                      y = median), colour = col_post) + ggplot2::geom_point(data = dat_,
  #                                                                                                                                                            aes(x = date, y = median, colour = "Case data Bern")) +
  #     ggplot2::scale_fill_manual(values = c(`Posterior predictive check` = col_post)) +
  #     ggplot2::scale_colour_manual(values = c(`Case data Bern` = col_simul)) +
  #     ggplot2::labs(x = "Time (weeks)", y = "# pos. tests",
  #                   colour = NULL, fill = NULL, title = "Negative binomial model") +
  #     ggplot2::theme(legend.position = "bottom", text = ggplot2::element_text(size = 16))
  #   prob_p = fit1$samples_posterior$summary(c("prob_infection")) %>%
  #     tidyr::separate(variable, "\\[|\\]", into = c("variable",
  #                                                   "time", "null"))
  #   prob_qp = fit2$samples_posterior$summary(c("prob_infection")) %>%
  #     tidyr::separate(variable, "\\[|\\]", into = c("variable",
  #                                                   "time", "null"))
  #   prob_nb = fit3$samples_posterior$summary(c("prob_infection")) %>%
  #     tidyr::separate(variable, "\\[|\\]", into = c("variable",
  #                                                   "time", "null"))
  #   prob_plot <- ggplot2::ggplot() + ggplot2::geom_ribbon(data = prob_p,
  #                                                         aes(x = as.numeric(time), ymin = q5, ymax = q95, fill = "Poisson model"),
  #                                                         alpha = 0.3) + ggplot2::geom_ribbon(data = prob_qp,
  #                                                                                             aes(x = as.numeric(time), ymin = q5, ymax = q95, fill = "Quasi Poisson model"),
  #                                                                                             alpha = 0.3) + ggplot2::geom_ribbon(data = prob_nb,
  #                                                                                                                                 aes(x = as.numeric(time), ymin = q5, ymax = q95, fill = "Negative Binomial model"),
  #                                                                                                                                 alpha = 0.3) + ggplot2::geom_line(data = prob_p, aes(x = as.numeric(time),
  #                                                                                                                                                                                      y = median, colour = "Poisson model")) + ggplot2::geom_line(data = prob_qp,
  #                                                                                                                                                                                                                                                  aes(x = as.numeric(time), y = median, colour = "Quasi Poisson model")) +
  #     ggplot2::geom_line(data = prob_nb, aes(x = as.numeric(time),
  #                                            y = median, colour = "Negative Binomial model")) +
  #     ggplot2::scale_fill_manual(values = c(`Poisson model` = "Orange",
  #                                           `Quasi Poisson model` = "Red", `Negative Binomial model` = "blue")) +
  #     ggplot2::scale_colour_manual(values = c(`Poisson model` = "Orange",
  #                                             `Quasi Poisson model` = "Red", `Negative Binomial model` = "blue")) +
  #     ggplot2::labs(x = "Time (weeks)", y = "Probability of transmission upon contact",
  #                   colour = NULL, fill = NULL) + ggplot2::theme(legend.position = "bottom",
  #                                                                text = ggplot2::element_text(size = 16))
  #   legend <- cowplot::get_legend(g1 + ggplot2::theme(legend.box.margin = ggplot2::margin(0,
  #                                                                                         0, 0, 12), legend.position = "bottom"))
  #   prow <- cowplot::plot_grid(g1 + ggplot2::theme(legend.position = "none"),
  #                              g2 + ggplot2::theme(legend.position = "none"), g3 +
  #                                ggplot2::theme(legend.position = "none"), nrow = 3,
  #                              labels = LETTERS)
  #   full_plot = cowplot::plot_grid(prow, legend, nrow = 2, rel_widths = c(3,
  #                                                                         0.4), rel_heights = c(3, 0.2))
  #   return(list(full_plot, prob_plot))
  # }
}

# summ = comparison_summary

#' Title
#'
#' @param summ
#'
#' @return
#' @export
#'
#' @examples
plot_unstrat_benchmark = function(summ) {
  cust_cols2 = c("dodgerblue","chartreuse3","orange")
  cust_shapes = c("1"=8,
                  "2"=7,
                  "3"=0,
                  "4"=6,
                  "5"=5,
                  "NA"=1)
  summ = summ %>%
    dplyr::mutate(type2=factor(type,levels=c("BM","spline","GP"),labels=c("Brownian motion","B-splines","Gaussian process")),
                  knot_sequence=ifelse(knot_sequence=="No knots","NA",knot_sequence))
  g1 = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = RMSE_prob, colour = type2, shape = knot_sequence),size=2) +
    scale_colour_manual(values = cust_cols2,guide="none")+
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes,guide="none") +
    labs(x="ESS per second", y="Error" , colour=NULL, fill=NULL, shape="Knot sequence") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position=c(.25,.8),
          legend.background = element_blank())

  g2 = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = CI_size_prob, colour = type2, shape = knot_sequence),size=2) +
    scale_colour_manual(values = cust_cols2,guide="none")+
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes) +
    labs(x="ESS per second", y="Sharpness" , colour=NULL, fill=NULL, shape="Knot sequence") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position=c(.25,.4),
          legend.background = element_blank())
  g = cowplot::plot_grid(g1,g2,nrow=1,labels=c("E","F"))

  return(g)
}

#' Title
#'
#' @param fit1
#' @param fit2
#' @param fit3
#' @param data_cases
#' @param data_prob
#' @param params
#'
#' @return
#' @export
#'
#' @examples
plot_unstrat_model_fit = function(fit1,fit2,fit3,data_cases,data_prob,params) {
  ggplot2::theme_set(ggplot2::theme_bw())
  cust_cols2 = c("dodgerblue","chartreuse3","orange")

  dat_ = data.frame(variable = "confirmed_cases_simulated", time = 1:length(data_cases),
                    null = "", median = data_cases,
                    prob_trans_sim=data_prob)

  # unstratified model fits

  pred_1 = fit1$samples_posterior$summary(c("confirmed_cases_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Brownian motion")
  pred_2 = fit2$samples_posterior$summary(c("confirmed_cases_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="B-splines")
  pred_3 = fit3$samples_posterior$summary(c("confirmed_cases_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Gaussian process")
  pred_all = dplyr::bind_rows(pred_1, pred_2, pred_3) %>%
    dplyr::mutate(time=as.numeric(time),
                  type=factor(type,levels=c("Brownian motion","B-splines","Gaussian process")))

  g1 = pred_all %>%
    ggplot2::ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin = q5, ymax = q95, fill = type), alpha = 0.3) +
    # geom_line(aes(y = q5, colour = type),linetype=2) +
    # geom_line(aes(y = q95, colour = type),linetype=2) +

    geom_line(aes(y = median, colour = type)) +
    ggplot2::geom_point(data = dat_, aes(y = median,shape="Simulated data"))  +
    scale_shape_manual(values=c("Simulated data"=4)) +
    scale_fill_manual(values=cust_cols2) +
    scale_colour_manual(values=cust_cols2) +
    guides(fill = "none", colour="none")  +
    labs(x="Time (weeks)",y="Confirmed cases",shape=NULL,fill=NULL,colour=NULL) +
    theme(legend.position=c(.3,.9),
          legend.background = element_blank())

  # probability of transmission
  dat_prob = data.frame(variable = "rho",
                        time = 1:length(data_prob), null = "", true = data_prob)
  pred_BM_prob = fit1$samples_posterior$summary(c("rho")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Brownian motion")
  pred_spline_prob = fit2$samples_posterior$summary(c("rho")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="B-splines")
  pred_GP_prob = fit3$samples_posterior$summary(c("rho")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Gaussian process")
  prob_all = dplyr::bind_rows(pred_BM_prob, pred_spline_prob, pred_GP_prob) %>%
    dplyr::mutate(time=as.numeric(time),
                  type=factor(type,levels=c("Brownian motion","B-splines","Gaussian process")))

  g2 = prob_all  %>%
    ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin=q5,ymax=q95,fill=type), alpha=.5) +
    geom_line(aes(y=median,colour=type)) +

    geom_point(data=dat_prob,aes(y=true,shape="True value")) +

    scale_shape_manual(values=c("True value"=4)) +
    scale_fill_manual(values=cust_cols2,guide="none") +
    scale_colour_manual(values=cust_cols2,guide="none") +
    theme(legend.position=c(.3,.9),
          legend.background = element_blank()) +
    labs(shape=NULL,x="Time (weeks)",y=expression(rho(t)))

  # ascertainment
  dat_asc = data.frame(variable = "Ascertainment rate",
                       time = c("Week 0 to 19", "Week 20 to 45") ,
                       null = "",
                       median = c(params$p_detect1, params$p_detect2))
  asc1 = fit1$samples_posterior$summary(c("pi_")) %>%
    dplyr::mutate(type="Brownian motion",
                  variable = "Ascertainment rate")
  asc2 = fit2$samples_posterior$summary(c("pi_")) %>%
    dplyr::mutate(type="B-splines",
                  variable = "Ascertainment rate")
  asc3 = fit3$samples_posterior$summary(c("pi_")) %>%
    dplyr::mutate(type="Gaussian process",
                  variable = "Ascertainment rate")
  asc_all = dplyr::bind_rows(asc1, asc2, asc3) %>%
    dplyr::mutate(time=rep(c("Week 0 to 19", "Week 20 to 45"), 3),
                  type=factor(type,levels=c("Brownian motion","B-splines","Gaussian process")))

  g3 = asc_all %>%
    left_join(dplyr::select(dat_asc,time,true=median),by = join_by(time)) %>%
    ggplot(aes(x=as.factor(time))) +
    geom_pointrange(aes(y=median,ymin=q5,ymax=q95,colour=type),
                    position=position_dodge(.8)) +
    geom_point(aes(y=true,shape="True value",group=type),
               position=position_dodge(.8)) +
    scale_shape_manual(values=c("True value"=4)) +
    scale_colour_manual(values=cust_cols2,guide="none") +
    scale_y_continuous(labels=scales::percent) +
    theme(legend.position=c(.3,.9),
          legend.background = element_blank(),
          axis.text.y = element_text(angle=90,hjust=.5)) +
    labs(shape=NULL,x=NULL,y=expression(pi[i])) +
    coord_flip()

  gx = prob_all  %>%
    ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin=q5,ymax=q95,fill=type), alpha=.5) +
    geom_line(aes(y=median,colour=type)) +
    scale_fill_manual(values=cust_cols2) +
    scale_colour_manual(values=cust_cols2) +
    theme(legend.position="right",
          legend.background = element_blank()) +
    labs(shape=NULL,x="Time (weeks)",y=expression(rho(t)),colour=
           "Time-varying transmission",fill="Time-varying transmission")

  leg = ggpubr::get_legend(gx)

  g = cowplot::plot_grid(
    cowplot::plot_grid(g1,g2,ncol=1,labels=c("B","C")),
    cowplot::plot_grid(g3,leg,ncol=1,rel_heights = c(3,1),labels=c("D","")),
    rel_widths = c(1.5,1),
    nrow=1)

  return(g)
}

#' Index finder
#'
#' @param comp the compartment of the S (1) I (2) E (3) R (4) of which the index is required
#' @param age_sex the required stratification class
#' @param num_class the total number of classes in the model (4)
#'
#' @return
#' @export
#'
#' @examples
ind <- function( comp, age_sex, num_class=4){
  return(age_sex+(comp-1)*num_class)
}

#' plot_fitsim_strat
#'
#' Function that produces Figure to visualize stratified model fit to simulated data.
#'
#' @param fit fit of stratified model, result of "sampling_model"
#' @param sim_data array of simulated positive PCR cases per week and age-class.
#' @param params list of parameters used to simulate the data, result of "set_parameters".
#' @param transmission_prob_sim array of simulated rho value per week per age-class.
#'
#' @return Figure of result of simulated, stratified model.
#' @export
#'
#' @examples
plot_fitsim_strat = function(fit, sim_data, params=params, transmission_prob_sim) {

  # aesthetics
  ggplot2::theme_set(ggplot2::theme_bw())
  cust_cols = c("cyan3","pink","purple")

  # data from stratified simulation
  popsize = params$popsize
  n_tested = params$n_tested_survey

  # posterior predictive check
  dat_ = rbind( data.frame(variable = "I_t_simulated",
                           time = 1:dim(sim_data)[1],
                           age_group = "0-19 years",
                           median = sim_data[,1] ),
                data.frame(variable = "I_t_simulated",
                           time = 1:dim(sim_data)[1],
                           age_group = "20-59 years",
                           median = sim_data[,2] ),
                data.frame(variable = "I_t_simulated",
                           time = 1:dim(sim_data)[1],
                           age_group = "60+ years",
                           median = sim_data[,3] ))

  post_ = fit$samples_posterior$summary(c("I_t_predicted")) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>% mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-59 years", "60+ years")))
  prior_ = fit$samples_prior$summary(c("I_t_predicted")) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time"))

  # cumulative infected
  comp_4 <- c(ind(4,1,3), ind(4,2,3), ind(4,3,3))
  cum_ = fit$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::mutate(comp = as.numeric( str_sub(d2, end = -2) ) ) %>%
    dplyr::filter(comp %in% comp_4 ) %>%
    dplyr::mutate( age_group = ifelse(comp==comp_4[1], "0-19 years", ifelse( comp==comp_4[2], "20-59 years", "60+ years") ),
                   time = rep(1:dim(sim_data)[1] , 3),
                   cum_inf_prop = ifelse(comp == comp_4[1], median/popsize[1], ifelse(comp == comp_4[2], median/popsize[2], median/popsize[3]) ),
                   cum_inf_prop_lwb = ifelse(comp == comp_4[1], q5/popsize[1], ifelse(comp == comp_4[2], q5/popsize[2], q5/popsize[3]) ),
                   cum_inf_prop_upb = ifelse(comp == comp_4[1], q95/popsize[1], ifelse(comp == comp_4[2], q95/popsize[2], q95/popsize[3]) ) )

  data(simulated_stratified)
  n_infected_survey_strat = simulated_stratified[[3]]
  sim_data_3 = simulated_stratified[[1]]

  serop_date = data.frame(time = c(rep(25, 3), rep(45,3)), age_group = c("0-19 years","20-59 years", "60+ years"))
  serop_date$prop = c(rep(NA, 3,) , n_infected_survey_strat/params$n_tested_survey)

  scaling = 5000
  g1_A = ggplot() +
    geom_ribbon(data=post_ , aes(x=as.numeric(time),ymin=q5,ymax=q95,fill=age_group),alpha=.3) +
    geom_line(data=post_ , aes(x=as.numeric(time),y=median,colour=age_group)) +
    geom_point(data=dat_, aes(x=as.numeric(time),y=median,colour=age_group,shape="Confirmed cases (simulated)") ) +
    geom_ribbon(data=cum_, aes(ymin=cum_inf_prop_lwb*scaling, ymax=cum_inf_prop_upb*scaling, x = time),fill="grey50",alpha=.3) +
    geom_line(data=cum_, aes(y=cum_inf_prop*scaling, x=time), linetype=2, colour="grey50") +
    facet_wrap(~age_group,ncol=1) +
    geom_point(data = serop_date, aes(x =time, y = prop*scaling,shape="Seroprevalence (simulated)")) +
    # add cumulative cases and the seroprevalence value
    scale_colour_manual(values=cust_cols) +
    scale_fill_manual(values=cust_cols) +
    scale_shape_manual(values=c("Confirmed cases (simulated)"=4,"Seroprevalence (simulated)"=3)) +
    labs(x="Time (weeks)",y="Laboratory-confirmed cases",shape=NULL,fill=NULL,colour=NULL) +
    theme(text=element_text(size=16), legend.position=c(.45,.25),legend.background = element_blank()) +
    scale_y_continuous(
      "Laboratory-confirmed cases",
      sec.axis = sec_axis(~ . /scaling, name = "Cumulative infected",labels=scales::percent)) +
    guides(colour=FALSE, fill=FALSE)

  dat_prob = rbind( data.frame(variable = "I_t_simulated",
                               time = 1:dim(sim_data)[1],
                               age_group = "0-19 years",
                               median = transmission_prob_sim[,1] ),
                    data.frame(variable = "I_t_simulated",
                               time = 1:dim(sim_data)[1],
                               age_group = "20-59 years",
                               median = transmission_prob_sim[,2] ),
                    data.frame(variable = "I_t_simulated",
                               time = 1:dim(sim_data)[1],
                               age_group = "60+ years",
                               median = transmission_prob_sim[,3] ))

  post_prob = fit$samples_posterior$summary(c("prob_infection")) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-59 years", "60+ years")))
  prior_prob = fit$samples_prior$summary(c("prob_infection")) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time"))

  # plot of probability
  g1_B = ggplot() +
    geom_ribbon(data=post_prob , aes(x=as.numeric(time),ymin=q5,ymax=q95,fill=age_group),alpha=.3) +
    geom_line(data=post_prob , aes(x=as.numeric(time),y=median,colour=age_group)) +
    geom_point(data=dat_prob, aes(x=as.numeric(time),y=median,colour=age_group, shape= "True value") ) +
    facet_wrap(~age_group,ncol=1) +
    # add cumulative cases and the seroprevalence value
    scale_colour_manual(values=cust_cols,guide="none") +
    scale_shape_manual(values=c("True value"=4)) +
    scale_fill_manual(values=cust_cols,guide="none") +
    labs(x="Time (weeks)",y=expression(rho(t)),shape=NULL,fill=NULL,colour=NULL) +
    theme(text=element_text(size=16), legend.position=c(.25,.25), legend.background = element_blank())

  # plot of ascertainment rate
  # ascertainment
  dat_asc = data.frame(variable = "Transmission_probability",
                       time = "Week 25 to 45" , null = "", median = c(params$p_detect2))
  asc = fit$samples_posterior$summary(c("p_detect2")) %>%
    dplyr::mutate(time = "Week 25 to 45" ,
                  age_group = "0.19 years")

  g3 = asc %>%
    left_join(dplyr::select(dat_asc,time,true=median),by = join_by(time)) %>%
    ggplot(aes(x=as.factor(time))) +
    geom_pointrange(aes(y=median,ymin=q5,ymax=q95,colour=age_group),
                    position=position_dodge(.8)) +
    geom_point(aes(y=true,shape="True value",group=age_group),
               position=position_dodge(.8)) +
    scale_shape_manual(values=c("True value"=4)) +
    scale_colour_manual(values=cust_cols,guide="none") +
    scale_y_continuous(labels=scales::percent) +
    theme(legend.position=c(0.15, 0.85),
          legend.background = element_blank(),
          axis.text.y = element_text(angle=90,hjust=.5),text=element_text(size=16)) +
    labs(shape=NULL,x=NULL,y=expression(pi[i]))

  full_plot <- cowplot::plot_grid(g1_A,g1_B,nrow=1,labels=c("A","B"), rel_widths = c(1.25,1))
  full_plot2 <- cowplot::plot_grid(full_plot,g3,nrow=2,labels=c("","C"), rel_heights = c(3,1))

  return(full_plot2)
}


#' Title
#'
#' @param fit
#' @param data
#' @param cum_inc
#' @param params
#'
#' @return
#' @export
#'
#' @examples
plot_fitsim_strat_GE <- function(fit, data, cum_inc, params){

  # aesthetics
  ggplot2::theme_set(ggplot2::theme_bw())
  cust_cols = c("cyan3","pink","purple")

  # data from stratified simulation
  popsize = params$popdist
  n_tested = params$n_tested_survey

  # posterior predictive check
  dat_ = rbind( data.frame(variable = "GE_data",
                           time = data$date,
                           age_group = "0-19 years",
                           median = data[,2] )%>% rename(median = cases_age1 ),
                data.frame(variable = "GE_data",
                           time = data$date,
                           age_group = "20-59 years",
                           median = data[,3] ) %>% rename(median = cases_age2 ),
                data.frame(variable = "GE_data",
                           time = data$date,
                           age_group = "60+ years",
                           median = data[,4] )%>% rename(median = cases_age3 ))

  post_ = fit$samples_posterior$summary(c("I_t_predicted")) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-59 years", "60+ years")))
  post_$time = dat_$time
  prior_ = fit$samples_prior$summary(c("I_t_predicted")) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time"))
  prior_$time = dat_$time

  # cumulative infected
  comp_4 <- c(ind(4,1,3), ind(4,2,3), ind(4,3,3))
  cum_ = fit$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::mutate(comp = as.numeric( str_sub(d2, end = -2) ) ) %>%
    dplyr::filter(comp %in% comp_4 ) %>%
    dplyr::mutate( age_group = ifelse(comp==comp_4[1], "0-19 years", ifelse( comp==comp_4[2], "20-59 years", "60+ years") ),
                   time = rep(data$date, 3),
                   cum_inf_prop = ifelse(comp == comp_4[1], median/popsize[1], ifelse(comp == comp_4[2], median/popsize[2], median/popsize[3]) ),
                   cum_inf_prop_lwb = ifelse(comp == comp_4[1], q5/popsize[1], ifelse(comp == comp_4[2], q5/popsize[2], q5/popsize[3]) ),
                   cum_inf_prop_upb = ifelse(comp == comp_4[1], q95/popsize[1], ifelse(comp == comp_4[2], q95/popsize[2], q95/popsize[3]) ) )

  n_infected_survey_strat = params$n_infected_survey

  # cumulative infected
  serop_date = tibble(pos=NA,tested=NA,date=NA, age_group = NA, prop = NA)
  age_groups = c("0-19 years", "20-59 years", "60+ years")
  check_group = c("group1", "group2", "group3")
  for(j in 1:3) { # for loop over age:groups
    for (i in 1:3){ # for loop over seroprevalence studies
      if ( dim(GE_data[[i]][[1]] %>% filter(age_group==check_group[j]))[1]>0 ){
        serop_date[(j-1)*3 +i,"pos"] = GE_data[[i]][[1]] %>% filter(age_group==check_group[j]) %>% select(num_pos_tests)
        serop_date[(j-1)*3 +i,"tested"] = GE_data[[i]][[1]]%>% filter(age_group==check_group[j]) %>% select(num_tested)
        serop_date[(j-1)*3 +i,"date"] = mean(c(lubridate::ymd(GE_data[[i]][[2]]),lubridate::ymd(GE_data[[i]][[3]])))
        serop_date[(j-1)*3 +i, "prop"] =  serop_date[(j-1)*3 +i,"pos"]/serop_date[(j-1)*3 +i,"tested"]
      } else {
        serop_date[(j-1)*3 +i,"pos"] = NA
        serop_date[(j-1)*3 +i,"tested"] = NA
        serop_date[(j-1)*3 +i,"date"] = NA
        serop_date[(j-1)*3 +i, "prop"] = NA
      }
      serop_date[(j-1)*3 +i, "age_group"] =  age_groups[j]
    }
  }

  scaling = 15000
  g1_A = ggplot() +
    geom_ribbon(data=post_ , aes(x=time,ymin=q5,ymax=q95,fill=age_group),alpha=.3) +
    geom_line(data=post_ , aes(x=time,y=median,colour=age_group)) +
    geom_point(data=dat_, aes(x=time,y=median,colour=age_group,shape="Confirmed cases (Geneva)") ) +
    geom_ribbon(data=cum_, aes(ymin=cum_inf_prop_lwb*scaling, ymax=cum_inf_prop_upb*scaling, x = time),fill="grey50",alpha=.3) +
    geom_line(data=cum_, aes(y=cum_inf_prop*scaling, x=time), linetype=2, colour="grey50") +
    facet_wrap(~age_group,ncol=1) +
    geom_point(data = serop_date, aes(x =date, y = prop*scaling,shape="Seroprevalence (Geneva)")) +
    # add cumulative cases and the seroprevalence value
    scale_colour_manual(values=cust_cols) +
    scale_fill_manual(values=cust_cols) +
    scale_shape_manual(values=c("Confirmed cases (Geneva)"=1,"Seroprevalence (Geneva)"=3)) +
    labs(x="Time",y="Laboratory-confirmed cases",shape=NULL,fill=NULL,colour=NULL) +
    theme(text=element_text(size=16), legend.position=c(.45,.95),legend.background = element_blank()) +
    scale_y_continuous(
      "Laboratory-confirmed cases",
      sec.axis = sec_axis(~ . /scaling, name = "Cumulative infected",labels=scales::percent)) +
    guides(colour=FALSE, fill=FALSE)

  post_prob = fit$samples_posterior$summary(c("prob_infection")) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-59 years", "60+ years")))
  post_prob$date = rep(data$date,each = 3)
  prior_prob = fit$samples_prior$summary(c("prob_infection")) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time"))
  prior_prob$date = rep(data$date,each = 3)

  # plot of probability
  g1_B = ggplot() +
    geom_ribbon(data=post_prob , aes(x=date,ymin=q5,ymax=q95,fill=age_group),alpha=.3) +
    geom_line(data=post_prob , aes(x=date,y=median,colour=age_group)) +
    facet_wrap(~age_group,ncol=1) +
    # add cumulative cases and the seroprevalence value
    scale_colour_manual(values=cust_cols,guide="none") +
    scale_fill_manual(values=cust_cols,guide="none") +
    labs(x=element_blank(),y=expression(rho(t)),shape=NULL,fill=NULL,colour=NULL) +
    theme(text=element_text(size=16), legend.position=c(.25,.25), legend.background = element_blank())

  # plot of ascertainment rate
  # ascertainment
  asc = fit$samples_posterior$summary(c("p_detect1","p_detect2","p_detect3")) %>%
    tidyr::separate(variable,sep=c(8,9,10,11),into=c("variable","seroprevalence", "d1", "age_group")) %>%
    dplyr::mutate( age_group = ifelse(age_group == 1,  "0-19 years", ifelse(age_group==2, "20-59 years", "60+ years") ) )
  asc$date = rep(serop_date$date[7:9], each = 3)

  g3 = asc %>% ggplot(aes(x=as.factor(date))) +
    geom_pointrange(aes(y=median,ymin=q5,ymax=q95,colour=age_group),
                    position=position_dodge(.8)) +
    scale_colour_manual(values=cust_cols,guide="none") +
    scale_y_continuous(labels=scales::percent) +
    theme(legend.position=c(0.15, 0.85),
          legend.background = element_blank(),
          axis.text.y = element_text(angle=90,hjust=.5),text=element_text(size=16)) +
    labs(shape=NULL,x=NULL,y=expression(pi[i]))

  full_plot <- cowplot::plot_grid(g1_A,g1_B,nrow=1,labels=c("A","B"), rel_widths = c(1.25,1))
  full_plot2 <- cowplot::plot_grid(full_plot,g3,nrow=2,labels=c("","C"), rel_heights = c(3,1))

  return(full_plot2)
}


#' Title
#'
#' @param summ
#'
#' @return
#' @export
#'
#' @examples
plot_single_benchmark = function(summ) {

  ggplot2::theme_set(ggplot2::theme_bw())
  cust_cols2 = c("cyan","violet","darkgrey", "lightgreen", "purple")
  cust_shapes = c("300"=2,
                  "500"=0)

  summ <- summ %>% rename(Tolerance = tol)

  g1 = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = RMSE_prob, colour = solver, shape = warmup_iter),size=2) +
    scale_colour_manual(values = cust_cols2, guide = "none", labels=c("0"="rk45","1"="Adams","2"="Bdf", "3"="ckrk","4"="Trapeziodal"))+
    facet_grid(rows = vars(Tolerance), labeller = label_both ) +
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes ) +
    labs(x="ESS per second", y="Error" , colour="ODE solver", fill=NULL, shape="Warmup iterations") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position="bottom",  #c(.38,.9),
          legend.background = element_blank(),text=element_text(size=16), strip.text.y = element_blank() )

  g2 = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = CI_size_prob, colour = solver , shape = warmup_iter),size=2) +
    scale_colour_manual(values = cust_cols2, labels=c("0"="rk45","1"="Adams","2"="Bdf", "3"="ckrk","4"="Trapeziodal"))+
    facet_grid(rows = vars(Tolerance), labeller = label_both ) +
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes , guide = "none") +
    labs(x="ESS per second", y="Sharpness" , colour="ODE solver", fill=NULL, shape="Warmup iterations") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position="bottom", #c(.28,.86),
          legend.background = element_blank(), text=element_text(size=16))

  g = cowplot::plot_grid(g1,g2,nrow=1,labels=c("A","B"))

  return(g)
}

plot_unstrat_GP = function(summ) {
  cust_cols2 = c("dodgerblue","chartreuse3","orange")
  cust_shapes = c("1"=8,
                  "2"=7,
                  "3"=0,
                  "4"=6,
                  "5"=5,
                  "NA"=1)
  summ = summ %>%
    dplyr::mutate(type2=factor(type,levels=c("BM","spline","GP"),labels=c("Brownian motion","B-splines","Gaussian process")),
                  knot_sequence=ifelse(knot_sequence=="No knots","NA",knot_sequence))
  g1 = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = RMSE_prob, colour = type2, shape = knot_sequence),size=2) +
    scale_colour_manual(values = cust_cols2,guide="none")+
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes,guide="none") +
    labs(x="ESS per second", y="Error" , colour=NULL, fill=NULL, shape="Knot sequence") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position=c(.25,.8),
          legend.background = element_blank(), text=element_text(size=16))

  g2 = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = CI_size_prob, colour = type2, shape = knot_sequence),size=2) +
    scale_colour_manual(values = cust_cols2,guide="none")+
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes) +
    labs(x="ESS per second", y="Sharpness" , colour=NULL, fill=NULL, shape="Knot sequence") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position=c(.25,.4),
          legend.background = element_blank(), text=element_text(size=16))
  g = cowplot::plot_grid(g1,g2,nrow=1,labels=c("E","F"))

  return(g)
}

