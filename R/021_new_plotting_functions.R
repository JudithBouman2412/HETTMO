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
  dat_ = fit1$samples_posterior$summary(c("I_t_simulated")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  dat_$date = data$date
  post_1 = fit1$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  post_1$date = data$date
  post_2 = fit2$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  post_2$date = data$date
  post_3 = fit3$samples_posterior$summary(c("I_t_predicted")) %>%
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
    serop_date[i,"date"] = mean(c(lubricate::ymd(GE_data[[i]][[2]]),lubricate::ymd(GE_data[[i]][[3]])))
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
    geom_point(data=dat_,aes(y=median,shape="Case counts (Geneva)")) +

    geom_point(data=serop_date,aes(y=prop*scaling,shape="Seroprevalence (Geneva)")) +
    geom_ribbon(data=cum_all,aes(ymin=cum_inf_prop_lwb*scaling,ymax=cum_inf_prop_upb*scaling),fill="grey50",alpha=.3) +
    geom_line(data=cum_all,aes(y=cum_inf_prop*scaling),linetype=2,colour="grey50") +

    facet_wrap(~type,ncol=1) +
    scale_shape_manual(values=c("Case counts (Geneva)"=1,"Seroprevalence (Geneva)"=3)) +
    scale_x_date(date_breaks="2 month", date_labels="%b. %Y") +
    scale_colour_manual(values=cust_cols) +
    scale_fill_manual(values=cust_cols) +
    guides(fill = "none", colour="none")  +
    labs(x="Time",y="Laboratory-confirmed cases",shape=NULL,fill=NULL,colour=NULL) +
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

  cust_cols2 = c("dodgerblue","chartreuse3","orange")

  dat_ = data.frame(variable = "I_t_simulated", time = 1:length(data_cases),
                    null = "", median = data_cases,
                    prob_trans_sim=data_prob)

  # unstratified model fits

  pred_1 = fit1$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Brownian motion")
  pred_2 = fit2$samples_posterior$summary(c("I_t_predicted")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="B-splines")
  pred_3 = fit3$samples_posterior$summary(c("I_t_predicted")) %>%
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
    labs(x="Time",y="Confirmed cases",shape=NULL,fill=NULL,colour=NULL) +
    theme(legend.position=c(.3,.9),
          legend.background = element_blank())

  # probability of transmission
  dat_prob = data.frame(variable = "Transmission_probability",
                        time = 1:length(data_prob), null = "", true = data_prob)
  pred_BM_prob = fit1$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Brownian motion")
  pred_spline_prob = fit2$samples_posterior$summary(c("prob_infection")) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="B-splines")
  pred_GP_prob = fit3$samples_posterior$summary(c("prob_infection")) %>%
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
    labs(shape=NULL,x="Time",y=expression(rho(t)))

  # ascertainment
  dat_asc = data.frame(variable = "Transmission_probability",
                       time = "Week 25 to 45" , null = "", median = c(params$p_detect2))
  asc1 = fit1$samples_posterior$summary(c("p_detect2")) %>%
    dplyr::mutate(type="Brownian motion")
  asc2 = fit2$samples_posterior$summary(c("p_detect2")) %>%
    dplyr::mutate(type="B-splines")
  asc3 = fit3$samples_posterior$summary(c("p_detect2")) %>%
    dplyr::mutate(type="Gaussian process")
  asc_all = dplyr::bind_rows(asc1, asc2, asc3) %>%
    dplyr::mutate(time="Week 25 to 45",
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
    labs(shape=NULL,x=NULL,y=expression(alpha[i])) +
    coord_flip()

  gx = prob_all  %>%
    ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin=q5,ymax=q95,fill=type), alpha=.5) +
    geom_line(aes(y=median,colour=type)) +
    scale_fill_manual(values=cust_cols2) +
    scale_colour_manual(values=cust_cols2) +
    theme(legend.position="right",
          legend.background = element_blank()) +
    labs(shape=NULL,x="Time",y=expression(rho(t)),colour=
           "Time-varying transmission",fill="Time-varying transmission")

  leg = ggpubr::get_legend(gx)

  g = cowplot::plot_grid(
    cowplot::plot_grid(g1,g2,ncol=1,labels=c("B","C")),
    cowplot::plot_grid(g3,leg,ncol=1,rel_heights = c(3,1),labels=c("D","")),
    rel_widths = c(1.5,1),
    nrow=1)

  return(g)
}


