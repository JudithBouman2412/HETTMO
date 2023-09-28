#' plot_fitsim
#'
#' Function that produces Figure to visualize model fit.
#'
#' @param fit modelfit, result of "sampling_model"
#' @param data data to which model was fit, vector with value per week
#' @param params list of parameters used to fit model
#'
#' @return
#' @export
#'
#' @examples
plot_fitsim = function(fit, data, params) {

  # aesthetics
  ggplot2::theme_set(ggplot2::theme_bw())
  bayesplot::color_scheme_set("teal")
  col_prior = "grey70"
    col_post = "orange"
      col_simul = "cyan3"
        col_cum = "pink"

        # posterior predictive check
        dat_ = data.frame(variable = "confirmed_cases_predicted", time = 1:length(data), null = "", median = data )
        pred_ = fit$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975)) ) |>
          tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
        post_ = fit$samples_prior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
          tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
        cum_post = fit$samples_posterior$summary("y", median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
          tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
          dplyr::filter(d2=="4]") %>%
          dplyr::mutate(time = 1:45,
                        cum_inf_prop = median/params$popsize,
                        cum_inf_prop_lwb = `2.5%`/params$popsize,
                        cum_inf_prop_upb = `97.5%`/params$popsize)

        data(simulated_nonstratified)
        data_sim_sero <- data.frame(time = c(20, 45), tot_tests = c(5000, 5000), pos_test = simulated_nonstratified[[3]] ) %>%
          mutate(recovered = pos_test/tot_tests*params$popsize)

        g1 = ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=post_,aes(x=as.numeric(time),ymin=`2.5%`,ymax=`97.5%`,fill="Prior predictive check"),alpha=.3) +
          ggplot2::geom_ribbon(data=pred_,aes(x=as.numeric(time),ymin=`2.5%`,ymax=`97.5%`,fill="Posterior predictive check"),alpha=.3) +
          ggplot2::geom_line(data=pred_,aes(x=as.numeric(time),y=median),colour=col_post) +
          ggplot2::geom_point(data=dat_,aes(x=as.numeric(time),y=median,colour="Simulated data")) +
          ggplot2::geom_ribbon(data = cum_post, aes(x = as.numeric(time), ymin = `2.5%`, ymax = `97.5%`, fill = "Recovered individuals")) +
          ggplot2::geom_line(data=cum_post, aes(x = as.numeric(time), y = median)) +
          ggplot2::geom_point(data = data_sim_sero, aes(x = time, y = recovered), shape = 6, col = "red" ) +
          ggplot2::scale_fill_manual(values=c("Prior predictive check"=col_prior,
                                              "Posterior predictive check"=col_post,
                                              "Recovered individuals" = col_cum)) +
          ggplot2::scale_colour_manual(values=c("Simulated data"=col_simul)) +
          ggplot2::labs(x="Time (weeks)",y="Incidence cases",colour=NULL,fill=NULL) +
          ggplot2::theme(legend.position = "bottom", text=ggplot2::element_text(size=16)) +
          ggplot2::ylim(c(0,27000))

        pred_prob = fit$samples_posterior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
          tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
        post_prob = fit$samples_prior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
          tidyr::separate(variable,"\\[|\\]",into=c("variable","time","null"))
        g2 = ggplot2::ggplot() +
          ggplot2::geom_ribbon(data=post_prob,aes(x=as.numeric(time),ymin=`2.5%`,ymax=`97.5%`,fill="Prior predictive check"),alpha=.3) +
          ggplot2::geom_ribbon(data=pred_prob,aes(x=as.numeric(time),ymin=`2.5%`,ymax=`97.5%`,fill="Posterior predictive check"),alpha=.3) +
          ggplot2::geom_line(data=pred_prob,aes(x=as.numeric(time),y=median),colour=col_post) +
          ggplot2::scale_fill_manual(values=c("Prior predictive check"=col_prior,
                                              "Posterior predictive check"=col_post)) +
          #scale_y_continuous(trans='log10') +
          ggplot2::labs(x="Time (weeks)",y="Incidence cases",colour=NULL,fill=NULL) +
          ggplot2::theme(legend.position = "bottom", text=ggplot2::element_text(size=16))

        post_param = fit$samples_posterior$summary(c("R0","I0_raw", "theta",  "pi_"),median, ~quantile(.x, probs = c(0.025, 0.975))) |>
          dplyr::mutate(type="Posterior distribution", median = median, `2.5%` = `2.5%`, `97.5%`=`97.5%`)

        return(list(g1, g2, post_param))
}

#' plot_compare_sampling_GE function to visualize the result of the model for
#' different sampling distributions
#'
#' @param fit1 model fit for poisson distribution
#' @param fit2 model fit for quasi-Poisson distribution
#' @param fit3 model fit for negative binomial distribution
#' @param data data frame with columns "weeks" and "cases_total", which are the
#' number of laboratory confirmed cases per week
#' @param GE_data list containing data about the canton of Geneva. Is included
#' as the data element "GE_data" in the package
#'
#' @return ggplot figure
#' @export
#'
#' @examples
plot_compare_sampling_GE = function (fit1, fit2, fit3, data, GE_data) {
  ggplot2::theme_set(ggplot2::theme_bw())
  cust_cols = c("chartreuse1","chartreuse3","darkgreen")

  # lab-confirmed cases
  dat_ = fit1$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975)) ) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  dat_$median = data$cases_total
  dat_$date = c(ISOweek::ISOweek2date(paste0("2020-W0",data$weeks[1], "-1")), ISOweek::ISOweek2date(paste0("2020-W",data$weeks[2:45], "-1"))) # as.Date( paste0("2020",data$weeks) , "%Y%w")

  post_1 = fit1$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975)) ) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  post_1$date = dat_$date

  post_2 = fit2$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975)) ) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  post_2$date = dat_$date
  post_3 = fit3$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975)) ) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable",
                                                  "time", "null"))
  post_3$date = dat_$date

  post_1$type = "Poisson"
  post_2$type = "Quasi-Poisson"
  post_3$type = "Negative-binomial"

  post_all = bind_rows(post_1,post_2,post_3)

  # cumulative infected
  popsize = sum(GE_data[[4]])
  serop_date = tibble(pos=NA,tested=NA,date=NA)
  for(i in 1:2) {
    serop_date[i,"pos"] = sum(GE_data[[i]][[1]]$num_pos_tests)
    serop_date[i,"tested"] = sum(GE_data[[i]][[1]]$num_tested)
    serop_date[i,"date"] = mean(c(lubridate::ymd(GE_data[[i]][[2]]),lubridate::ymd(GE_data[[i]][[3]])))
  }
  serop_date$prop = serop_date$pos/serop_date$tested

  # test correction for plotting
  sens = 0.93
  spec = 1.0

  serop_date$cor_prop <- (as.numeric(serop_date$prop) +spec -1) /(sens + spec -1)

  cum1 = fit1$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::filter(d2=="4]") %>%
    dplyr::mutate(date=dat_$date,
                  type="Poisson",
                  cum_inf_prop = median/popsize,
                  cum_inf_prop_lwb = q5/popsize,
                  cum_inf_prop_upb = q95/popsize)
  cum2 = fit2$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::filter(d2=="4]") %>%
    dplyr::mutate(date=dat_$date,
                  type="Quasi-Poisson",
                  cum_inf_prop = median/popsize,
                  cum_inf_prop_lwb = q5/popsize,
                  cum_inf_prop_upb = q95/popsize)
  cum3 = fit3$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::filter(d2=="4]") %>%
    dplyr::mutate(date=dat_$date,
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
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=type), alpha=.5) +
    geom_line(aes(y=median,colour=type)) +
    geom_point(data=dat_,aes(y=median,shape="Confirmed cases (Geneva)")) +
    geom_point(data=serop_date,aes(y=cor_prop*scaling,shape="Seroprevalence (Geneva)")) +
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

}

#' plot_unstrat_benchmark function that visualizes the differences in accuaracy
#' and performance of the various model implementations
#'
#' @param summ a dataframe with all summary statistics of the model runs from
#' Ubelix
#'
#' @return ggplot figure
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
    labs(x="ESS per second", y="Width of 95% CrI" , colour=NULL, fill=NULL, shape="Knot sequence") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position=c(.25,.4),
          legend.background = element_blank())
  g = cowplot::plot_grid(g1,g2,nrow=1,labels=c("E","F"))

  return(g)
}

#' plot_unstrat_model_fit function for the visual comparison of the model fits to simulated data
#' using the three different methods for implementing time-varying transmission.
#'
#' @param fit1 model result of Brownian motion based time-varying transmission
#' @param fit2 model result of B-spine based time-varying transmission
#' @param fit3 model result of Gaussian processes based time-varying transmission
#' @param data_cases vector with simulated laboratory confirmed case data per week
#' @param data_prob vector with simulated rho(t) value per week
#' @param params list of parameters used for simulation of data
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

  pred_1 = fit1$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Brownian motion")
  pred_2 = fit2$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="B-splines")
  pred_3 = fit3$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Gaussian process")
  pred_all = dplyr::bind_rows(pred_1, pred_2, pred_3) %>%
    dplyr::mutate(time=as.numeric(time),
                  type=factor(type,levels=c("Brownian motion","B-splines","Gaussian process")))

  g1 = pred_all %>%
    ggplot2::ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = type), alpha = 0.3) +
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
  pred_BM_prob = fit1$samples_posterior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Brownian motion")
  pred_spline_prob = fit2$samples_posterior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="B-splines")
  pred_GP_prob = fit3$samples_posterior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    tidyr::separate(variable, "\\[|\\]", into = c("variable", "time", "null")) %>%
    dplyr::mutate(type="Gaussian process")
  prob_all = dplyr::bind_rows(pred_BM_prob, pred_spline_prob, pred_GP_prob) %>%
    dplyr::mutate(time=as.numeric(time),
                  type=factor(type,levels=c("Brownian motion","B-splines","Gaussian process")))

  g2 = prob_all  %>%
    ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`,fill=type), alpha=.5) +
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
  asc1 = fit1$samples_posterior$summary(c("pi_"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    dplyr::mutate(type="Brownian motion",
                  variable = "Ascertainment rate")
  asc2 = fit2$samples_posterior$summary(c("pi_"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    dplyr::mutate(type="B-splines",
                  variable = "Ascertainment rate")
  asc3 = fit3$samples_posterior$summary(c("pi_"), median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    dplyr::mutate(type="Gaussian process",
                  variable = "Ascertainment rate")
  asc_all = dplyr::bind_rows(asc1, asc2, asc3) %>%
    dplyr::mutate(time=rep(c("Week 0 to 19", "Week 20 to 45"), 3),
                  type=factor(type,levels=c("Brownian motion","B-splines","Gaussian process")))

  g3 = asc_all %>%
    left_join(dplyr::select(dat_asc,time,true=median),by = join_by(time)) %>%
    ggplot(aes(x=as.factor(time))) +
    geom_pointrange(aes(y=median,ymin=`2.5%`,ymax=`97.5%`,colour=type),
                    position=position_dodge(.8)) +
    geom_point(aes(y=true,shape="True value",group=type),
               position=position_dodge(.8)) +
    scale_shape_manual(values=c("True value"=4)) +
    scale_colour_manual(values=cust_cols2,guide="none") +
    scale_y_continuous(labels=scales::percent) +
    theme(legend.position=c(.3,.9),
          legend.background = element_blank(),
          axis.text.y = element_text(angle=90,hjust=.5)) +
    labs(shape=NULL,x=NULL,y= "Ascertainment rate") +#expression(pi[i])) +
    coord_flip() + ylim(c(0,0.75))

  gx = prob_all  %>%
    ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=type), alpha=.5) +
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
  dat_ = rbind( data.frame(variable = "confirmed_cases_predicted",
                           time = 1:dim(sim_data)[1],
                           age_group = "0-19 years",
                           median = sim_data[,1] ),
                data.frame(variable = "confirmed_cases_predicted",
                           time = 1:dim(sim_data)[1],
                           age_group = "20-64 years",
                           median = sim_data[,2] ),
                data.frame(variable = "confirmed_cases_predicted",
                           time = 1:dim(sim_data)[1],
                           age_group = "65+ years",
                           median = sim_data[,3] ))

  post_ = fit$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>% mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-64 years", "65+ years")))
  prior_ = fit$samples_prior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time"))

  # cumulative infected
  comp_4 <- c(ind(4,1,3), ind(4,2,3), ind(4,3,3))
  cum_ = fit$samples_posterior$summary("y") %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::mutate(comp = as.numeric( str_sub(d2, end = -2) ) ) %>%
    dplyr::filter(comp %in% comp_4 ) %>%
    dplyr::mutate( age_group = ifelse(comp==comp_4[1], "0-19 years", ifelse( comp==comp_4[2], "20-64 years", "65+ years") ),
                   time = rep(1:dim(sim_data)[1] , 3),
                   cum_inf_prop = ifelse(comp == comp_4[1], median/popsize[1], ifelse(comp == comp_4[2], median/popsize[2], median/popsize[3]) ),
                   cum_inf_prop_lwb = ifelse(comp == comp_4[1], q5/popsize[1], ifelse(comp == comp_4[2], q5/popsize[2], q5/popsize[3]) ),
                   cum_inf_prop_upb = ifelse(comp == comp_4[1], q95/popsize[1], ifelse(comp == comp_4[2], q95/popsize[2], q95/popsize[3]) ) )

  data(simulated_stratified)
  n_infected_survey_strat = simulated_stratified[[3]]
  sim_data_3 = simulated_stratified[[1]]

  serop_date = data.frame(time = c(rep(25, 3), rep(45,3)), age_group = c("0-19 years","20-64 years", "65+ years"))
  serop_date$prop = c( n_infected_survey_strat[,1]/params$n_tested_survey,
                       n_infected_survey_strat[,2]/params$n_tested_survey )

  scaling = 10000
  g1_A = ggplot() +
    geom_ribbon(data=post_ , aes(x=as.numeric(time),ymin=`2.5%`,ymax=`97.5%`,fill=age_group),alpha=.3) +
    geom_line(data=post_ , aes(x=as.numeric(time),y=median,colour=age_group)) +
    geom_point(data=dat_, aes(x=as.numeric(time),y=median,shape="Confirmed cases (simulated)"), colour = "black" ) +
    geom_ribbon(data=cum_, aes(ymin=cum_inf_prop_lwb*scaling, ymax=cum_inf_prop_upb*scaling, x = time),fill="grey50",alpha=.3) +
    geom_line(data=cum_, aes(y=cum_inf_prop*scaling, x=time), linetype=2, colour="grey50") +
    facet_wrap(~age_group,ncol=1) +
    geom_point(data = serop_date, aes(x =time, y = prop*scaling,shape="Seroprevalence (simulated)")) +
    # add cumulative cases and the seroprevalence value
    scale_colour_manual(values=cust_cols) +
    scale_fill_manual(values=cust_cols) +
    scale_shape_manual(values=c("Confirmed cases (simulated)"=4,"Seroprevalence (simulated)"=3)) +
    labs(x="Time (weeks)",y="Laboratory-confirmed cases",shape=NULL,fill=NULL,colour=NULL) +
    theme(text=element_text(size=16), legend.position=c(.45,.95),legend.background = element_blank()) +
    scale_y_continuous(
      "Laboratory-confirmed cases",
      sec.axis = sec_axis(~ . /scaling, name = "Cumulative infected",labels=scales::percent)) +
    guides(colour=FALSE, fill=FALSE)

  dat_prob = rbind( data.frame(variable = "rho_simulated",
                               time = 1:dim(sim_data)[1],
                               age_group = "0-19 years",
                               median = transmission_prob_sim[,1] ),
                    data.frame(variable = "rho_simulated",
                               time = 1:dim(sim_data)[1],
                               age_group = "20-64 years",
                               median = transmission_prob_sim[,2] ),
                    data.frame(variable = "rho_simulated",
                               time = 1:dim(sim_data)[1],
                               age_group = "65+ years",
                               median = transmission_prob_sim[,3] ))

  post_prob = fit$samples_posterior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-64 years", "65+ years")))
  prior_prob = fit$samples_prior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time"))%>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-64 years", "65+ years")))

  # plot of probability
  g1_B = ggplot() +
    geom_ribbon(data=post_prob , aes(x=as.numeric(time),ymin=`2.5%`,ymax=`97.5%`,fill=age_group),alpha=.3) +
    geom_ribbon(data=prior_prob , aes(x=as.numeric(time),ymin=`2.5%`,ymax=`97.5%`),alpha=.3,fill="grey80") +
    geom_line(data=post_prob , aes(x=as.numeric(time),y=median,colour=age_group)) +
    geom_point(data=dat_prob, aes(x=as.numeric(time),y=median, shape= "True value"),colour="black" ) +
    facet_wrap(~age_group,ncol=1) +
    # add cumulative cases and the seroprevalence value
    scale_colour_manual(values=cust_cols,guide="none") +
    scale_shape_manual(values=c("True value"=4)) +
    scale_fill_manual(values=cust_cols,guide="none") +
    labs(x="Time (weeks)",y=expression(rho(t)),shape=NULL,fill=NULL,colour=NULL) +
    theme(text=element_text(size=16), legend.position=c(.25,.95), legend.background = element_blank())

  # plot of ascertainment rate
  # ascertainment
  dat_asc = data.frame(variable = "Transmission_probability",
                       time = rep(c("Week 0 to 19", "Week 20 to 45"), each = 3) ,
                       age_group = rep(c("0-19 years", "20-64 years", "65+ years"), 2),
                       null = "", median = c(params$p_detect1, params$p_detect2)) %>%
    mutate(link = paste(time, age_group))
  asc = fit$samples_posterior$summary(c("pi_")) %>%
    dplyr::mutate(time =rep(c("Week 0 to 19", "Week 20 to 45"), each = 3),
                  age_group = rep(c("0-19 years", "20-64 years", "65+ years"), 2),
                  link = paste(time, age_group))

  g3 = asc %>%
    left_join(dplyr::select(dat_asc,time,true=median, link)) %>%
    ggplot(aes(x=as.factor(time))) +
    geom_pointrange(aes(y=median,ymin=q5,ymax=q95,colour=age_group),
                    position=position_dodge(.8)) +
    geom_point(aes(y=true,shape="True value",group=age_group),
               position=position_dodge(.8)) +
    scale_shape_manual(values=c("True value"=4)) +
    scale_colour_manual(values=cust_cols) +
    scale_y_continuous(labels=scales::percent) +
    theme(legend.position=c(0.15, 0.65),
          legend.title = element_blank(),
          legend.background = element_blank(),
          axis.text.y = element_text(angle=90,hjust=.5),text=element_text(size=16)) +
    labs(shape=NULL,x=NULL,y=expression(pi[i]))

  full_plot <- cowplot::plot_grid(g1_A,g1_B,nrow=1,labels=c("A","B"), rel_widths = c(1.25,1))
  full_plot2 <- cowplot::plot_grid(full_plot,g3,nrow=2,labels=c("","C"), rel_heights = c(3,1))

  return(full_plot2)
}


#' plot_fitsim_strat_GE function to visualize the model fit to the stratified
#' data of the canton of Geneva
#'
#' @param fit model fit
#' @param data dataframe with weekly laboratory confirmed cases per week per age-class
#' @param params list of parameters fixed during the analyses
#'
#' @return
#' @export
#'
#' @examples
plot_fitsim_strat_GE <- function(fit, data, params){

  # aesthetics
  ggplot2::theme_set(ggplot2::theme_bw())
  cust_cols = c("cyan3","pink","purple")

  # data from stratified simulation
  popsize = params$popdist
  n_tested = params$n_tested_survey

  dates = c(ISOweek::ISOweek2date(paste0("2020-W0",data$weeks[1], "-1")), ISOweek::ISOweek2date(paste0("2020-W",data$weeks[2:45], "-1")))

  # posterior predictive check
  dat_ = rbind( data.frame(variable = "GE_data",
                           time = dates,
                           age_group = "0-19 years",
                           median = data[,2] ) ,
                data.frame(variable = "GE_data",
                           time = dates,
                           age_group = "20-64 years",
                           median = data[,3] ) ,
                data.frame(variable = "GE_data",
                           time = dates,
                           age_group = "65+ years",
                           median = data[,4] ) )

  post_ = fit$samples_posterior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-64 years", "65+ years")))
  post_$date = rep(dates, each = 3)
  prior_ = fit$samples_prior$summary(c("confirmed_cases_predicted"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-64 years", "65+ years")))
  prior_$date = rep(dates, each = 3)

  # cumulative infected
  comp_4 <- c(ind(4,1,3), ind(4,2,3), ind(4,3,3))
  cum_ = fit$samples_posterior$summary("y", median, ~quantile(.x, probs = c(0.025, 0.975))) %>%
    tidyr::separate(variable,sep=",",into=c("d1","d2")) %>%
    dplyr::mutate(comp = as.numeric( str_sub(d2, end = -2) ) ) %>%
    dplyr::filter(comp %in% comp_4 ) %>%
    dplyr::mutate( age_group = ifelse(comp==comp_4[1], "0-19 years", ifelse( comp==comp_4[2], "20-64 years", "65+ years") ),
                   time = rep(data$date, 3),
                   cum_inf_prop = ifelse(comp == comp_4[1], median/popsize[1], ifelse(comp == comp_4[2], median/popsize[2], median/popsize[3]) ),
                   cum_inf_prop_lwb = ifelse(comp == comp_4[1], `2.5%`/popsize[1], ifelse(comp == comp_4[2], `2.5%`/popsize[2], `2.5%`/popsize[3]) ),
                   cum_inf_prop_upb = ifelse(comp == comp_4[1], `97.5%`/popsize[1], ifelse(comp == comp_4[2], `97.5%`/popsize[2], `97.5%`/popsize[3]) ) )

  cum_$time <- rep(dates , 3)

  n_infected_survey_strat = params$n_infected_survey

  # cumulative infected
  serop_date = tibble(pos=NA,tested=NA,date=NA, age_group = NA, prop = NA)
  age_groups = c("0-19 years", "20-64 years", "65+ years")
  check_group = c("group1", "group2", "group3")
  for(j in 1:3) { # for loop over age:groups
    for (i in 1:2){ # for loop over seroprevalence studies
      if ( dim(GE_data[[i]][[1]] %>% filter(age_group==check_group[j]))[1]>0 ){
        serop_date[(j-1)*2 +i,"pos"] = GE_data[[i]][[1]] %>% filter(age_group==check_group[j]) %>% select(num_pos_tests)
        serop_date[(j-1)*2 +i,"tested"] = GE_data[[i]][[1]]%>% filter(age_group==check_group[j]) %>% select(num_tested)
        serop_date[(j-1)*2 +i,"date"] = mean(c(lubridate::ymd(GE_data[[i]][[2]]),lubridate::ymd(GE_data[[i]][[3]])))
        serop_date[(j-1)*2 +i, "prop"] = ( serop_date[(j-1)*2 +i,"pos"]/serop_date[(j-1)*2 +i,"tested"] + params$spec -1 ) / (params$sens + params$spec -1)
      } else {
        serop_date[(j-1)*2 +i,"pos"] = NA
        serop_date[(j-1)*2 +i,"tested"] = NA
        serop_date[(j-1)*2 +i,"date"] = NA
        serop_date[(j-1)*2 +i, "prop"] = NA
      }
      serop_date[(j-1)*2 +i, "age_group"] =  age_groups[j]
    }
  }

  scaling = 20000
  g1_A = ggplot() +
    geom_ribbon(data=post_ , aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=age_group),alpha=.3) +
    geom_line(data=post_ , aes(x=date,y=median,colour=age_group)) +
    geom_point(data=dat_, aes(x=time,y=median,shape="Confirmed cases (Geneva)"),colour="black" ) +
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
    geom_vline(xintercept = dat_$time[10], color = "grey70") +
    guides(colour=FALSE, fill=FALSE)

  post_prob = fit$samples_posterior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time")) %>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-64 years", "65+ years")))
  post_prob$date = rep(dates,each = 3)
  prior_prob = fit$samples_prior$summary(c("rho"), median, ~quantile(.x, probs = c(0.025, 0.975))) |>
    tidyr::separate(variable,"\\[|\\]",into=c("variable","step", "NULL")) |>
    tidyr::separate(step,",",into=c("age_group","time"))%>%
    mutate(age_group = ifelse(age_group==1, "0-19 years", ifelse(age_group==2, "20-64 years", "65+ years")))
  prior_prob$date = rep(dates,each = 3)

  # plot of probability
  g1_B = ggplot() +
    #geom_ribbon(data=prior_prob , aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill="Prior"),alpha=.3) +
    geom_ribbon(data=post_prob , aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=age_group),alpha=.3) +
    geom_line(data=post_prob , aes(x=date,y=median,colour=age_group)) +
    facet_wrap(~age_group,ncol=1) +
    geom_vline(xintercept = dat_$time[10], color = "grey70") +
    # add cumulative cases and the seroprevalence value
    scale_colour_manual(values=c(cust_cols, "grey80"),guide="none") +
    scale_fill_manual(values=c(cust_cols, "grey80"),guide="none") +
    labs(x=element_blank(),y=expression(rho(t)),shape=NULL,fill=NULL,colour=NULL) +
    theme(text=element_text(size=16), legend.position=c(.25,.25), legend.background = element_blank(),
          legend.title = element_blank() )

  # plot of ascertainment rate
  # ascertainment
  asc = fit$samples_posterior$summary(c("pi_")) %>%
    tidyr::separate(variable,sep=c(4,5,6,7),into=c("variable", "age_group", "d1" ,"Time period")) %>%
    dplyr::mutate( age_group = ifelse(age_group == 1,  "0-19 years", ifelse(age_group==2, "20-64 years", "65+ years") ),
                   date = ifelse(`Time period`==1, "Spring 2020" , "Fall/Winter 2020" ))

  asc$date = factor(asc$date, levels = c("Spring 2020", "Fall/Winter 2020"))

  g3 = asc %>% ggplot(aes(x=as.factor(date))) +
    geom_pointrange(aes(y=median,ymin=q5,ymax=q95,colour=age_group),
                    position=position_dodge(.8)) +
    scale_colour_manual(values=cust_cols) +
    scale_y_continuous(labels=scales::percent) +
    theme(legend.position=c(0.15, 0.85),
          legend.background = element_blank(),
          axis.text.y = element_text(angle=90,hjust=.5),text=element_text(size=16),
          legend.title = element_blank()) +
    labs(shape=NULL,x=NULL,y=expression(pi[i]))

  full_plot <- cowplot::plot_grid(g1_A,g1_B,nrow=1,labels=c("A","B"), rel_widths = c(1.25,1))
  full_plot2 <- cowplot::plot_grid(full_plot,g3,nrow=2,labels=c("","C"), rel_heights = c(3,1))

  return(full_plot2)
}


#' plot_single_benchmark function to visualize the result on the model runs for
#' a single time-varyiation method
#'
#' @param summ dataframe with results of individual model runs from Ubelix
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
    theme(legend.position="none",  #c(.38,.9),
          legend.background = element_blank(),text=element_text(size=16), strip.text.y = element_blank() ) +
    guides(color=guide_legend(nrow=2, byrow=TRUE))

  g2 = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = CI_size_prob, colour = solver , shape = warmup_iter),size=2) +
    scale_colour_manual(values = cust_cols2, labels=c("0"="rk45","1"="Adams","2"="Bdf", "3"="ckrk","4"="Trapeziodal"))+
    facet_grid(rows = vars(Tolerance), labeller = label_both ) +
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes , guide = "none") +
    labs(x="ESS per second", y="Sharpness" , colour="ODE solver", fill=NULL, shape="Warmup iterations") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position="none", #c(.28,.86),
          legend.background = element_blank(), text=element_text(size=16))

  g3 = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = max_rhat, colour = solver , shape = warmup_iter),size=2) +
    scale_colour_manual(values = cust_cols2, labels=c("0"="rk45","1"="Adams","2"="Bdf", "3"="ckrk","4"="Trapeziodal"))+
    facet_grid(rows = vars(Tolerance), labeller = label_both ) +
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes , guide = "none") +
    labs(x="ESS per second", y="Sharpness" , colour="ODE solver", fill=NULL, shape="Warmup iterations") +
    # guides(shape=guide_legend(ncol=2)) +
    theme(legend.position="none", #c(.28,.86),
          legend.background = element_blank(), text=element_text(size=16))

  gx = summ %>%
    ggplot() +
    geom_point(aes(x = 1/time_per_ESS, y = RMSE_prob, colour = solver, shape = warmup_iter),size=2) +
    scale_colour_manual(values = cust_cols2, guide = "none", labels=c("0"="rk45","1"="Adams","2"="Bdf", "3"="ckrk","4"="Trapeziodal"))+
    facet_grid(rows = vars(Tolerance), labeller = label_both ) +
    scale_x_continuous(trans='log10') +
    scale_shape_manual(values=cust_shapes ) +
    labs(x="ESS per second", y="Error" , colour="ODE solver", fill=NULL, shape="Warmup iterations") +
    theme(legend.position="bottom",
          legend.background = element_blank(),text=element_text(size=16), strip.text.y = element_blank() ) +
    guides(color=guide_legend(nrow=2, byrow=TRUE))

  leg = ggpubr::get_legend(gx)

  g = cowplot::plot_grid(g1,g2, nrow=1,labels=c("A","B"))
  gfull = cowplot::plot_grid(g, leg, nrow=2, rel_heights = c(6,1))

  return(gfull)
}
