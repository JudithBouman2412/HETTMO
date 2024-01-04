# HETTMO
HETerogeneous Transmission MOdel -- R package

This package contains a Bayesian workflow for modelling respiratory infectious diseases with time-varying transmission in Stan. 

It includes three methods for implementing the time-varying transmission: Brownian motion, B-splines and Gaussian processes. 

The package relies on functionality of cmdstanR.

More information can be found in [this pre-print](https://doi.org/10.1101/2023.10.09.23296742).


TO DO:

add instructions on how to change the age-groups. 

add instructions on how to adapt the modelled ODE-system. 

## Install package

To install the package directly from Github in R, you can use:

library(devtools)

install_github("JudithBouman2412/HETTMO")

## Generalizability of code 

### SEIR backbone 

Here, we indicate the parts of the code that need to be adjusted to change the compartmental model that is the backbone of the workflow. 

For the unstratified B-spline based model, adjust:
+ the ODE-system: [for standard ODE-solvers](https://github.com/JudithBouman2412/HETTMO/blob/83e897ffa062ee2270ec68342bc725e4a75ed59a/inst/functions/SEIR_ode_spline.stan#L1C1-L33C1), and [for the trapezoidal solver](https://github.com/JudithBouman2412/HETTMO/blob/83e897ffa062ee2270ec68342bc725e4a75ed59a/inst/functions/SEIR_ode_spline.stan#L78C1-L138C2). 

For the stratified B-spline based model, adjust:
+ the ODE-system: [for standard ODE-solvers](https://github.com/JudithBouman2412/HETTMO/blob/83e897ffa062ee2270ec68342bc725e4a75ed59a/inst/functions/SEIR_ode_spline.stan#L34C1-L76C2), and [for the trapezoidal solver](https://github.com/JudithBouman2412/HETTMO/blob/83e897ffa062ee2270ec68342bc725e4a75ed59a/inst/functions/SEIR_ode_spline.stan#L140C1-L202C2)

For the unstratified aGP based model, adjust:
+ the ODE-system: [for standard ODE-solvers](https://github.com/JudithBouman2412/HETTMO/blob/83e897ffa062ee2270ec68342bc725e4a75ed59a/inst/functions/SEIR_ode_GP.stan#L1-L23) and [for the trapezoidal solver](https://github.com/JudithBouman2412/HETTMO/blob/83e897ffa062ee2270ec68342bc725e4a75ed59a/inst/functions/SEIR_ode_GP.stan#L26C1-L89C2).

For the unstratified BM based model, adjust:
+ the ODE-system: [for standard ODE-solvers](https://github.com/JudithBouman2412/HETTMO/blob/83e897ffa062ee2270ec68342bc725e4a75ed59a/inst/functions/SEIR_ode_BM.stan#L1C1-L31C2) and for [the trapezoidal solver](https://github.com/JudithBouman2412/HETTMO/blob/83e897ffa062ee2270ec68342bc725e4a75ed59a/inst/functions/SEIR_ode_BM.stan#L33C1-L103C2). 



