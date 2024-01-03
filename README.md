# HETTMO
HETerogeneous Transmission MOdel -- R package

This package contains a Bayesian workflow for modelling respiratory infectious diseases with time-varying transmission in Stan. 

It includes three methods for implementing the time-varying transmission: Brownian motion, B-splines and Gaussian processes. 

The package relies on functionality of cmdstanR.

More information can be found in [this pre-print](https://doi.org/10.1101/2023.10.09.23296742).


TO DO:

add instructions on how to change the age-groups. 

add instructions on how to adapt the modelled ODE-system. 

## Install package

To install the package directly from Github, you can use:

library(devtools)
install_github("JudithBouman2412/HETTMO")
