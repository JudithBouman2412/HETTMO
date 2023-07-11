#' build_b_spline
#'
#' Construct basic b spline of (ind) that indicates the values at time t for this
#' spline. Used from Milad Kharratzadeh (https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html))
#'
#' @param t time points at which to calculate the spline values
#' @param ext_knots extended series of knots
#' @param ind index over the number of basic splines
#' @param order order of the spline = degree + 1
#'
#' @return a vector of the value of the b spline of ind at times t
#' @export
#'
#' @examples
build_b_spline <-  function(t, ext_knots, ind, order) {
  b_spline = rep(NA, length(t))
  w1 = rep(0, length(t)) # weight on the first part of the interval
  w2 = rep(0, length(t)) # weight on the second part of the interval
  if (order==1){
    for (i in 1:length(t)){ ## B-splines of order 1 are piece-wise constant
      b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]) }
  } else {
    if (ext_knots[ind] != ext_knots[ind+order-1]) {
      w1 = (t - rep(ext_knots[ind], length(t))) / (ext_knots[ind+order-1] - ext_knots[ind])}
    if (ext_knots[ind+1] != ext_knots[ind+order]) {
      w2 = 1 - (t - rep(ext_knots[ind+1], length(t))) / (ext_knots[ind+order] - ext_knots[ind+1])}
    ## Calculating the B-spline recursively as linear interpolation of two lower-order splines
    b_spline = w1 * build_b_spline(t, ext_knots, ind, order-1) + w2 * build_b_spline(t, ext_knots, ind+1, order-1) }
  return(b_spline)
}

#' Create_B
#'
#' Function to run the algorithm (inspired from Milad Kharratzadeh (https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html))
#' and create matrix B to define a basic spline
#'
#' @param t time points at which to calculate the b spline
#' @param order order of the spline = degree + 1
#' @param knots set of knots
#'
#' @return matrix B with each row representing a basic spline and each column a value of t
#' @export
#'
#' @examples
create_B <- function(t, order, knots){

  degree = order - 1
  num_knots = length(knots)
  num_basis = num_knots + degree - 1

  ext_knots = c(rep(knots[1],degree), knots, rep(knots[num_knots], degree))

  B = matrix(NA, nrow = num_basis, ncol = length(t))

  for (ind in 1:num_basis){
    B[ind,] <-  build_b_spline(t, ext_knots, ind, order)
  }
  B[num_basis, length(t)] <- 1

  return(B)
}

#' ana_B
#'
#' From a matrix B, containing the values of each basic spline (rows) for the knots
#' and (degree - 1) points between the knots, this function calculates the support
#' of each basic spline and the coefficients of the polynomial of each basic spline
#' between each set of knots.
#'
#' @param B matrix with value of basic splines (rows) for given set of time points (columns)
#' @param ext_knots set of extended knots
#' @param int_knots set of time points, defined such that there are degree -1 points between
#' each consecutive knots
#' @param order degree of the basic spline polynomial + 1
#' @param knots set of knots
#' @param num_basis number of basic splines (number of knots + degree - 1)
#'
#' @return matrix b.hat with for each b_spline the coefficients of the polynomial
#' for one combination of knots
#' matrix support_bs which indicates on which knots a basic spline is defined
#' @export
#'
#' @examples
ana_B <- function(B, ext_knots, int_knots, order, knots, num_basis ){
  # t is the t used to set-up matrix B

  # empty matrix to save estimated coefficients for each b_spline and each area at which it is defined
  num_interval <- length(knots)-1
  b.hat <- matrix(NA, nrow = num_basis, ncol= order*num_interval )

  # define t for which polynomal is well defined for each t
  support_bs <- matrix(NA, nrow = num_basis, ncol = num_interval)
  for (i in 1:num_basis){
    support <- c()
    for (j in 1:num_interval){
      support <- c(support,  ((knots[j]>=ext_knots[i]) & (knots[j+1]<=ext_knots[i+order]))  )
    }
    support_bs[i,] <- support
  }

  # fit polynomal for each b_spline at the interval at which it is defined
  for (i in 1:num_basis){
    for (j in 1:num_interval){
      if (support_bs[i,j]){
        t_support <- knots[j]<=int_knots & knots[j+1]>=int_knots
        b.hat[i,((j-1)*order+1):(j*order)] <- stats::coef(stats::lm( B[i,t_support]~ int_knots[t_support] + I(int_knots[t_support]^2)+ I(int_knots[t_support]^3))) #extend to use any degree/order
      }
    }
  }

  return(list(b.hat, support_bs))
}


#' analytical_bspline
#'
#'function to calculate full spline for new values, given the coefficients of the
#'polynomials that define the spline (b_hat) the coefficients of the spline (a)
#'and the support of each polynomial (b_support)
#'
#' @param b_hat matrix with coefficients of the b_spline for each interval between knots
#' @param a coefficients of the b spline
#' @param b_support matrix indicationg for which interval between knots a b spline is defined
#' @param t_new the value for which the b spline is to be calculated
#' @param knots set of knots
#'
#' @return spline value at t_new
#' @export
#'
#' @examples
analytical_bspline <-  function( b_hat, a, b_support, t_new, knots ) {

  num_interval <- length(knots)-1
  order = 4

  interval = c()
  for (i in 1:num_interval){
    interval = c(interval,
                 t_new>=knots[i] & t_new<knots[i+1] )
  }

  interval.int <- (1:num_interval)[interval]

  b.hat.now <- b_hat[,(( interval.int -1)*order+1):( interval.int *order)]

  b_spline_full <- a*apply(b.hat.now, MARGIN = 1, FUN = function(x, t_new) {
    x[1] + x[2]*t_new + x[3]*t_new^2 + x[4]*t_new^3    }, t_new=t_new )  # extend to use any degree/order

  b_spline_out <- sum(b_spline_full, na.rm = TRUE)

  return(b_spline_out)
}

#' find_poly
#'
#' Use the Lagrange algorithm to fit a polynomial of degree (number of points -1)
#' through a set of points.
#' Code inspired by: https://en.wikibooks.org/wiki/Algorithm_Implementation/Mathematics/Polynomial_interpolation#:~:text=Lagrange%20interpolation%20is%20an%20algorithm,xi%2C%20yi).
#'
#' @param x x values
#' @param y y values
#'
#' @return coefficients of the fitted polynomial
#' @export
#'
#' @examples
find_poly <- function(x,y){

  m = length(x)

  thepoly <- rep(0,m)

  for (i in 1:m){
    theterm <- rep(0,m)
    r <- 1
    for (j in 1:m){
      if (i!=j){
        r <- r * (x[i]-x[j])
      }
      tmp <- y[i]/r
    }

    theterm[1] = tmp

    theterm <- rev(theterm)

    for (j in 1:m){
      if (i!=j){
        for (k in 2:(m)){
          theterm[k-1] <- theterm[k-1]+theterm[k]
          theterm[k] <- theterm[k]*(-x[j])
        }
      }
    }

    for (j in 1:m){
      thepoly[j] <- thepoly[j]+theterm[j]
    }
  }

  thepoly <- rev(thepoly)

  return(thepoly)
}

#' probability_transmission_spline
#'
#' Function that produces the probability of transmission for times t for a spline with
#' given degree, knots and coefficients a
#'
#' @param t time
#' @param order degree of the basic spline polynomial + 1
#' @param knots set of knots
#' @param a coefficients for the spline
#' @param beta intercept of the spline function
#'
#' @return
#' @export
#'
#' @examples
probability_transmission_spline <- function(t, order, knots, a, beta){

  B = create_B(t, order, knots)

  a_tot <- c(beta, a)

  prob <- boot::inv.logit(a_tot%*%B)

  return(prob[1,])
}



