vector build_b_spline(array[] real t, data vector ext_knots, data int index, data int order);
vector build_b_spline(array[] real t, data vector ext_knots, data int index, data int order){
  // INPUTS:
  //    t:          the points at which the b_spline is calculated
  //    ext_knots:  the set of extended knots
  //    index:        the indexex of the b_spline
  //    order:      the order of the b-spline
  vector[size(t)] b_spline;
  vector[size(t)] w1 = rep_vector(0, size(t));
  vector[size(t)] w2 = rep_vector(0, size(t));
  if (order==1)
    for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
      b_spline[i] = (ext_knots[index] <= t[i]) && (t[i] < ext_knots[index+1]);
  else {
    if (ext_knots[index] != ext_knots[index+order-1])
      w1 = (to_vector(t) - rep_vector(ext_knots[index], size(t))) /
           (ext_knots[index+order-1] - ext_knots[index]);
    if (ext_knots[index+1] != ext_knots[index+order])
      w2 = 1 - (to_vector(t) - rep_vector(ext_knots[index+1], size(t))) /
               (ext_knots[index+order] - ext_knots[index+1]);
    // Calculating the B-spline recursively as linear interpolation of two lower-order splines
    b_spline = w1 .* build_b_spline(t, ext_knots, index, order-1) +
               w2 .* build_b_spline(t, ext_knots, index+1, order-1);
  }
  return b_spline;
}
matrix find_support( data vector ext_knots, data int order, data vector knots, data int num_basis ){
    int num_interval = size(knots)-1;
    // define t for which polynomal is well defined for each b-spline
    matrix[num_basis, num_interval] support_bs;

    for (i in 1:num_basis){
      vector[num_interval] support = rep_vector(0, num_interval);

      for (j in 1:num_interval){
        if (knots[j]>=ext_knots[i]){
          if (knots[j+1]<=ext_knots[i+order]){
            support[j] = 1;
          }
        }
      }
      support_bs[i,] = to_row_vector(support);
    }

    return support_bs;
}
vector find_poly( data array[] real x, data row_vector y){
  // given m points x,y find polynomial of degree m-1 through these points
  int m = size(x);

  // vector initializing the polynomial
  vector[m] thepoly = rep_vector(0, m);

  for (i in 1:m){
    vector[m] theterm = rep_vector(0, m);

    real r = 1;
    for (j in 1:m){
      if (i!=j){
        r = r * (x[i]-x[j]);
      }
    }

    theterm[m] = y[i]/r;

    for (j in 1:m){
      if (i!=j){
        for (k in 2:m){
          theterm[k-1] += theterm[k];
          theterm[k] *= (-x[j]);
        }
      }
    }

    for (j in 1:m){
      thepoly[j] = thepoly[j]+theterm[j];
    }
  }

  thepoly = reverse(thepoly);

  return thepoly;
}

matrix ana_B( data matrix B, data vector ext_knots, data array[] real int_knots, data int order, data vector knots, data int num_basis, data matrix support_bs){
    // empty matrix to save estimated coefficients for each b_spline and each area at which it is defined
    int num_interval = size(knots)-1;
    matrix[num_basis, num_interval*order] b_hat = rep_matrix(0.0, num_basis, num_interval*order);

    // fit polynomal for each b_spline at the interval at which it is defined
    for (i in 1:num_basis){
      for (j in 1:num_interval){
        if (support_bs[i,j]){
          b_hat[i,((j-1)*order+1):(j*order)] = to_row_vector( find_poly(int_knots[((j-1)*(order-1)+1):((order-1)*j+1)], B[i,((j-1)*(order-1)+1):((order-1)*j+1)]) ); //extend to use any degree/order
        }
      }
    }

    return b_hat;
  }

// function to calculate full spline for new values
real analytical_bspline( data matrix b_hat, vector a, real t_new, data vector knots, data int order ) {

  int num_interval = size(knots)-1;
  int n_b_spline = dims(b_hat)[1];

  matrix[n_b_spline, order] b_hat_now;

  for (i in 1:num_interval){
    if (t_new>=knots[i]){
      if (t_new<knots[i+1]){
        b_hat_now = b_hat[,((i-1)*order+1):(i*order)];
      }
    }
  }

  vector[n_b_spline] b_spline;

  for (i in 1:n_b_spline){
    b_spline[i] = a[i] * (b_hat_now[i,1] +b_hat_now[i,2]*t_new +b_hat_now[i,3]*t_new^2 + b_hat_now[i,4]*t_new^3);
  }

  real b_spline_1 = sum(b_spline);

  return(b_spline_1);
}
