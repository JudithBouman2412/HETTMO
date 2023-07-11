// transformed data block for splines: create b_hat

int num_basis = num_knots + spline_degree - 1; // total number of B-splines
int order = spline_degree + 1;

// create extended knots for fitting the splineC
vector[spline_degree + num_knots] ext_knots_temp;
vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));

// Create internal knots to have exactly enough points for each interval to fit a polynomial of degree 3
int num_int = (num_knots-1)*(order-1)+1;
array[num_int] real int_knots;
array[order-1] real int_now;
for (i in 2:num_knots){
  real dif = floor( (knots[i]-knots[i-1])/(order-1));
  for (j in 1:(order-1)){
    int_now[j] = knots[i-1] + (j-1)*dif;
  }
  int_knots[((i-2)*(order-1)+1):((i-1)*(order-1))] = int_now;
}
int_knots[num_int] = knots[num_knots];

// Build the spline --> this spline stays the same now
matrix[num_basis, num_int] B;

for (indicator in 1:num_basis){
  B[indicator,:] = to_row_vector(build_b_spline( int_knots, ext_knots, indicator, spline_degree + 1));
}
B[num_knots + spline_degree - 1, num_int] = 1; // --> this is in the manual, but I do not understand why...

// find support for each spline
matrix[num_basis, (num_knots-1) ] support_bs = find_support( ext_knots, order, knots, num_basis );

// find coefficients for each support of each spline
matrix[num_basis, (num_knots-1)*order] b_hat = ana_B( B, ext_knots, int_knots, order, knots, num_basis, support_bs);
