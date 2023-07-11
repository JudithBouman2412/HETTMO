// compute boundary value
real L_f1 = c_f1*num_t;

// compute basis functions for f1
matrix[num_t, M_f1] PHI_f1 = PHI_EQ(num_t, M_f1, L_f1, ts); // use transformed ts data instead of the real ts
