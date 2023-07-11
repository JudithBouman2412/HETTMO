// data objects specific for stratified model 
int num_age;
int num_sex;
vector[num_age*num_sex] popdist;

// contact matrix
matrix[num_age*num_sex,num_age*num_sex] contact;