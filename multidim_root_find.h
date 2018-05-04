#include <stdio.h>
#include <gsl/gsl_multiroots.h>

int run_rosenbrock(const gsl_multiroot_fsolver_type *T, double * x_init, struct rosenbrock_params p);
int run_five(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p);
int run_powell(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct powell_params pw);
