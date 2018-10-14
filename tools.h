#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

int print_state_f(size_t iter, gsl_multiroot_fsolver * s);
int print_state_fdf(size_t iter, gsl_multiroot_fdfsolver *s);
int print_state_five(size_t iter, gsl_multiroot_fsolver *s);
int print_state_five_fdf(size_t iter, gsl_multiroot_fdfsolver *s);

