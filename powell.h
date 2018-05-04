#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>


struct powell_params {
	double A;

} powell_params;

int powell_f (const gsl_vector * x, void * p, gsl_vector * f);
int powell_df (const gsl_vector *x, void * p, gsl_matrix * J);
int powell_fdf (const gsl_vector *x, void *p, gsl_vector * f, gsl_matrix * J);
