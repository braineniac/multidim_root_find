#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

struct five_params
{
	double a;
	double b;
	double c;
	double d;
	double e;

} five_params;

int five_f  (const gsl_vector * x,  void * params, gsl_vector * f);
int five_df (const gsl_vector *x, void * p, gsl_matrix * J);
int five_fdf (const gsl_vector *x, void *p, gsl_vector * f, gsl_matrix * J);
