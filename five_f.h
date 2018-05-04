#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

struct five_f_params
{
	double a;
	double b;
	double c;
	double d;
	double e;

} five_f_params;

int five_f  (const gsl_vector * x,  void * params, gsl_vector * f);
