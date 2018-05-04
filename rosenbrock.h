#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

struct rosenbrock_params
{       
	double a;
	double b;

} rosenbrock_params;

	
int rosenbrock_f (const gsl_vector * x, void * params, gsl_vector * f);
