#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include "rosenbrock.h"


int rosenbrock_f (const gsl_vector * x, void * params, gsl_vector * f)
{       
	double a = ((struct rosenbrock_params *) params)->a;
	double b = ((struct rosenbrock_params *) params)->b;
	
	const double x0 = gsl_vector_get(x,0);
	const double x1 = gsl_vector_get(x,1);
	
	const double y0 = a * (1 - x0);
	const double y1 = b* (x1-x0 *x0);
	
	gsl_vector_set(f,0,y0);
	gsl_vector_set(f,1,y1);
	
	return GSL_SUCCESS;
}
