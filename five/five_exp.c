#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include "five.h"

int five_exp_f (const gsl_vector * x, void * params, gsl_vector * f) 
{

	double a = ((struct five_params *) params)->a;
	double b = ((struct five_params *) params)->b;
	double c = ((struct five_params *) params)->c;
	double d = ((struct five_params *) params)->d;
	double e = ((struct five_params *) params)->e;

	const double x0 = gsl_vector_get(x,0);
	const double x1 = gsl_vector_get(x,1);
	const double x2 = gsl_vector_get(x,2);
	const double x3 = gsl_vector_get(x,3);
	const double x4 = gsl_vector_get(x,4);

	const double y0 = c*exp(x1) + b*x0 + d*x3 +e*x4;
	const double y1 = b*x1 -d*x4 + e*x3;
	const double y2 = a*x0 - b*x3 - x2;
	const double y3 = b*d*x4 - c*x2;
	const double y4 = a*e*x0 - b*x3;

	gsl_vector_set(f,0,y0);
	gsl_vector_set(f,1,y1);
	gsl_vector_set(f,2,y2);
	gsl_vector_set(f,3,y3);
	gsl_vector_set(f,4,y4);

	return GSL_SUCCESS;

}


