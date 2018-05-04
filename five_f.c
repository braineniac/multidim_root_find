#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include "five_f.h"

int five_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	double a = ((struct five_f_params *) params)->a;
	double b = ((struct five_f_params *) params)->b;
	double c = ((struct five_f_params *) params)->c;
	double d = ((struct five_f_params *) params)->d;
	double e = ((struct five_f_params *) params)->e;

	const double x0 = gsl_vector_get(x,0);
	const double x1 = gsl_vector_get(x,1);
	const double x2 = gsl_vector_get(x,2);
	const double x3 = gsl_vector_get(x,3);
	const double x4 = gsl_vector_get(x,4);

	const double y0 = a*x0*x1 + c*x4 + d*x4*x2;
	const double y1 = b*x1 -d*x4*x3 + e*x2;
	const double y2 = a*b*x2*x3 + d*b*x1 - x4*x2;
	const double y3 = b*d*x1*x2 - x2 + x3*x4;
	const double y4 = a*x2 - c*e*x4*x1;

	gsl_vector_set(f,0,y0);
	gsl_vector_set(f,1,y1);
	gsl_vector_set(f,2,y2);
	gsl_vector_set(f,3,y3);
	gsl_vector_set(f,4,y4);

	return GSL_SUCCESS;

}
