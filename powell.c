#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include "powell.h"

int powell_f (const gsl_vector * x, void * p, gsl_vector * f) {
	 struct powell_params * params = (struct powell_params *)p;
	 const double A = (params->A);
	 const double x0 = gsl_vector_get(x,0);
	 const double x1 = gsl_vector_get(x,1);
	 
	 gsl_vector_set (f, 0, A * x0 * x1 - 1);
	 gsl_vector_set (f, 1, (exp(-x0) + exp(-x1) - (1.0 + 1.0/A)));
	 return GSL_SUCCESS;
}

int powell_df (const gsl_vector *x, void * p, gsl_matrix * J) {
	const double x0 = gsl_vector_get(x,0);
	const double x1 = gsl_vector_get(x,1);
	struct powell_params * params = (struct powell_params *) p;
	const double A = (params->A);
	gsl_matrix_set(J,0,0, A*x1);
	gsl_matrix_set(J,0,1, A*x0);
	gsl_matrix_set(J,1,0, -exp(-x0));
	gsl_matrix_set(J,1,1, -exp(-x1));
	return GSL_SUCCESS;
}

int powell_fdf (const gsl_vector *x, void *p, gsl_vector * f, gsl_matrix * J) {
	struct powell_params * params = (struct powell_params *) p;
	const double A = (params->A);
	const double x0 = gsl_vector_get(x,0);
	const double x1 = gsl_vector_get(x,1);
	
	const double u0 = exp(-x0);
	const double u1 = exp(-x1);
	
	gsl_vector_set (f, 0, A * x0 * x1 - 1);
	gsl_vector_set (f, 1, u0 + u1 - (1 + 1/A));
	
	gsl_matrix_set (J, 0, 0, A * x1);
	gsl_matrix_set (J, 0, 1, A * x0);
	gsl_matrix_set (J, 1, 0, -u0);
	gsl_matrix_set (J, 1, 1, -u1);
	return GSL_SUCCESS;
}


