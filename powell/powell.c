#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

int run_powell(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct powell_params pw) {
	
	int status;
	size_t iter=0;

	gsl_multiroot_fdfsolver *s;

	const size_t n = 2;
	gsl_multiroot_function_fdf f = {&powell_f,&powell_df,&powell_fdf, n, &pw};

	gsl_vector *x = gsl_vector_alloc(n);

	gsl_vector_set(x,0,x_init[0]);
	gsl_vector_set(x,1,x_init[1]);

	s = gsl_multiroot_fdfsolver_alloc(T,2);
	gsl_multiroot_fdfsolver_set(s,&f,x);
	
	print_state_fdf(iter,s);

	do {
		iter++;
		status = gsl_multiroot_fdfsolver_iterate(s);

		print_state_fdf(iter,s);

		if(status)
			break;

		status = gsl_multiroot_test_residual(s->f, 1e-7);
	}

	while(status==GSL_CONTINUE && iter<1000);

	printf("status =%s\n", gsl_strerror(status));

	gsl_multiroot_fdfsolver_free(s);
	gsl_vector_free(x);
		
	return 0;
}


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
