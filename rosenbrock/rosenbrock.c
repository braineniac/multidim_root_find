#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

int run_rosenbrock(const gsl_multiroot_fsolver_type *T, double * x_init, struct rosenbrock_params p)
{
	gsl_multiroot_fsolver *s;

	int status;
	size_t iter = 0;

	const size_t n =2;
	gsl_multiroot_function f = {&rosenbrock_f, n, &p};

	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set(x, 0, x_init[0]);
	gsl_vector_set(x,1, x_init[1]);

	s = gsl_multiroot_fsolver_alloc (T,2);
	gsl_multiroot_fsolver_set(s, &f,x);

	print_state_f(iter,s);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);

		print_state_f(iter,s);

		if (status)
			break;

		status = gsl_multiroot_test_residual(s->f, 1e-7);
	}

	while (status == GSL_CONTINUE && iter <1000);

	printf("status =%s\n", gsl_strerror(status));

	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);

	return 0;
}


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
