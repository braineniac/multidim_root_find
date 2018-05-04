#include <stdio.h>
#include <string.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include "powell.c"
#include "rosenbrock.c"
#include "tools.c"

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

int run_powell(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct powell_params pw)
{
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

	do
	{
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

int main(int argc, char *argv[]) 
{
	double x_init_rosenbrock[2] = {-10.0,-5.0};
	double x_init_powell[2] = {10.0,-20.0};
	struct rosenbrock_params p_rosenbrock = {.a=1.0, .b=10};
	struct powell_params p_powell = {.A=10000};
	
	if (argc==2 && strcmp(argv[1],"0") == 0 )
	{
		printf("Running Rosenbrock with the GSL newton solver!\n");
		run_rosenbrock(gsl_multiroot_fsolver_dnewton,x_init_rosenbrock, p_rosenbrock);

		printf("Running Powell with the GSL newton solver!\n");
		run_powell(gsl_multiroot_fdfsolver_newton,x_init_powell,  p_powell);
	}
	else if (argc == 2 && strcmp(argv[1], "1") == 0)
	{
		printf("Running Rosenbrock with the custom newton solver!\n");
		run_rosenbrock(gsl_multiroot_fsolver_dnewton,x_init_rosenbrock, p_rosenbrock);

		printf("Running Powell with the custom newton solver!\n");
		run_powell(gsl_multiroot_fdfsolver_newton,x_init_powell,p_powell);
	}

}
