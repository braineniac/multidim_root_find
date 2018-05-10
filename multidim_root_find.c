#include <stdio.h>
#include <string.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include "powell.c"
#include "rosenbrock.c"
#include "tools.c"
#include "five.c"
#include "multidim_root_find.h"

#include "newton_custom_f.c"
#include "newton_custom_fdf.c"

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

int run_five_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p)
{
	gsl_multiroot_fsolver *s;
	
	int status;
	size_t iter = 0;
	
	const size_t n = 5;
	gsl_multiroot_function f = {&five_f, n, &p};
	
	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set(x, 0, x_init[0]);
	gsl_vector_set(x, 1, x_init[1]);
	gsl_vector_set(x, 2, x_init[2]);
	gsl_vector_set(x, 3, x_init[3]);
	gsl_vector_set(x, 4, x_init[4]);
	
	s = gsl_multiroot_fsolver_alloc (T,5);
	gsl_multiroot_fsolver_set(s, &f,x);

	print_state_five_f(iter,s);

	do
	{
		 iter++;
		 status = gsl_multiroot_fsolver_iterate (s);

		 print_state_five_f(iter,s);

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

int run_five_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p)
{
	int status;
	size_t iter =0;

	gsl_multiroot_fdfsolver *s;

	const size_t n =5;

	gsl_multiroot_function_fdf f = {&five_f,&five_df,&five_fdf,n,&p};

	gsl_vector *x = gsl_vector_alloc(n);
	
	gsl_vector_set(x,0,x_init[0]);
	gsl_vector_set(x,1,x_init[1]);
	gsl_vector_set(x,2,x_init[2]);
	gsl_vector_set(x,3,x_init[3]);
	gsl_vector_set(x,4,x_init[4]);
	
	s = gsl_multiroot_fdfsolver_alloc(T,5);
	gsl_multiroot_fdfsolver_set(s,&f,x);

	print_state_five_fdf(iter,s);

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
	double x_init_five_f[5] = {10.0,5.0,-4.0,-3.0,-8.0};
	struct rosenbrock_params p_rosenbrock = {.a=1.0, .b=10};
	struct powell_params p_powell = {.A=10000};
	struct five_params p_five_f = {.a=5.64, .b=-6.32, .c=3.14259, .d=2.78, .e=-8.523};
	
	if (argc==2 && strcmp(argv[1],"0") == 0 )
	{
		printf("Running Rosenbrock with the GSL newton solver!\n");
		run_rosenbrock(gsl_multiroot_fsolver_dnewton,x_init_rosenbrock, p_rosenbrock);

		printf("Running Powell with the GSL newton solver!\n");
		run_powell(gsl_multiroot_fdfsolver_newton,x_init_powell,  p_powell);

		printf("Running Five with the GSL newton solver, no derivative!\n");
		run_five_f(gsl_multiroot_fsolver_dnewton,x_init_five_f, p_five_f);

		printf("Running Five with the GSL newton solver, with derivative!\n");
		run_five_fdf(gsl_multiroot_fdfsolver_newton,x_init_five_f,p_five_f);
	}
	else if (argc == 2 && strcmp(argv[1], "1") == 0)
	{
		printf("Running Rosenbrock with the custom newton solver!\n");
		run_rosenbrock(gsl_multiroot_fsolver_newton_custom_f,x_init_rosenbrock, p_rosenbrock);

		printf("Running Powell with the custom newton solver!\n");
		run_powell(gsl_multiroot_fdfsolver_newton_custom_fdf,x_init_powell,p_powell);
	
		printf("Running Five with the custom newton solver!\n");
		run_five_f(gsl_multiroot_fsolver_newton_custom_f,x_init_five_f, p_five_f);
	
	
		printf("Running Five with the custom newton solver, with derivative!\n");
		run_five_fdf(gsl_multiroot_fdfsolver_newton_custom_fdf,x_init_five_f,p_five_f);	
	}

}
