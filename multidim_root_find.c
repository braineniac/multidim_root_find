#include <stdio.h>
#include <string.h>
#include <argp.h>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include "powell/powell.c"
#include "rosenbrock/rosenbrock.c"

#include "tools.c"

#include "five/five.h"
#include "five/five_lin.c"
#include "five/five_sq.c"
#include "five/five_trig.c"
#include "five/five_exp.c"

#include "multidim_root_find.h"

#include "solvers/newton_custom_f.c"
#include "solvers/newton_custom_fdf.c"

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

int run_five_lin_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p)
{
	gsl_multiroot_fsolver *s;
	
	int status;
	size_t iter = 0;
	
	const size_t n = 5;
	gsl_multiroot_function f = {&five_lin_f, n, &p};
	
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

int run_five_sq_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p)
{
	gsl_multiroot_fsolver *s;
	
	int status;
	size_t iter = 0;
	
	const size_t n = 5;
	gsl_multiroot_function f = {&five_sq_f, n, &p};
	
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

int run_five_trig_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p)
{
	gsl_multiroot_fsolver *s;
	
	int status;
	size_t iter = 0;
	
	const size_t n = 5;
	gsl_multiroot_function f = {&five_lin_f, n, &p};
	
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

int run_five_exp_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p)
{
	gsl_multiroot_fsolver *s;
	
	int status;
	size_t iter = 0;
	
	const size_t n = 5;
	gsl_multiroot_function f = {&five_exp_f, n, &p};
	
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

int run_five_sq_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p)
{
	int status;
	size_t iter =0;

	gsl_multiroot_fdfsolver *s;

	const size_t n =5;

	gsl_multiroot_function_fdf f = {&five_sq_f,&five_sq_df,&five_sq_fdf,n,&p};

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

int run_five_lin_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p)
{
	int status;
	size_t iter =0;

	gsl_multiroot_fdfsolver *s;

	const size_t n =5;

	gsl_multiroot_function_fdf f = {&five_lin_f,&five_lin_df,&five_lin_fdf,n,&p};

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


int run_five_exp_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p)
{
	int status;
	size_t iter =0;

	gsl_multiroot_fdfsolver *s;

	const size_t n =5;

	gsl_multiroot_function_fdf f = {&five_exp_f,&five_exp_df,&five_exp_fdf,n,&p};

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

int run_five_trig_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p)
{
	int status;
	size_t iter =0;

	gsl_multiroot_fdfsolver *s;

	const size_t n =5;

	gsl_multiroot_function_fdf f = {&five_trig_f,&five_trig_df,&five_trig_fdf,n,&p};

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

static char doc[] = 
"This program finds the roots of non-linear systems of equations using GSL.\n";

static char args_doc[] = "Args: --solver(-s), --output(-o), --quiet(-q), --with-derivative(-d)";

static struct argp_option options[] = {
	{"solver",		's',	"SOLVER",	0,	"Solver to use"},
	{"output",		'o',	"FILE",		0,	"Output to a file"},
	{"quiet",		'q',	0,		0,	"No screen output"},
	{"with-derivative",	'd',	0,		0,	"Switch to fdf solver"},
	{0}
};

struct arguments {
	char *args[1];
	int verbose;
	char *solver;
	char *output_file;
	int derivative;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
	
	struct arguments *arguments = state->input;

	switch(key) {
		case 'q':
			arguments->verbose = 0;
			break;
		case 'o':
			arguments->output_file=arg;
			break;
		case 's':
			arguments->solver=arg;
			break;
		case 'd':
			arguments->derivative=1;
			break;
		case ARGP_KEY_ARG:
			if (state->arg_num >= 1)
				argp_usage(state);
			arguments->args[state->arg_num] = arg;
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 1)
				argp_usage(state);
			break;
		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}


static struct argp argp = { options, parse_opt, args_doc, doc};

int main(int argc, char *argv[]) {

	double x_init_rosenbrock[] = {-10.0,-5.0};
	double x_init_powell[] = {10.0,-20.0};
	double x_init_five_lin_f[] = {10,5.0,-4.0,-3.0,-8.0};
	double x_init_five_trig_f[] = {10,5.0,-4.0,-3.0,-8.0};
	double x_init_five_sq_f[] = {10,5.0,-4.0,-3.0,-8.0};
	double x_init_five_exp_f[] = {10,5.0,-4.0,-3.0,-8.0};

	
	struct rosenbrock_params p_rosenbrock = {.a=1.0, .b=10};
	struct powell_params p_powell = {.A=10000};
	struct five_params p_five_f = {.a=5.64, .b=-6.32, .c=3.14259, .d=2.78, .e=-8.523};
	
	struct arguments arguments;

	arguments.verbose=1;
	arguments.output_file="-";
	arguments.solver = "GSL_newton";
	arguments.derivative = 0;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	printf( "ARG=%s\n"
		"OUTPUT_FILE=%s\n"
		"VERBOSE=%s\n"
		"SOLVER=%s\n"
		"DERIVATIVE=%s\n",
		arguments.args[0],
		arguments.output_file,
		arguments.verbose ? "yes" : "no",
		arguments.solver,
		arguments.derivative ? "yes": "no");
	
	if (strcmp(arguments.args[0],"rosenbrock") == 0) {
		printf("Running Rosenbrock with the GSL newton solver!\n");
		run_rosenbrock(gsl_multiroot_fsolver_dnewton,x_init_rosenbrock, p_rosenbrock);
	}
	else if (strcmp(arguments.args[0],"powell") == 0) {
		printf("Running Powell with the GSL newton solver, with fdf!\n");
		run_powell(gsl_multiroot_fdfsolver_newton,x_init_powell,  p_powell);
	}
	
	else if (strcmp(arguments.args[0],"five_lin") == 0) {
		if (arguments.derivative) {
			printf("Running Five lin with the GSL newton solver, with derivative!\n");
			run_five_lin_fdf(gsl_multiroot_fdfsolver_newton,x_init_five_lin_f,p_five_f);
		}
		else {
			printf("Running Five lin with the GSL newton solver, no derivative!\n");
			run_five_lin_f(gsl_multiroot_fsolver_dnewton,x_init_five_lin_f,p_five_f);
		}
	}
	else if (strcmp(arguments.args[0],"five_sq") == 0) {
		if (arguments.derivative) {
			printf("Running Five sq with the GSL newton solver, with derivative!\n");
			run_five_sq_fdf(gsl_multiroot_fdfsolver_newton,x_init_five_sq_f,p_five_f);
		}
		else {
			printf("Running Five lin with the GSL newton solver, no derivative!\n");
			run_five_sq_f(gsl_multiroot_fsolver_dnewton,x_init_five_sq_f,p_five_f);
		}
	}
	else if (strcmp(arguments.args[0],"five_exp") == 0) {
		if (arguments.derivative) {
			printf("Running Five exp with the GSL newton solver, with derivative!\n");
			run_five_exp_fdf(gsl_multiroot_fdfsolver_newton,x_init_five_exp_f,p_five_f);
		}
		else {
			printf("Running Five exp with the GSL newton solver, no derivative!\n");
			run_five_exp_f(gsl_multiroot_fsolver_dnewton,x_init_five_exp_f,p_five_f);
		}
	}
	else if (strcmp(arguments.args[0],"five_trig") == 0) {
		if (arguments.derivative) {
			printf("Running Five trig with the GSL newton solver, with derivative!\n");
			run_five_trig_fdf(gsl_multiroot_fdfsolver_newton,x_init_five_trig_f,p_five_f);
		}
		else {
			printf("Running Five lin with the GSL newton solver, no derivative!\n");
			run_five_trig_f(gsl_multiroot_fsolver_dnewton,x_init_five_trig_f,p_five_f);
		}
	}
	else
		exit (1);
	
	exit (0);
}
