#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

int run_five_trig_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p) {
	
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

	do {
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


int run_five_trig_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p) {
	
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

int five_trig_f (const gsl_vector * x, void * params, gsl_vector * f) {

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

	const double y0 = c*cos(x1) + b*x0 + d*x3 + e*x4;
	const double y1 = b*x1 - d*x4 + e*x3;
	const double y2 = a*x0 - b*x3 + x2;
	const double y3 = b*d*x4 + c*x2;
	const double y4 = a*e*x0 + b*x3;

	gsl_vector_set(f,0,y0);
	gsl_vector_set(f,1,y1);
	gsl_vector_set(f,2,y2);
	gsl_vector_set(f,3,y3);
	gsl_vector_set(f,4,y4);

	return GSL_SUCCESS;

}

int five_trig_df (const gsl_vector *x, void *params, gsl_matrix * J) {

	//const double x0 = gsl_vector_get(x,0);
	const double x1 = gsl_vector_get(x,1);
	//const double x2 = gsl_vector_get(x,2);
	//const double x3 = gsl_vector_get(x,3);
	//const double x4 = gsl_vector_get(x,4);

	double a = ((struct five_params *) params)->a;
	double b = ((struct five_params *) params)->b;
	double c = ((struct five_params *) params)->c;
	double d = ((struct five_params *) params)->d;
	double e = ((struct five_params *) params)->e;
	
	gsl_matrix_set(J,0,0, b);
	gsl_matrix_set(J,0,1, -c*sin(x1));
	gsl_matrix_set(J,0,2, 0);
	gsl_matrix_set(J,0,3, d);
	gsl_matrix_set(J,0,4, e);

	gsl_matrix_set(J,1,0, 0);
	gsl_matrix_set(J,1,1, b);
	gsl_matrix_set(J,1,2, 0);
	gsl_matrix_set(J,1,3, e);
	gsl_matrix_set(J,1,4, -d);

	gsl_matrix_set(J,2,0, a);
	gsl_matrix_set(J,2,1, 0);
	gsl_matrix_set(J,2,2, 1);
	gsl_matrix_set(J,2,3, -b);
	gsl_matrix_set(J,2,4, 0);
	
	gsl_matrix_set(J,3,0, 0);
	gsl_matrix_set(J,3,1, 0);
	gsl_matrix_set(J,3,2, c);
	gsl_matrix_set(J,3,3, 0);
	gsl_matrix_set(J,3,4, b*d);

	gsl_matrix_set(J,4,0, a*e);
	gsl_matrix_set(J,4,1, 0);
	gsl_matrix_set(J,4,2, 0);
	gsl_matrix_set(J,4,3, b);
	gsl_matrix_set(J,4,4, 0);

	return GSL_SUCCESS;
}

int five_trig_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {

	five_trig_f(x, params, f);
	five_trig_df(x, params, J);

	return GSL_SUCCESS;

}
