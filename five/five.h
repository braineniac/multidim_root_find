#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

struct five_params
{
	double a;
	double b;
	double c;
	double d;
	double e;

} five_params;

int five_sq_f  (const gsl_vector * x,  void * params, gsl_vector * f);
int five_sq_df (const gsl_vector *x, void * p, gsl_matrix * J);
int five_sq_fdf (const gsl_vector *x, void *p, gsl_vector * f, gsl_matrix * J);
int run_five_sq_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p);
int run_five_sq_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p);

int five_trig_f  (const gsl_vector * x,  void * params, gsl_vector * f);
int five_trig_df (const gsl_vector *x, void * p, gsl_matrix * J);
int five_trig_fdf (const gsl_vector *x, void *p, gsl_vector * f, gsl_matrix * J);
int run_five_trig_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p);
int run_five_trig_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p);

int five_lin_f  (const gsl_vector * x,  void * params, gsl_vector * f);
int five_lin_df (const gsl_vector *x, void * p, gsl_matrix * J);
int five_lin_fdf (const gsl_vector *x, void *p, gsl_vector * f, gsl_matrix * J);
int run_five_lin_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p);
int run_five_lin_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p);

int five_exp_f  (const gsl_vector * x,  void * params, gsl_vector * f);
int five_exp_df (const gsl_vector *x, void * p, gsl_matrix * J);
int five_exp_fdf (const gsl_vector *x, void *p, gsl_vector * f, gsl_matrix * J);
int run_five_exp_f(const gsl_multiroot_fsolver_type *T, double * x_init, struct five_params p);
int run_five_exp_fdf(const gsl_multiroot_fdfsolver_type *T, double  * x_init, struct five_params p);
