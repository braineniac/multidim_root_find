#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

typedef struct
{
	gsl_matrix * J;
	gsl_matrix * lu;
	gsl_permutation * permutation;	

} newton_custom_state_f_t;

static int newton_custom_f_alloc (void * vstate, size_t n);
static int newton_custom_f_set (void * vstate, gsl_multiroot_function * function, gsl_vector * x, gsl_vector * f, gsl_vector * dx);
static int newton_custom_f_iterate (void * vstate, gsl_multiroot_function* function, gsl_vector * x, gsl_vector * f, gsl_vector * dx);
static void newton_custom_f_free (void * vstate);

static const gsl_multiroot_fsolver_type newton_custom_f_type =
{	"newton_custom_f",
	sizeof (newton_custom_state_f_t),
	&newton_custom_f_alloc,
	&newton_custom_f_set,
	&newton_custom_f_iterate,
	&newton_custom_f_free
};

const gsl_multiroot_fsolver_type * gsl_multiroot_fsolver_newton_custom_f = &newton_custom_f_type;

