#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

typedef struct
  {
    gsl_matrix * lu;
    gsl_permutation * permutation;
  }
newton_custom_state_fdf_t;

static int newton_custom_fdf_alloc (void * vstate, size_t n);
static int newton_custom_fdf_set (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx);
static int newton_custom_fdf_iterate (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx);
static void newton_custom_fdf_free (void * vstate);

static const gsl_multiroot_fdfsolver_type newton_custom_fdf_type = 
{
	"newton_custom_fdf",
	sizeof(newton_custom_state_fdf_t),
	&newton_custom_fdf_alloc,
	&newton_custom_fdf_set,
	&newton_custom_fdf_iterate,
	&newton_custom_fdf_free
};

const gsl_multiroot_fdfsolver_type * gsl_multiroot_fdfsolver_newton_custom_fdf = &newton_custom_fdf_type;
