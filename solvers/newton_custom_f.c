#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include "newton_custom_f.h"

static int newton_custom_f_alloc (void * vstate, size_t n)
{
	newton_custom_state_f_t * state = (newton_custom_state_f_t *) vstate;
	gsl_permutation * p;
	gsl_matrix * m, *J;

	m = gsl_matrix_calloc(n,n);

	if (m==0)
		GSL_ERROR ("Failed to allocate space for an n x n matrix!\n", GSL_ENOMEM);
	
	state->lu = m;

	p = gsl_permutation_calloc (n);
	if (p == 0)
	{
		gsl_matrix_free(m);
		GSL_ERROR ("failed to allocate space for permutation", GSL_ENOMEM);
	}
	
	state->permutation = p ;

	J = gsl_matrix_calloc (n,n);
	
	if (J == 0)
	{
		gsl_permutation_free(p);
		gsl_matrix_free(m);
		
		GSL_ERROR ("failed to allocate space for d", GSL_ENOMEM);
	}
	
	state->J = J;

	return GSL_SUCCESS;
}

static int newton_custom_f_set (void * vstate,gsl_multiroot_function * function, gsl_vector * x, gsl_vector * f, gsl_vector * dx)
{
	newton_custom_state_f_t * state = (newton_custom_state_f_t *) vstate;

	size_t i, n = function->n ;
	int status;

	status = GSL_MULTIROOT_FN_EVAL (function, x, f);

	if (status)
		return status;

	status = gsl_multiroot_fdjacobian (function, x, f, GSL_SQRT_DBL_EPSILON,state->J);

	if (status)
		return status;
       
	for (i = 0; i < n; i++)
    	{
		gsl_vector_set (dx, i, 0.0);
	}
	
	return GSL_SUCCESS;
}

static int newton_custom_f_iterate (void * vstate, gsl_multiroot_function * function, gsl_vector * x, gsl_vector * f, gsl_vector * dx)
{
	newton_custom_state_f_t * state = (newton_custom_state_f_t *) vstate;
	
	int signum;
	size_t i;
	size_t n = function->n ;

	gsl_matrix_memcpy (state->lu, state->J);
	
    	int status = gsl_linalg_LU_decomp (state->lu, state->permutation, &signum);

    	if (status)
      		return status;   
      
	status = gsl_linalg_LU_solve (state->lu, state->permutation, f, dx);
	if (status)
		return status;


  	for (i = 0; i < n; i++)
  	{
	 	double e = gsl_vector_get (dx, i);
	  	double y = gsl_vector_get (x, i);
	  	gsl_vector_set (dx, i, -e);
	  	gsl_vector_set (x, i, y - e);
  	}

   	status = GSL_MULTIROOT_FN_EVAL (function, x, f);

    	if (status != GSL_SUCCESS)
		return GSL_EBADFUNC;
    	
	gsl_multiroot_fdjacobian (function, x, f, GSL_SQRT_DBL_EPSILON, state->J);


    	return GSL_SUCCESS;
}


static void newton_custom_f_free (void * vstate)
{
  	newton_custom_state_f_t * state = (newton_custom_state_f_t *) vstate;
	
	gsl_matrix_free(state->J);

  	gsl_matrix_free(state->lu);

  	gsl_permutation_free(state->permutation);
}
/*
static const gsl_multiroot_fsolver_type newton_custom_type =
{       "newton_custom",
	sizeof (newton_custom_state_f_t),
	&newton_custom_f_alloc,
	&newton_custom_f_set,
	&newton_custom_f_iterate,
	&newton_custom_f_free
};

const gsl_multiroot_fsolver_type * gsl_multiroot_fsolver_newton_custom = &newton_custom_type;
*/
