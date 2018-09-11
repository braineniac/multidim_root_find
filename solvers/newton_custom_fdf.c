#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include "newton_custom_fdf.h"
	
static int newton_custom_fdf_alloc (void * vstate, size_t n)
{
	newton_custom_state_fdf_t * state = (newton_custom_state_fdf_t *) vstate;
	gsl_permutation * p;
	gsl_matrix * m;
	
	m = gsl_matrix_calloc (n,n);
	
	if (m == 0)
		GSL_ERROR ("failed to allocate space for lu", GSL_ENOMEM);
	
	state->lu = m ;
	
	p = gsl_permutation_calloc (n);
	
	if (p == 0)
	{
		gsl_matrix_free(m);
		
		GSL_ERROR ("failed to allocate space for permutation", GSL_ENOMEM);
	}
	
	state->permutation = p ;
	
	return GSL_SUCCESS;
}

static int newton_custom_fdf_set (void * vstate, gsl_multiroot_function_fdf * FDF, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx)
{
	newton_custom_state_fdf_t * state = (newton_custom_state_fdf_t *) vstate;
	
	size_t i, n = FDF->n ;
	
	state = 0 ; /* avoid warnings about unused parameters */
	
	GSL_MULTIROOT_FN_EVAL_F_DF (FDF, x, f, J);
	
	for (i = 0; i < n; i++)
	{
		gsl_vector_set (dx, i, 0.0);
	}
	
	return GSL_SUCCESS;
}

static int newton_custom_fdf_iterate (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx)
{
	newton_custom_state_fdf_t * state = (newton_custom_state_fdf_t *) vstate;
  	
	int signum;
	size_t i;
	
	size_t n = fdf->n ;
	
	gsl_matrix_memcpy (state->lu, J);
	
	gsl_linalg_LU_decomp (state->lu, state->permutation, &signum);
	
	
	int status = gsl_linalg_LU_solve (state->lu, state->permutation, f, dx);
	
	if (status)
		return status;
	
	for (i = 0; i < n; i++)
	{
		double e = gsl_vector_get (dx, i);
		double y = gsl_vector_get (x, i);
		gsl_vector_set (dx, i, -e);
		gsl_vector_set (x, i, y - e);
	}
	
	status = GSL_MULTIROOT_FN_EVAL_F_DF (fdf, x, f, J);
	
	if (status != GSL_SUCCESS) 
	{
		return GSL_EBADFUNC;
	}
	
	return GSL_SUCCESS;
}


static void newton_custom_fdf_free (void * vstate)
{
	newton_custom_state_fdf_t * state = (newton_custom_state_fdf_t *) vstate;
	
	gsl_matrix_free(state->lu);

	gsl_permutation_free(state->permutation);
}

/*
static const gsl_multiroot_fdfsolver_type newton_type =
{"newton",                             
 sizeof (newton_custom_state_fdf_t),
 &newton_alloc,
 &newton_set,
 &newton_iterate,
 &newton_free};

const gsl_multiroot_fdfsolver_type * gsl_multiroot_fdfsolver_newton = &newton_type;
*/
