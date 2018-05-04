#include <stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include "tools.h"

int print_state_f(size_t iter, gsl_multiroot_fsolver * s)
{
	printf("iter = %zu  x = % .3f % .3f "
			 "f(s) = % .3e % .3e\n",
			 iter,
			 gsl_vector_get(s->x,0),
			 gsl_vector_get(s->x,1),
			 gsl_vector_get(s->f,0),
			 gsl_vector_get(s->f,1));
	return 0;
}

int print_state_fdf(size_t iter, gsl_multiroot_fdfsolver *s)
{
	printf("iter = %zu  x = % .3f % .3f "
			"f(s) = % .3e % .3e\n",
			iter,
			gsl_vector_get(s->x,0),
			gsl_vector_get(s->x,1),
			gsl_vector_get(s->f,0),
			gsl_vector_get(s->f,1));
	return 0;
}

int print_state_five(size_t iter, gsl_multiroot_fsolver *s)
{
	 printf("iter = %zu  x = % .3f % .3f % .3f % .3f % .3f "
			 "f(s) = % .3e % .3e % .3f % .3f % .3f\n",
			 iter,
			 gsl_vector_get(s->x,0),
			 gsl_vector_get(s->x,1),
			 gsl_vector_get(s->x,2),
			 gsl_vector_get(s->x,3),
			 gsl_vector_get(s->x,4),
			 gsl_vector_get(s->f,0),
			 gsl_vector_get(s->f,1),
			 gsl_vector_get(s->f,2),
			 gsl_vector_get(s->f,3),
			 gsl_vector_get(s->f,4));
	 return 0;
}
