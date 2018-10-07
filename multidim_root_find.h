#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include "tools.h"
#include "tools.c"

#include "powell/powell.h"
#include "powell/powell.c"

#include "rosenbrock/rosenbrock.h"
#include "rosenbrock/rosenbrock.c"

#include "five/five.h"
#include "five/five_lin.c"
#include "five/five_sq.c"
#include "five/five_trig.c"
#include "five/five_exp.c"

#include "solvers/newton_custom_f.h"
#include "solvers/newton_custom_f.c"

#include "solvers/newton_custom_fdf.h"
#include "solvers/newton_custom_fdf.c"
