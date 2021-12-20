#pragma once
#include "real.h"

typedef struct SolverParameters
{
	real CFL;
	real tol_dry;
	real g;
	real epsilon;
	int L;

} SolverParameters;