#pragma once
#include "real.h"

typedef struct SolverParameters
{
	real CFL;
	real tolDry;
	real g;
	real epsilon;
	int L;

} SolverParameters;