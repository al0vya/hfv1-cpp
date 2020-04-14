#pragma once
#include "real.h"
#include <cmath>

typedef struct SolverParameters
{
	real CFL;
	real tolDry;
	real g;
	real epsilon;
	int L;

} SolverParameters;