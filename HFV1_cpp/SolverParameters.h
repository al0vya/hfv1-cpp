#pragma once
#include "real.h"
#include <cmath>

typedef struct SolverParameters
{
	int cells;
	real CFL;
	real tolDry;
	real g;
	real epsilon;
	int L;

} SolverParameters;