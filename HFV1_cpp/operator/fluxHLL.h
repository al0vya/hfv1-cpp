#pragma once
#include <cmath>
#include <algorithm>

using namespace std;

#include "../classes/real.h"
#include "../classes/AssembledSolution.h"
#include "../classes/SolverParameters.h"
#include "../classes/StarValues.h"
#include "../classes/Fluxes.h"

void fluxHLL
(
	AssembledSolution& assem_sol, 
	SolverParameters&  solverParameters,
	StarValues&        star_vals,
	Fluxes&            fluxes
);