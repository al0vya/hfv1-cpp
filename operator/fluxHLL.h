#pragma once
#include <cmath>
#include <cmath>

using namespace std;

#include "real.h"
#include "AssembledSolution.h"
#include "SolverParameters.h"
#include "StarValues.h"
#include "Fluxes.h"

void fluxHLL
(
	AssembledSolution& assem_sol, 
	SolverParameters&  solverParameters,
	StarValues&        star_vals,
	Fluxes&            fluxes
);