#pragma once

#include <cmath>

#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "Fluxes.h"
#include "AssembledSolution.h"
#include "BarValues.h"

void fv1_operator
(
	int*&              dry, 
	Fluxes&            fluxes, 
	SolverParameters&  solver_params, 
	BarValues&         bar_vals, 
	AssembledSolution& assem_sol, 
	real&              dt
);