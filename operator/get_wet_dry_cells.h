#pragma once

#include <algorithm>

#include "SimulationParameters.h"
#include "AssembledSolution.h"
#include "SolverParameters.h"

void get_wet_dry_cells
(
	int*               dry, 
	AssembledSolution& assem_sol, 
	SolverParameters&  solver_params
);