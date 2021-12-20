#pragma once

#include <algorithm>

#include "../classes/AssembledSolution.h"
#include "../classes/SimulationParameters.h"
#include "../classes/SolverParameters.h"

void get_dt_CFL
(
	SolverParameters&  solver_params, 
	AssembledSolution& assem_sol, 
	real&              dt, 
	real&              total_mass
);