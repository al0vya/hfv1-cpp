#pragma once

#include <cmath>

#include "AssembledSolution.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"

void get_dt_CFL
(
	SolverParameters&  solver_params, 
	AssembledSolution& assem_sol, 
	real&              dt, 
	real&              total_mass
);