#pragma once

#include "../classes/SimulationParameters.h"
#include "../classes/SolverParameters.h"

void regularisation
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	int&                  details_per_cell, 
	int*                  sig_details
);