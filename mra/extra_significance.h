#pragma once

#include <cmath>

#include "../classes/SimulationParameters.h"
#include "../classes/SolverParameters.h"

void extra_significance
(
	SimulationParameters& sim_params,
	int&                  details_per_cell, 
	SolverParameters&     solver_params,
	int*&                 sig_details, 
	real*&                norm_details
);