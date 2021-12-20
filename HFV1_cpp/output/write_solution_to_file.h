#pragma once

#include <fstream>

#include "../classes/SimulationParameters.h"
#include "../classes/AssembledSolution.h"

void write_solution_to_file
(
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol
);