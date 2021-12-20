#pragma once

#include <fstream>

#include "SimulationParameters.h"
#include "AssembledSolution.h"

void write_solution_to_file
(
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol
);