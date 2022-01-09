#pragma once

#include <fstream>
#include <string>

#include "../classes/SimulationParameters.h"
#include "../classes/AssembledSolution.h"
#include "../classes/SaveInterval.h"

void write_solution_to_file
(
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol,
	SaveInterval&         saveint
);