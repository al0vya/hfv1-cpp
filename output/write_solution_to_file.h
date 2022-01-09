#pragma once

#include <fstream>
#include <string>

#include "SimulationParameters.h"
#include "AssembledSolution.h"
#include "SaveInterval.h"

void write_solution_to_file
(
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol,
	SaveInterval&         saveint
);