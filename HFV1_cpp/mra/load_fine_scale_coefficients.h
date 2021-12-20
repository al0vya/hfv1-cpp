#pragma once

#include "../classes/SolverParameters.h"
#include "../classes/SimulationParameters.h"
#include "../classes/FlattenedScaleCoeffs.h"
#include "../classes/AssembledSolution.h"

void load_fine_scale_coefficients
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	int&                  scale_coeffs_per_cell, 
	AssembledSolution&    assem_sol, 
	FlattenedScaleCoeffs& scale_coeffs
);