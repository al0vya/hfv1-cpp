#pragma once

#include "../classes/AssembledSolution.h"
#include "../classes/SimulationParameters.h"
#include "../classes/SolverParameters.h"
#include "../classes/FlattenedDetails.h"
#include "traverse_tree_decode.h"
#include "../classes/FlattenedScaleCoeffs.h"
#include "../classes/FlattenedDetails.h"

void decoding
(
	AssembledSolution&    assem_sol, 
	SimulationParameters& sim_params, 
	int&                  scale_coeffs_per_cell, 
	int&                  details_per_cell,
	real*&                xFlattened, 
	real*&                xIntCoarse, 
	real*&                dxFlattened, 
	real&                 dx_coarse, 
	int*&                 levelIndicesFlattened, 
	int*&                 sig_details, 
	SolverParameters&     solver_params, 
	FlattenedScaleCoeffs  scale_coeffs, 
	FlattenedDetails      details
);