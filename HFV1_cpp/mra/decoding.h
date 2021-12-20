#pragma once

#include "AssembledSolution.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "FlattenedDetails.h"
#include "traverse_tree_decode.h"
#include "FlattenedScaleCoeffs.h"
#include "FlattenedDetails.h"

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