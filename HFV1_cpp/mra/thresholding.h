#pragma once

#include "SimulationParameters.h"
#include "SolverParameters.h"


void thresholding
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	int&                  scale_coeffs_per_cell, 
	int&                  details_per_cell, 
	real*&                norm_details, 
	int*&                 sig_details
);