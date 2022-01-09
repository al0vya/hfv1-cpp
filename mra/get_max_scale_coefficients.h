#pragma once

#include <cmath>

#include "Maxes.h"
#include "AssembledSolution.h"
#include "FlattenedScaleCoeffs.h"

Maxes get_max_scale_coefficients
(
	Maxes&                maxes, 
	int&                  num_fine_cells, 
	AssembledSolution&    assem_sol, 
	bool&                 first_time_step, 
	FlattenedScaleCoeffs& scale_coeffs
);