#pragma once

#include <algorithm>

#include "../classes/Maxes.h"
#include "../classes/AssembledSolution.h"
#include "../classes/FlattenedScaleCoeffs.h"

Maxes get_max_scale_coefficients
(
	Maxes&                maxes, 
	int&                  num_fine_cells, 
	AssembledSolution&    assem_sol, 
	bool&                 first_time_step, 
	FlattenedScaleCoeffs& scale_coeffs
);