#pragma once

#include <algorithm>

#include "../classes/SimulationParameters.h"
#include "../classes/SolverParameters.h"
#include "../classes/FlattenedScaleCoeffs.h"
#include "../classes/FlattenedDetails.h"
#include "encode_detail.h"
#include "encode_scale.h"
#include "../classes/Maxes.h"

void encoding
(
	SimulationParameters& sim_params, 
	int&                  scale_coeffs_per_cell, 
	int&                  details_per_cell,
	int&                  num_details,
	SolverParameters&     solver_params, 
	int*&                 sig_details, 
	FlattenedScaleCoeffs& scale_coeffs, 
	FlattenedDetails&     details,
	real*&                norm_details,
	bool&                 first_time_step,
	Maxes&                maxes
);