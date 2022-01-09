#pragma once

#include "SimulationParameters.h"
#include "NodalValues.h"
#include "BoundaryConditions.h"
#include "../operator/h_init_c_property.h"
#include "../operator/h_init_overtop.h"
#include "../operator/bed_data_triangle.h"
#include "../operator/bed_data_c_property.h"

void get_nodal_values
(
	NodalValues&          nodal_vals,
	SimulationParameters& sim_params,
	int&                  num_fine_cells,
	BoundaryConditions&   bcs,
	real&                 dx,
	int                   test_case
);