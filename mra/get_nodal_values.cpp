#include "get_nodal_values.h"

void get_nodal_values
(
	NodalValues&          nodal_vals,
	SimulationParameters& sim_params,
	int&                  num_fine_cells,
	BoundaryConditions&   bcs,
	real&                 dx,
	int                   test_case
)
{
	for (int i = 0; i < num_fine_cells + 1; i++) nodal_vals.x[i] = sim_params.xmin + i * dx;
	
	switch (test_case)
	{
	case 1:
	case 2:
	case 3:
		for (int i = 0; i < num_fine_cells + 1; i++)
		{
			nodal_vals.z[i] = C(0.0);
			nodal_vals.h[i] = h_init_overtop(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = C(0.0);
		}
		break;
	case 4:
	case 5:
		for (int i = 0; i < num_fine_cells + 1; i++)
		{
			nodal_vals.z[i] = bed_data_c_property(nodal_vals.x[i]);
			nodal_vals.h[i] = h_init_c_property(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = C(0.0);
		}
		break;
	case 6:
		for (int i = 0; i < num_fine_cells + 1; i++)
		{
			nodal_vals.z[i] = bed_data_c_property(nodal_vals.x[i]);
			nodal_vals.h[i] = h_init_overtop(bcs, nodal_vals.z[i], nodal_vals.x[i]);
			nodal_vals.q[i] = C(0.0);
		}
		break;
	case 7:
		for (int i = 0; i < num_fine_cells + 1; i++)
		{
			nodal_vals.z[i] = bed_data_triangle(nodal_vals.x[i]);
			nodal_vals.h[i] = ( nodal_vals.x[i] < C(15.5) ) ? C(0.75) : C(0.0);;
			nodal_vals.q[i] = C(0.0);
		}
		break;
	default:
		break;
	}
}