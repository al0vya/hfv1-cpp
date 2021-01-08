#include "malloc_nodal_values.h"

void malloc_nodal_values
(
	NodalValues&          nodal_values,
	int                   num_fine_cells
)
{
	nodal_values.q = new real[num_fine_cells + 1];
	nodal_values.h = new real[num_fine_cells + 1];
	nodal_values.z = new real[num_fine_cells + 1];
	nodal_values.x = new real[num_fine_cells + 1];
}