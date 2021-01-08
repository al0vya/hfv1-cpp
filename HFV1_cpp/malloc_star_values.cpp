#include "malloc_star_values.h"

void malloc_star_values
(
	StarValues& star_values,
	int&        num_fine_cells
)
{
	star_values.q_east = new real[num_fine_cells + 1];
	star_values.h_east = new real[num_fine_cells + 1];

	star_values.q_west = new real[num_fine_cells + 1];
	star_values.h_west = new real[num_fine_cells + 1];

	star_values.z_east = new real[num_fine_cells + 1];
	star_values.z_west = new real[num_fine_cells + 1];
}