#include "malloc_bar_values.h"

void malloc_bar_values
(
	BarValues& bar_vals,
	int&       num_fine_cells
)
{
	bar_vals.h = new real[num_fine_cells];
	bar_vals.z = new real[num_fine_cells];
}