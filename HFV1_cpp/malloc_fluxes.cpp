#include "malloc_fluxes.h"

void malloc_fluxes
(
	Fluxes& fluxes,
	int&    num_fine_cells
)
{
	fluxes.mass     = new real[num_fine_cells + 1];
	fluxes.momentum = new real[num_fine_cells + 1];
}