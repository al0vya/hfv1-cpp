#include "malloc_assembled_solution.h"

void malloc_assembled_solution
(
	AssembledSolution&    assem_sol,
	int&                  num_fine_cells
)
{
	assem_sol.q_BC          = new real[num_fine_cells + 2];
	assem_sol.h_BC          = new real[num_fine_cells + 2];
	assem_sol.z_BC          = new real[num_fine_cells + 2];
	assem_sol.x             = new real[num_fine_cells + 2];
	assem_sol.dx_BC         = new real[num_fine_cells + 2];
	assem_sol.activeIndices = new int[num_fine_cells + 2];
}