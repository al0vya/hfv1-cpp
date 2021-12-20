#include "get_wet_dry_cells.h"

void get_wet_dry_cells
(
	int*               dry, 
	AssembledSolution& assem_sol, 
	SolverParameters&  solver_params
)
{
	dry[0] = false;
	dry[assem_sol.length + 1] = false;

	// initialising dry vs wet cells, ignore ghost cells
	for (int i = 1; i < assem_sol.length + 1; i++)
	{
		real h_loc = assem_sol.h_BC[i];
		real h_back = assem_sol.h_BC[i - 1];
		real h_forward = assem_sol.h_BC[i + 1];

		real h_max = std::max(h_loc, h_back);
		     h_max = std::max(h_forward, h_max);

		dry[i] = (h_max <= solver_params.tol_dry);
	}
}