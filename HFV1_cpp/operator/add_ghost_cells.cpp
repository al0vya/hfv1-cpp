#include "add_ghost_cells.h"

void add_ghost_cells
(
	AssembledSolution&    assem_sol, 
	BoundaryConditions&   bcs
)
{
	assem_sol.q_BC[0] = assem_sol.q_BC[1];
	assem_sol.q_BC[assem_sol.length + 1] = assem_sol.q_BC[assem_sol.length];

	assem_sol.h_BC[0] = assem_sol.h_BC[1];
	assem_sol.h_BC[assem_sol.length + 1] = assem_sol.h_BC[assem_sol.length];

	assem_sol.z_BC[0] = assem_sol.z_BC[1];
	assem_sol.z_BC[assem_sol.length + 1] = assem_sol.z_BC[assem_sol.length];
}