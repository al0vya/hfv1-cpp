#include "get_dt_CFL.h"

void get_dt_CFL
(
	SolverParameters&  solver_params, 
	AssembledSolution& assem_sol, 
	real&              dt, 
	real&              total_mass
)
{
	for (int i = 1; i < assem_sol.length + 1; i++)
	{
		if (assem_sol.h_BC[i] > solver_params.tol_dry)
		{
			real u = assem_sol.q_BC[i] / assem_sol.h_BC[i];
			real dtCFL = solver_params.CFL * assem_sol.dx_BC[i] / (abs(u) + sqrt(solver_params.g * assem_sol.h_BC[i]));
			dt = std::min(dt, dtCFL);
		}

		total_mass += assem_sol.h_BC[i] * assem_sol.dx_BC[i];
	}
}