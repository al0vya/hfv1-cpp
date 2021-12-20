#include "fv1_operator.h"

void fv1_operator
(
	int*&              dry, 
	Fluxes&            fluxes, 
	SolverParameters&  solver_params, 
	BarValues&         bar_vals, 
	AssembledSolution& assem_sol, 
	real&              dt
)
{
	for (int i = 1; i < assem_sol.length + 1; i++)
	{
		// skip increment in dry cells
		if (!dry[i])
		{
			real mass_increment = -(1 / assem_sol.dx_BC[i]) * (fluxes.mass[i] - fluxes.mass[i - 1]);
			real momentum_increment = -(1 / assem_sol.dx_BC[i]) * (fluxes.momentum[i] - fluxes.momentum[i - 1] + 2 * sqrt(C(3.0)) * solver_params.g * bar_vals.h[i - 1] * bar_vals.z[i - 1]);

			assem_sol.h_BC[i] += dt * mass_increment;
			assem_sol.q_BC[i] = (assem_sol.h_BC[i] <= solver_params.tol_dry) ? 0 : assem_sol.q_BC[i] + dt * momentum_increment;
		}
	}
}