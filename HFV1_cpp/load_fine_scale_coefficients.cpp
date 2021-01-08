#include "load_fine_scale_coefficients.h"

void load_fine_scale_coefficients
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	int&                  scale_coeffs_per_cell, 
	AssembledSolution&    assem_sol, 
	FlattenedScaleCoeffs& scale_coeffs
)
{
	for (int cell = 0; cell < sim_params.cells; cell++)
	{
		int scaleStep = cell * scale_coeffs_per_cell;

		for (int k = 0; k < (1 << solver_params.L); k++)
		{
			scale_coeffs.q[scaleStep + (1 << solver_params.L) - 1 + k] = assem_sol.q_BC[cell * (1 << solver_params.L) + k + 1];
			scale_coeffs.eta[scaleStep + (1 << solver_params.L) - 1 + k] = assem_sol.h_BC[cell * (1 << solver_params.L) + k + 1] + assem_sol.z_BC[cell * (1 << solver_params.L) + k + 1];
			scale_coeffs.z[scaleStep + (1 << solver_params.L) - 1 + k] = assem_sol.z_BC[cell * (1 << solver_params.L) + k + 1];
		}
	}
}