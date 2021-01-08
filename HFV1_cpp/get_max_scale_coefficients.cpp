#include "get_max_scale_coefficients.h"

Maxes get_max_scale_coefficients
(
	Maxes&                maxes, 
	int&                  num_fine_cells, 
	AssembledSolution&    assem_sol, 
	bool&                 first_time_step, 
	FlattenedScaleCoeffs& scale_coeffs
)
{
	maxes.q = 0;
	maxes.eta = 0;

	if (first_time_step)
	{
		for (int i = 0; i < num_fine_cells; i++)
		{
			maxes.q = std::max(maxes.q, abs(assem_sol.q_BC[i + 1]));
			maxes.z = std::max(maxes.z, abs(assem_sol.z_BC[i + 1]));
			maxes.eta = std::max(maxes.eta, abs(assem_sol.h_BC[i + 1] + assem_sol.z_BC[i + 1]));
		}
	}
	else
	{
		for (int i = 0; i < assem_sol.length; i++)
		{
			real q = assem_sol.q_BC[i + 1];
			real eta = assem_sol.h_BC[i + 1] + assem_sol.z_BC[i + 1];

			scale_coeffs.q[assem_sol.activeIndices[i]] = q;
			scale_coeffs.eta[assem_sol.activeIndices[i]] = eta;

			maxes.q = std::max(maxes.q, abs(q));
			maxes.eta = std::max(maxes.eta, abs(eta));
		}
	}

	maxes.q = std::max(maxes.q, C(1.0));
	maxes.eta = std::max(maxes.eta, C(1.0));
	maxes.z = std::max(maxes.z, C(1.0));

	return maxes;
}