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
			maxes.q = std::fmax(maxes.q, std::abs(assem_sol.q_BC[i + 1]));
			maxes.z = std::fmax(maxes.z, std::abs(assem_sol.z_BC[i + 1]));
			maxes.eta = std::fmax(maxes.eta, std::abs(assem_sol.h_BC[i + 1] + assem_sol.z_BC[i + 1]));
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

			maxes.q = std::fmax(maxes.q, std::abs(q));
			maxes.eta = std::fmax(maxes.eta, std::abs(eta));
		}
	}

	maxes.q = std::fmax(maxes.q, C(1.0));
	maxes.eta = std::fmax(maxes.eta, C(1.0));
	maxes.z = std::fmax(maxes.z, C(1.0));

	return maxes;
}