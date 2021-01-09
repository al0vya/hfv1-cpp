#include "encoding.h"

void encoding
(
	SimulationParameters& sim_params, 
	int&                  scale_coeffs_per_cell, 
	int&                  details_per_cell,
	int&                  num_details,
	SolverParameters&     solver_params, 
	int*&                 sig_details, 
	FlattenedScaleCoeffs& scale_coeffs, 
	FlattenedDetails&     details,
	real*&                norm_details,
	bool&                 first_time_step,
	Maxes&                maxes

)
{
	// thresholding details to zero for next step
	for (int i = 0; i < num_details; i++)
	{
		details.q[i] = 0;
		details.eta[i] = 0;
	}
	
	for (int cell = 0; cell < sim_params.cells; cell++)
	{
		int scaleStep = cell * scale_coeffs_per_cell;
		int detailStep = cell * details_per_cell;

		for (int n = solver_params.L - 1; n >= 0; n--)
		{
			int currentLevStart = (1 << n) - 1;
			int currentLevEnd = (1 << (n + 1)) - 2;
			int kHigher = currentLevEnd + 1;

			for (int k = currentLevStart; k <= currentLevEnd; k++)
			{
				if (sig_details[detailStep + k] || first_time_step)
				{
					real q1 = scale_coeffs.q[scaleStep + kHigher];
					real eta1 = scale_coeffs.eta[scaleStep + kHigher];

					real q2 = scale_coeffs.q[scaleStep + kHigher + 1];
					real eta2 = scale_coeffs.eta[scaleStep + kHigher + 1];

					scale_coeffs.q[scaleStep + k] = encodeScale(q1, q2);
					scale_coeffs.eta[scaleStep + k] = encodeScale(eta1, eta2);

					details.q[detailStep + k] = encodeDetail(q1, q2);
					details.eta[detailStep + k] = encodeDetail(eta1, eta2);

					if (first_time_step)
					{
						real z1 = scale_coeffs.z[scaleStep + kHigher];
						real z2 = scale_coeffs.z[scaleStep + kHigher + 1];

						scale_coeffs.z[scaleStep + k] = encodeScale(z1, z2);
						details.z[detailStep + k] = encodeDetail(z1, z2);
					}
				}

				kHigher += 2;
			}
		}
	}

	first_time_step = false;

	for (int i = 0; i < num_details; i++)
	{
		sig_details[i] = false;

		real a = abs(details.q[i]) / maxes.q;
		real b = abs(details.eta[i]) / maxes.eta;
		real c = abs(details.z[i]) / maxes.z;

		norm_details[i] = std::max(a, std::max(b, c));
	}
}