#include "extra_significance.h"

void extra_significance
(
	SimulationParameters& sim_params,
	int&                  details_per_cell, 
	SolverParameters&     solver_params,
	int*&                 sig_details, 
	real*&                norm_details
)
{
	for (int cell = 0; cell < sim_params.cells; cell++)
	{
		int detailStep = cell * details_per_cell;

		for (int n = 0; n < solver_params.L; n++)
		{
			int currentLevStart = (1 << n) - 1;
			int currentLevEnd = (1 << (n + 1)) - 2;
			int kHigher = currentLevEnd + 1; // index of child element

			for (int k = currentLevStart; k <= currentLevEnd; k++)
			{
				if (sig_details[detailStep + k])
				{
					real mBar = 1.5;
					real epsilonLocal = solver_params.epsilon * pow(C(2.0), n - solver_params.L);

					if (norm_details[detailStep + k] >= epsilonLocal * pow(C(2.0), mBar + 1) && n + 1 != solver_params.L)
					{
						// if extra signficant child elements marked as active
						sig_details[detailStep + kHigher] = true;
						sig_details[detailStep + kHigher + 1] = true;
					}
				}

				kHigher += 2;
			}
		}
	}
}