#include "regularisation.h"

void regularisation
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	int&                  details_per_cell, 
	int*                  sig_details
)
{
	for (int cell = 0; cell < sim_params.cells; cell++)
	{
		int detailStep = cell * details_per_cell;

		for (int n = solver_params.L; n > 1; n--)
		{
			int k = (1 << (n - 1)) - 1;
			int kLower = (1 << (n - 2)) - 1;
			int currentLevEnd = (1 << n) - 2;

			for (k; k < currentLevEnd; k += 2)
			{
				if (sig_details[detailStep + k] || sig_details[detailStep + k + 1])
				{
					sig_details[detailStep + kLower] = true;
				}

				kLower++; // step along only the one parent element
			}
		}
	}
}