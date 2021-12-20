#include "thresholding.h"

void thresholding
(
	SimulationParameters& sim_params, 
	SolverParameters&     solver_params, 
	int&                  scale_coeffs_per_cell, 
	int&                  details_per_cell, 
	real*&                norm_details, 
	int*&                 sig_details
)
{
	for (int cell = 0; cell < sim_params.cells; cell++)
	{
		int scaleStep = cell * scale_coeffs_per_cell;
		int detailStep = cell * details_per_cell;

		for (int n = 0; n < solver_params.L; n++)
		{
			int currentLevStart = (1 << n) - 1;
			int currentLevEnd = (1 << (n + 1)) - 2;

			for (int k = currentLevStart; k <= currentLevEnd; k++)
			{
				real epsilonLocal = solver_params.epsilon / (1 << (solver_params.L - n));

				if (norm_details[detailStep + k] >= epsilonLocal)
				{
					sig_details[detailStep + k] = true;

					if (cell + 1 < sim_params.cells && k == currentLevEnd) // if it's not the last cell and at the right-most edge
					{
						sig_details[(cell + 1) * details_per_cell + currentLevStart] = true; // make the subelement to the right cell also significant
					}

					if (cell > 0 && k == currentLevStart) // if it's not the first cell and at the left-most edge
					{
						sig_details[(cell - 1) * details_per_cell + currentLevEnd] = true; // the subelement on the left cell is also significant
					}
				}
			}
		}
	}
}