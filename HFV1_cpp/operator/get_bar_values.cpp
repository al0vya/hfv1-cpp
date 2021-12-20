#include "get_bar_values.h"

void get_bar_values
(
	StarValues&        star_vals, 
	AssembledSolution& assem_sol, 
	BarValues&         bar_vals
)
{
	for (int i = 0; i < assem_sol.length; i++)
	{
		// essentially 0th order projection but taking into account east/west locality
		bar_vals.h[i] = (star_vals.h_west[i + 1] + star_vals.h_east[i]) / 2;

		// 1st order projection
		bar_vals.z[i] = (star_vals.z_west[i + 1] - star_vals.z_east[i]) / (2 * sqrt(C(3.0)));
	}
}