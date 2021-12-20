#include "write_solution_to_file.h"

void write_solution_to_file
(
	SimulationParameters& sim_params, 
	AssembledSolution&    assem_sol
)
{
	real x1 = C(0.0);
	real x2 = C(0.0);

	std::ofstream test;

	test.open("solution_data.csv");

	test << "x1,x2,q,z,eta" << std::endl;

	for (int i = 0; i < assem_sol.length; i++)
	{
		x2 += assem_sol.dx_BC[i + 1];

		test << x1 << "," 
		     << x2 << "," 
			 << assem_sol.q_BC[i + 1] << "," 
			 << assem_sol.z_BC[i + 1] << "," 
			 << std::max(assem_sol.z_BC[i + 1], assem_sol.h_BC[i + 1] + assem_sol.z_BC[i + 1]) 
			 << std::endl;

		x1 += assem_sol.dx_BC[i + 1];
	}

	test.close();
}